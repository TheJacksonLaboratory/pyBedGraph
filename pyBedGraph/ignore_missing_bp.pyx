import numpy as np
import cython
from libc.float cimport DBL_MAX
from libc.math cimport ceil, sqrt

def load_smallest_bins(double[:] value_map, int[:] index_list, unsigned int size,
                       unsigned int[:] interval_start, unsigned int[:] interval_end,
                       unsigned int bin_size):
    cdef size_t bin_index

    cdef unsigned int numb_bins = <int>ceil(size / bin_size)

    bins = np.full(numb_bins, -1, dtype=np.float64)
    bins_coverage = np.zeros(numb_bins, dtype=np.uint32)
    cdef double[:] bins_view = bins
    cdef unsigned int[:] bins_coverage_view = bins_coverage

    cdef unsigned int start, end, coverage
    cdef double value

    for bin_index in range(numb_bins):
        start = bin_index * bin_size
        end = start + bin_size

        if end > size:
            end = size

        value, coverage = get_values(value_map, index_list, interval_start,
                                     interval_end, start, end)
        if coverage > 0:
            bins_view[bin_index] = value
            bins_coverage_view[bin_index] = coverage

    return bins, bins_coverage

def load_bins(double[:] prev_bin_level_mean, unsigned int[:] prev_bin_level_coverage):

    assert tuple(prev_bin_level_mean.shape) == tuple(prev_bin_level_coverage.shape)

    cdef size_t prev_bin_level_size = prev_bin_level_mean.size
    cdef size_t bin_index

    # just take the average of two bins from prev level
    cdef char bin_size = 2

    cdef unsigned int numb_bins = <int>ceil(prev_bin_level_size / bin_size)

    bins = np.full(numb_bins, -1, dtype=np.float64)
    bins_coverage = np.zeros(numb_bins, dtype=np.uint32)
    cdef double[:] bins_view = bins
    cdef unsigned int[:] bins_coverage_view = bins_coverage

    cdef unsigned int prev_bin_index, coverage
    cdef double value

    for bin_index in range(numb_bins):
        prev_bin_index = bin_index * bin_size

        # just add them up
        value = prev_bin_level_mean[prev_bin_index]
        coverage = prev_bin_level_coverage[prev_bin_index]
        if prev_bin_index + 1 < prev_bin_level_size:
            value += prev_bin_level_mean[prev_bin_index + 1]
            coverage += prev_bin_level_coverage[prev_bin_index + 1]

        if coverage > 0:
            bins_view[bin_index] = value
            bins_coverage_view[bin_index] = coverage

    return bins, bins_coverage

cdef get_values(double[:] value_map, int[:] index_list,
                unsigned int[:] interval_start, unsigned int[:] interval_end,
                unsigned int start, unsigned int end):

    cdef double total = 0, value
    cdef unsigned int coverage = 0, value_index, temp_end, interval_size
    cdef size_t i, numb_intervals = interval_start.size

    # get to an interval
    while start < end and index_list[start] == -1:
        start += 1

    if start == end:
        return total, coverage

    value_index = index_list[start]
    while start < end and start < interval_end[value_index]:
        temp_end = interval_end[value_index]
        if temp_end > end:
            temp_end = end
        interval_size = temp_end - start

        total += value_map[value_index] * interval_size
        coverage += interval_size

        value_index += 1
        if value_index == numb_intervals:
            break
        start = interval_start[value_index]

    return total, coverage

@cython.cdivision(True)
@cython.boundscheck(False)
def get_approx_means(double[:] bin_list, unsigned int[:] bin_coverage_list,
                     int max_bin_size, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, start, end, bin_end, bin_index
    cdef size_t num_tests = start_list.size
    cdef double total, numb_value, fraction
    cdef unsigned int weight

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]

        bin_index = <unsigned int>(start / max_bin_size)
        bin_end = <unsigned int>(end / max_bin_size)

        # special case where interval is within a single bin
        if bin_index == bin_end:
            if bin_coverage_list[bin_index] == 0:
                continue
            result_view[i] = bin_list[bin_index] / bin_coverage_list[bin_index]
            continue

        total = 0
        numb_value = 0

        # first bin
        weight = bin_coverage_list[bin_index]
        if weight > 0:
            fraction = <double>(max_bin_size - start % max_bin_size) / max_bin_size
            total += bin_list[bin_index] * fraction
            numb_value += weight * fraction
        bin_index += 1

        # middle bins
        while bin_index < bin_end:
            weight = bin_coverage_list[bin_index]
            if weight > 0:
                total += bin_list[bin_index]
                numb_value += weight

            bin_index += 1

        # last bin
        weight = bin_coverage_list[bin_index]
        if weight > 0:
            fraction = <double>(end % max_bin_size) / max_bin_size
            total += bin_list[bin_index] * fraction
            numb_value += weight * fraction

        if numb_value == 0:
            continue

        result_view[i] = total / numb_value

    return result

@cython.cdivision(True)
@cython.boundscheck(False)
cpdef get_exact_means(double[:] value_map, int[:] index_list,
                 unsigned int[:] interval_start, unsigned int[:] interval_end,
                 int[:] start_list, int[:] end_list):
    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, value_index
    cdef unsigned int numb_value, temp_end, interval_size, numb_intervals = interval_start.size
    cdef double total

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        total = 0
        numb_value = 0
        start = start_list[i]
        end = end_list[i]

        # get to an interval
        while index_list[start] == -1 and start < end:
            start += 1

        if start == end:
            continue

        value_index = index_list[start]
        while start < end:
            temp_end = interval_end[value_index]
            if temp_end > end:
                temp_end = end
            interval_size = temp_end - start

            total += value_map[value_index] * interval_size
            numb_value += interval_size

            value_index += 1
            if value_index == numb_intervals:
                break
            start = interval_start[value_index]

        if numb_value != 0:
            result_view[i] = total / numb_value

    return result

def get_minimums(double[:] value_map, int[:] index_list,
                 unsigned int[:] interval_start, unsigned int[:] interval_end,
                 int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, value_index
    cdef size_t numb_intervals = interval_start.size
    cdef double minimum

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        minimum = DBL_MAX
        start = start_list[i]
        end = end_list[i]

        # get to an interval
        while index_list[start] == -1 and start < end:
            start += 1

        if start == end:
            continue

        value_index = index_list[start]
        while interval_start[value_index] < end:
            if value_map[value_index] < minimum:
                minimum = value_map[value_index]

            value_index += 1
            if value_index == numb_intervals:
                break

        if minimum != DBL_MAX:
            result_view[i] = minimum

    return result

def get_maximums(double[:] value_map, int[:] index_list,
                 unsigned int[:] interval_start, unsigned int[:] interval_end,
                 int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, value_index
    cdef size_t numb_intervals = interval_start.size
    cdef double maximum, value

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        maximum = -1
        start = start_list[i]
        end = end_list[i]

        # get to an interval
        while index_list[start] == -1 and start < end:
            start += 1

        if start == end:
            continue

        value_index = index_list[start]
        while interval_start[value_index] < end:
            value = value_map[value_index]
            if value > maximum:
                maximum = value

            value_index += 1
            if value_index == numb_intervals:
                break

        if maximum != -1:
            result_view[i] = maximum

    return result

def get_coverages(int[:] index_list, unsigned int[:] interval_start,
                  unsigned int[:] interval_end, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, current_start
    cdef size_t numb_intervals = interval_start.size
    cdef unsigned int numb_covered, temp_end, value_index

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        numb_covered = 0
        start = start_list[i]
        end = end_list[i]
        current_start = start

        # get to an interval
        while index_list[current_start] == -1 and current_start < end:
            current_start += 1

        if current_start == end:
            continue

        value_index = index_list[current_start]
        while current_start < end:
            temp_end = interval_end[value_index]
            if temp_end > end:
                temp_end = end
            numb_covered += (temp_end - current_start)

            value_index += 1
            if value_index == numb_intervals:
                break
            current_start = interval_start[value_index]

        if numb_covered > 0:
            result_view[i] = <double>numb_covered / (end - start)

    return result

def get_stds(double[:] value_map, int[:] index_list,
                 unsigned int[:] interval_start, unsigned int[:] interval_end,
                 int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, j

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result
    cdef double std, difference, value
    cdef unsigned int numb_value, interval_size, value_index, temp_end
    cdef size_t numb_intervals = interval_start.size

    cdef double[:] means = get_exact_means(value_map, index_list, interval_start,
                            interval_end, start_list, end_list)

    for i in range(num_tests):
        mean = means[i]
        if mean == -1:
            continue

        start = start_list[i]
        end = end_list[i]

        # get to an interval
        while index_list[start] == -1 and start < end:
            start += 1

        if start == end:
            continue

        std = 0
        numb_value = 0

        value_index = index_list[start]
        while start < end:
            temp_end = interval_end[value_index]
            if temp_end > end:
                temp_end = end

            interval_size = temp_end - start
            difference = value_map[value_index] - mean
            std += difference * difference * interval_size
            numb_value += interval_size

            value_index += 1
            if value_index == numb_intervals:
                break
            start = interval_start[value_index]

        if numb_value == 0:
            continue

        std /= (numb_value)
        result_view[i] = sqrt(std)

    return result
