import numpy as np
from libc.float cimport DBL_MAX, DBL_MIN
from libc.math cimport ceil, sqrt

def load_smallest_bins(double[:] value_list, unsigned int bin_size):
    cdef size_t value_list_size = value_list.size
    cdef size_t bin_index

    cdef unsigned int numb_bins = <int>ceil(value_list_size / bin_size)

    bins = np.full(numb_bins, -1, dtype=np.float64)
    bins_coverage = np.zeros(numb_bins, dtype=np.uint32)
    cdef double[:] bins_view = bins
    cdef unsigned int[:] bins_coverage_view = bins_coverage

    cdef unsigned int start, end, coverage
    cdef double value

    for bin_index in range(numb_bins):
        start = bin_index * bin_size
        end = start + bin_size

        if end > value_list_size:
            end = value_list_size

        value, coverage = get_values(value_list, start, end)
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

cdef get_bin_mean(double[:] value_list, size_t start, size_t end):
    cdef double total = 0
    cdef unsigned int num_counted = 0
    cdef size_t i
    cdef size_t value_list_size = value_list.size

    if end > value_list_size:
        end = value_list_size

    for i in range(start, end, 1):
        if value_list[i] < 0:
            continue

        num_counted += 1
        total += value_list[i]

    if num_counted == 0:
        return -1

    return total / num_counted

cdef get_values(double[:] value_list, start, end):

    cdef double total = 0, value
    cdef unsigned int coverage = 0
    cdef size_t i

    for i in range(start, end, 1):
        value = value_list[i]
        if value < 0:
            continue

        total += value
        coverage += 1

    return total, coverage

cpdef get_exact_means(double[:] value_list, double[:] bin_list, int bin_size,
                         unsigned int[:] bin_coverage_list, int[:] start_list,
                         int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, start, end, bin_end, bin_index
    cdef size_t num_tests = start_list.size
    cdef double total, value
    cdef unsigned int numb_value, weight, prev_length, coverage

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]

        bin_index = <unsigned int>ceil(start / bin_size)
        bin_end = <unsigned int>(end / bin_size)

        # special case where interval is within a single bin
        if bin_index >= bin_end:
            value = get_bin_mean(value_list, start, end)
            if value != -1:
                result_view[i] = value
                continue

        total = 0
        numb_value = 0

        # first bin
        if start % bin_size != 0:
            prev_length = bin_size - start % bin_size
            value, weight = get_values(value_list, start, start + prev_length)
            if value != -1:
                total += value
                numb_value += weight

        # middle bins
        while bin_index < bin_end:
            coverage = bin_coverage_list[bin_index]
            if coverage > 0:
                total += bin_list[bin_index]
                numb_value += coverage

            bin_index += 1

        # last bin
        prev_length = end % bin_size
        if prev_length != 0:
            value, weight = get_values(value_list, end - prev_length, end)
            if value != -1:
                total += value
                numb_value += weight

        if numb_value == 0:
            continue

        result_view[i] = total / numb_value

    return result

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
            fraction = (max_bin_size - start % max_bin_size) / max_bin_size
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
            fraction = (end % max_bin_size) / max_bin_size
            total += bin_list[bin_index] * fraction
            numb_value += weight * fraction

        if numb_value == 0:
            continue

        result_view[i] = total / numb_value

    return result


def get_minimums(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef double minimum

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result, wanted_list, interval

    for i in range(num_tests):
        minimum = DBL_MAX
        start = start_list[i]
        end = end_list[i]
        for j in range(start, end, 1):
            if value_list[j] < 0:
                continue

            if value_list[j] < minimum:
                minimum = value_list[j]

        if minimum != DBL_MAX:
            result_view[i] = minimum

    return result

def get_maximums(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef double maximum

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        maximum = DBL_MIN
        start = start_list[i]
        end = end_list[i]
        for j in range(start, end, 1):
            if value_list[j] < 0:
                continue

            if value_list[j] > maximum:
                maximum = value_list[j]

        if maximum != DBL_MIN:
            result_view[i] = maximum

    return result

cpdef get_coverages(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef unsigned int numb_covered

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        numb_covered = 0
        start = start_list[i]
        end = end_list[i]
        for j in range(start, end, 1):
            # if 0 was given as value, it is "covered" in the bedGraph
            if value_list[j] >= 0:
                numb_covered += 1

        result_view[i] = numb_covered / (end - start)

    return result

def get_stds(double[:] value_list, double[:] bin_list, unsigned int bin_size,
            unsigned int[:] bin_coverage_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, j

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result
    cdef double std, difference, value
    cdef unsigned int total

    cdef double[:] means = get_exact_means(value_list, bin_list, bin_size,
                            bin_coverage_list, start_list, end_list)

    for i in range(num_tests):
        mean = means[i]
        if mean == -1:
            continue

        start = start_list[i]
        end = end_list[i]

        std = 0
        total = 0
        for j in range(start, end, 1):
            value = value_list[j]

            if value == -1:
                continue

            difference = value - mean
            std += difference * difference
            total += 1

        std /= total
        result_view[i] = sqrt(std)

    return result
