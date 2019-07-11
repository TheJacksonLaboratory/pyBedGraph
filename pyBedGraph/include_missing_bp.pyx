import numpy as np
cimport cython
from libc.float cimport DBL_MAX, DBL_MIN
from libc.math cimport ceil, sqrt

#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)   # Deactivate negative indexing.

def load_smallest_bins(double[:] value_list, unsigned int bin_size):
    cdef size_t value_list_size = value_list.size
    cdef size_t bin_index

    cdef unsigned int numb_bins = <int>ceil(value_list_size / bin_size)

    bins = np.zeros(numb_bins, dtype=np.float64)
    cdef double[:] bins_view = bins

    cdef unsigned int start, end
    cdef double value

    for bin_index in range(numb_bins):
        start = bin_index * bin_size
        end = start + bin_size

        if end > value_list_size:
            end = value_list_size

        value = get_total(value_list, start, end)
        if value > 0:
            bins_view[bin_index] = value

    return bins

def load_bins(double[:] prev_bin_level_mean):

    cdef size_t prev_bin_level_size = prev_bin_level_mean.size
    cdef size_t bin_index

    # just take the average of two bins from prev level
    cdef char bin_size = 2

    cdef unsigned int numb_bins = <int>ceil(prev_bin_level_size / bin_size)

    bins = np.zeros(numb_bins, dtype=np.float64)
    cdef double[:] bins_view = bins

    cdef unsigned int prev_bin_index
    cdef double value

    for bin_index in range(numb_bins):
        prev_bin_index = bin_index * bin_size

        # just add them up
        value = prev_bin_level_mean[prev_bin_index]
        if prev_bin_index + 1 < prev_bin_level_size:
            value += prev_bin_level_mean[prev_bin_index + 1]

        if value > 0:
            bins_view[bin_index] = value

    return bins

cdef get_total(double[:] value_list, start, end):

    cdef double total = 0
    cdef size_t i

    for i in range(start, end, 1):
        total += value_list[i]

    return total

# contiguous array
#def mean(double[::1] values, int[::1] start_list, int[::1] end_list):
cpdef get_exact_means(double[:] value_list, double[:] bin_list, int bin_size,
                       int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, start, end, bin_end, bin_index
    cdef size_t num_tests = start_list.size
    cdef double total, value
    cdef unsigned int numb_value, weight, prev_length

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]
        numb_value = end - start

        bin_index = <unsigned int>ceil(start / bin_size)
        bin_end = <unsigned int>(end / bin_size)

        # special case where interval is within a single bin
        if bin_index >= bin_end:
            value = get_total(value_list, start, end)
            if value > 0:
                result_view[i] = value / (end - start)
                continue

        total = 0

        # first bin
        if start % bin_size != 0:
            prev_length = bin_size - start % bin_size
            value = get_total(value_list, start, start + prev_length)
            if value > 0:
                total += value

        # middle bins
        while bin_index < bin_end:
            value = bin_list[bin_index]
            if value > 0:
                total += value

            bin_index += 1

        # last bin
        prev_length = end % bin_size
        if prev_length != 0:
            value = get_total(value_list, end - prev_length, end)
            if value > 0:
                total += value

        if total == 0:
            continue

        result_view[i] = total / numb_value

    return result

def get_approx_means(double[:] bin_list, int max_bin_size, int[:] start_list,
                    int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, start, end, bin_end, bin_index
    cdef size_t num_tests = start_list.size
    cdef double total, numb_value, fraction, value

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]

        bin_index = <unsigned int>(start / max_bin_size)
        bin_end = <unsigned int>(end / max_bin_size)

        # special case where interval is within a single bin
        if bin_index == bin_end:
            if bin_list[bin_index] == 0:
                continue
            result_view[i] = bin_list[bin_index] / max_bin_size
            continue

        total = 0
        numb_value = 0

        # first bin
        value = bin_list[bin_index]
        if value > 0:
            fraction = (max_bin_size - start % max_bin_size) / max_bin_size
            total += bin_list[bin_index] * fraction
            numb_value += max_bin_size * fraction
        bin_index += 1

        # middle bins
        while bin_index < bin_end:
            value = bin_list[bin_index]
            if value > 0:
                total += bin_list[bin_index]
                numb_value += max_bin_size

            bin_index += 1

        # last bin
        value = bin_list[bin_index]
        if value > 0:
            fraction = (end % max_bin_size) / max_bin_size
            total += bin_list[bin_index] * fraction
            numb_value += max_bin_size * fraction

        if numb_value == 0:
            continue

        result_view[i] = total / numb_value

    return result

def get_medians(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t num_tests = start_list.size, i, start, end
    result = np.zeros(num_tests, dtype=np.float64)

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]

        result[i] = np.median(value_list[start:end])

    return result

def get_minimums(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef double minimum

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        minimum = DBL_MAX
        start = start_list[i]
        end = end_list[i]
        for j in range(start, end, 1):
            if value_list[j] < minimum:
                minimum = value_list[j]

        result_view[i] = minimum

    return result

def get_maximums(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef double maximum

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        maximum = DBL_MIN
        start = start_list[i]
        end = end_list[i]
        for j in range(start, end, 1):
            if value_list[j] > maximum:
                maximum = value_list[j]

        result_view[i] = maximum

    return result

def get_coverages(double[:] value_list, int[:] start_list, int[:] end_list):

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
            if value_list[j] > 0:
                numb_covered += 1

        result_view[i] = numb_covered / (end - start)

    return result

def get_medians(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result
    cdef double median

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]
        median = np.median(value_list[start:end])

        if median > 0:
            result_view[i] = median

    return result

def get_stds(double[:] value_list, double[:] bin_list, int bin_size,
                int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end, j

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result
    cdef double std, difference

    cdef double[:] means = get_exact_means(value_list, bin_list, bin_size,
                            start_list, end_list)

    for i in range(num_tests):
        mean = means[i]
        if mean == 0:
            continue

        start = start_list[i]
        end = end_list[i]

        std = 0
        for j in range(start, end, 1):
            difference = value_list[j] - mean
            std += difference * difference

        std /= (end - start)
        result_view[i] = sqrt(std)

    return result
