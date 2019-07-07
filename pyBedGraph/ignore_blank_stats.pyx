import numpy as np
cimport cython
from libc.float cimport DBL_MAX, DBL_MIN
from libc.math cimport ceil
from libc.stdlib cimport malloc, free

def load_smallest_bins(double[:] value_list, unsigned int bin_size):
    cdef size_t value_list_size = value_list.size
    cdef size_t bin_index

    cdef unsigned int numb_bins = <int>ceil(value_list_size / bin_size)

    bins = np.full(numb_bins, -1, dtype=np.float64)
    bins_coverage = np.zeros(numb_bins, dtype=np.float64)
    cdef double[:] bins_view = bins
    cdef double[:] bins_coverage_view = bins_coverage

    cdef unsigned int start, end
    cdef double mean, coverage

    for bin_index in range(numb_bins):
        start = bin_index * bin_size
        end = start + bin_size

        bins_view[bin_index] = get_bin_mean(value_list, start, end)
        bins_coverage_view[bin_index] = get_coverage(value_list, start, end)

        '''if mean != -1:
            bins_view[bin_index] = mean
            bins_coverage_view[bin_index] = coverage'''

    return bins, bins_coverage

def load_bins(double[:] prev_bin_level_mean, double[:] prev_bin_level_coverage):

    assert tuple(prev_bin_level_mean.shape) == tuple(prev_bin_level_coverage.shape)

    cdef size_t prev_bin_level_size = prev_bin_level_mean.size
    cdef size_t bin_index

    # just take the average of two bins from prev level
    cdef char bin_size = 2

    cdef unsigned int numb_bins = <int>ceil(prev_bin_level_size / bin_size)

    bins = np.full(numb_bins, -1, dtype=np.float64)
    bins_coverage = np.zeros(numb_bins, dtype=np.float64)
    cdef double[:] bins_view = bins
    cdef double[:] bins_coverage_view = bins_coverage

    cdef unsigned int prev_bin_index
    cdef double mean1, mean2, cov1, cov2, avg_cov, avg_mean, cov_sum

    for bin_index in range(numb_bins):
        prev_bin_index = bin_index * bin_size

        # a slightly more complicated average of the previous bin level
        mean1 = prev_bin_level_mean[prev_bin_index]
        if prev_bin_index + 1 < prev_bin_level_size:
            mean2 = prev_bin_level_mean[prev_bin_index + 1]
            cov2 = prev_bin_level_coverage[prev_bin_index + 1]
        else:
            mean2 = -1
            cov2 = 0

        cov1 = prev_bin_level_coverage[prev_bin_index]
        cov_sum = cov1 + cov2
        avg_cov = cov_sum / 2

        # both mean1 and mean2 are -1
        # this bin does not cover anything
        if avg_cov == 0:
            continue

        if mean1 != -1 and mean2 != -1:
            avg_mean = (mean1 * cov1 + mean2 * cov2) / cov_sum
        elif mean2 == -1:
            avg_mean = mean1
        else:
            avg_mean = mean2

        bins_view[bin_index] = avg_mean
        bins_coverage_view[bin_index] = avg_cov

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
    cdef unsigned int num_counted = 0
    cdef size_t i

    for i in range(start, end, 1):
        value = value_list[i]
        if value < 0:
            continue

        total += value
        num_counted += 1

    return total, num_counted

def get_exact_means(double[:] value_list, double[:] bin_list, int bin_size,
                         double[:] bin_coverage_list, int[:] start_list,
                         int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, start, end, bin_end, bin_index
    cdef size_t num_tests = start_list.size
    cdef double total, value, coverage
    cdef unsigned int numb_value, weight, prev_length

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
            value = bin_list[bin_index]
            if value != -1:
                coverage = bin_coverage_list[bin_index]
                weight = <int>(coverage * bin_size)
                total += weight * value
                numb_value += weight

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

def get_approx_means(double[:] bin_list, int max_bin_size,
        int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, start, end, bin_end, bin_index
    cdef size_t num_tests = start_list.size
    cdef double total, value
    cdef unsigned int numb_value, weight

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]

        bin_index = <unsigned int>(start / max_bin_size)
        bin_end = <unsigned int>(end / max_bin_size)

        # special case where interval is within a single bin
        if bin_index == bin_end:
            if bin_list[bin_index] == -1:
                continue
            result_view[i] = bin_list[bin_index]

        total = 0
        numb_value = 0

        # first bin
        value = bin_list[bin_index]
        if value != -1:
            weight = (bin_index + 1) * max_bin_size - start
            total += weight * value
            numb_value += weight

        # middle bins
        weight = max_bin_size
        bin_index += 1
        while bin_index < bin_end:
            value = bin_list[bin_index]
            if value != -1:
                total += weight * value
                numb_value += weight

            bin_index += 1

        # last bin
        value = bin_list[bin_index]
        if value != -1:
            weight = end - bin_index * max_bin_size
            total += weight * value
            numb_value += weight

        if numb_value == 0:
            continue

        result_view[i] = total / numb_value

    return result

def get_mod_approx_means(list bins_list, unsigned int max_bin_size,
                         unsigned int bin_list_numb, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef double[:] max_size_bin_list = bins_list[bin_list_numb - 1]
    cdef double[:] bin_list = bins_list[0]
    cdef double* bin_array

    cdef size_t i, num_tests = start_list.size, start, end
    cdef size_t max_bin_index
    cdef double total
    cdef unsigned int numb_value, weight, bin_end, bin_index, current_bin_size
    cdef unsigned int bin_size_index, arr_size = bin_list.size

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

    cdef double ** bins_array = <double **>malloc(bin_list_numb * sizeof(double*))
    for i in range(bin_list_numb):
        bins_array[i] = <double *>malloc(arr_size * sizeof(double))
        bin_list = bins_list[i]
        for j in range(arr_size):
            bins_array[i][j] = bin_list[j]
        arr_size /= 2

    for i in range(num_tests):
        start = start_list[i]
        end = end_list[i]

        max_bin_index = <unsigned int>(start / max_bin_size)
        bin_end = <unsigned int>(end / max_bin_size)

        # special case where interval is within a single bin
        # TODO: Look at smaller bins if possible?
        if max_bin_index == bin_end:
            if max_size_bin_list[max_bin_index] == -1:
                continue
            result_view[i] = max_size_bin_list[max_bin_index]

        total = 0
        numb_value = 0

        # first bin
        if max_size_bin_list[max_bin_index] != -1:
            weight = (max_bin_index + 1) * max_bin_size - start

            current_bin_size = max_bin_size
            bin_size_index = bin_list_numb - 1
            bin_index = max_bin_index
            while weight < current_bin_size / 2 and bin_size_index > 0:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2 + 1

            bin_array = bins_array[bin_size_index]
            if bin_array[bin_index] != -1:
                total += weight * bin_array[bin_index]
                numb_value += weight

        # middle bins
        weight = max_bin_size
        max_bin_index += 1
        while max_bin_index < bin_end:
            if max_size_bin_list[max_bin_index] != -1:
                total += weight * max_size_bin_list[max_bin_index]
                numb_value += weight

            max_bin_index += 1

        # last bin
        if max_size_bin_list[max_bin_index] != -1:
            weight = end - max_bin_index * max_bin_size

            current_bin_size = max_bin_size
            bin_size_index = bin_list_numb - 1
            bin_index = max_bin_index
            while weight < current_bin_size / 2 and bin_size_index > 0:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2

            bin_array = bins_array[bin_size_index]
            if bin_array[bin_index] != -1:
                total += weight * bin_array[bin_index]
                numb_value += weight

        if numb_value == 0:
            continue

        result_view[i] = total / numb_value

    for i in range(bin_list_numb):
        free(bins_array[i])
    free(bins_array)

    return result

def get_minimums(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef double minimum

    result = np.full(num_tests, -1, dtype=np.float64)
    cdef double[:] result_view = result

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

cdef get_coverage(double[:] value_list, size_t start, size_t end):
    cdef unsigned int num_covered = 0
    cdef size_t i
    cdef double fraction_covered
    cdef size_t value_list_size = value_list.size
    cdef size_t orig_end = end

    if end > value_list_size:
        end = value_list_size

    for i in range(start, end, 1):
        if value_list[i] < 0:
            continue

        num_covered += 1

    fraction_covered = num_covered / (orig_end - start)
    return fraction_covered
