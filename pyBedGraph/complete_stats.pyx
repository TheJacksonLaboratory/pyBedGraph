import numpy as np
cimport cython
from libc.float cimport DBL_MAX, DBL_MIN

#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)   # Deactivate negative indexing.

# contiguous array
#def mean(double[::1] values, int[::1] start_list, int[::1] end_list):
def get_exact_mean(double[:] value_list, int[:] start_list, int[:] end_list):

    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t i, num_tests = start_list.size, start, end
    cdef double total

    result = np.zeros(num_tests, dtype=np.float64)
    cdef double[:] result_view = result

    for i in range(num_tests):
        total = 0
        start = start_list[i]
        end = end_list[i]
        for j in range(start, end, 1):
            total += value_list[j]

        result_view[i] = total / (end - start)

    return result

def get_min(double[:] value_list, int[:] start_list, int[:] end_list):

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

def get_min(double[:] value_list, int[:] start_list, int[:] end_list):

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

def get_max(double[:] value_list, int[:] start_list, int[:] end_list):

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

def get_coverage(double[:] value_list, int[:] start_list, int[:] end_list):

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
