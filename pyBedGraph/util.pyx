def fill_value_array(unsigned int[:] start_list, unsigned int[:] end_list,
                    double[:] value_map, double[:] value_list):
    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t size = start_list.size, i, start, end

    for i in range(size):
        start = start_list[i]
        end = end_list[i]
        value_list[start:end] = value_map[i]

def fill_index_array(unsigned int[:] start_list, unsigned int[:] end_list,
                    double[:] value_map, int[:] index_list):
    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t size = start_list.size, i, start, end

    for i in range(size):
        start = start_list[i]
        end = end_list[i]
        index_list[start:end]= i

cpdef get_bin_value(double[:] value_map, int[:] index_list,
                 unsigned int[:] interval_start, unsigned int[:] interval_end,
                 int start, int end):

    cdef size_t i, value_index
    cdef unsigned int numb_value, temp_end, interval_size, numb_intervals = interval_start.size
    cdef double total

    total = 0
    numb_value = 0

    # get to an interval
    while start < index_list.size and index_list[start] == -1 and start < end:
        start += 1

    if start == end or start == index_list.size:
        return -1

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

    return total
