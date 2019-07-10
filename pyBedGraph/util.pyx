def fill_value_array(unsigned int[:] start_list, unsigned int[:] end_list,
                    double[:] value_map, double[:] value_list):
    assert tuple(start_list.shape) == tuple(end_list.shape)

    cdef size_t size = start_list.size, i, start, end

    for i in range(size):
        start = start_list[i]
        end = end_list[i]
        value_list[start:end] = value_map[i]
