import math
from .ignore_missing_bp import *
from .util import *

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

MAX_NUMB_BIN_LIST = 1
MIN_BIN_SIZE = 2


class Chrom_Data:

    def __init__(self, name, size):
        self.name = name
        self.size = size

        # don't use this until user loads this chromosome for searching
        self.loaded_chrom = False
        self.index_list = None
        self.value_list = None
        self.total_size = 0

        # starting length is the size of chromosome
        # later shorten to save memory
        self.value_map = np.zeros(size, dtype=np.float64)

        # TODO:
        # make intervals unsigned shorts to use 75% memory when loading file
        self.intervals = [
            np.zeros(size, dtype=np.uint32),
            np.zeros(size, dtype=np.uint32)
        ]
        self.current_index = 0  # number of intervals in bedGraph file
        self.max_index = 0  # end of last interval in bedGraph file

        self.loaded_bins = False
        self.max_bin_size = None
        self.min_bin_size = None
        self.bins_list = []
        self.bins_list_coverages = []
        self.bin_list_numb = 0

        print(f"Reading in {name} ...")

    def add_data(self, data):
        start = int(data[START_INDEX])
        end = int(data[END_INDEX])
        value = float(data[VALUE_INDEX])

        self.intervals[0][self.current_index] = start
        self.intervals[1][self.current_index] = end
        self.total_size += (end - start)

        self.value_map[self.current_index] = value
        self.current_index += 1

        self.max_index = end

    # shorten length to # of intervals in bedGraph file for the chromosome
    def trim_extra_space(self):
        self.value_map = self.value_map[:self.current_index]
        self.intervals[0] = self.intervals[0][:self.current_index]
        self.intervals[1] = self.intervals[1][:self.current_index]

    # assume that missing space in bedGraph file is different from value of 0
    def initialize_index_array(self):
        # self.value_list = np.full(self.size, -1, dtype=np.float64)
        self.index_list = np.full(self.size, -1, dtype=np.int32)

    # fill the value array for fast indexing
    def load_index_array(self):

        print(f"Loading {self.name} ...")

        if self.loaded_chrom:
            print(f"{self.name} is already loaded")
            return

        self.initialize_index_array()
        if self.value_list is not None:
            fill_value_array(self.intervals[0], self.intervals[1], self.value_map,
                             self.value_list)

        if self.index_list is not None:
            fill_index_array(self.intervals[0], self.intervals[1], self.value_map,
                             self.index_list)

        self.loaded_chrom = True

    def free_value_array(self):
        self.value_list = None
        self.index_list = None
        self.loaded_chrom = False
        print(f"Freed memory for {self.name}'s value_list and index_list")

        self.bins_list.clear()
        if self.bins_list_coverages is not None:
            self.bins_list_coverages.clear()
        self.max_bin_size = None
        self.min_bin_size = None
        self.bin_list_numb = 0
        self.loaded_bins = False
        print(f"Freed memory for {self.name}'s bins")

    # load only set of bins for now
    def load_bins(self, max_bin_size):

        if max_bin_size is None:
            print("Did not specify max_bin_size")
            return

        if max_bin_size == self.max_bin_size and self.loaded_bins is True:
            print(f"Already loaded bins for: {max_bin_size}")
            return

        self.max_bin_size = max_bin_size

        self.bins_list.clear()
        self.bins_list_coverages.clear()

        # could make more sets of bins by changing constant
        bin_size = max_bin_size
        for i in range(MAX_NUMB_BIN_LIST - 1):
            if bin_size % 2 == 0:
                bin_size /= 2
            else:
                break
        bin_size = int(bin_size)
        self.min_bin_size = bin_size

        # Loading smallest bins
        print(f"Loading bins of size: {bin_size} for {self.name} ...")
        print(f"Number of bins: {math.ceil(self.size / bin_size)}")
        prev_bins_list, prev_bins_coverage_list = \
            load_smallest_bins(self.value_map, self.index_list, self.size,
                               self.intervals[0], self.intervals[1], bin_size)

        self.bins_list.append(prev_bins_list)
        self.bins_list_coverages.append(prev_bins_coverage_list)
        bin_size *= 2

        # Load larger bins
        while bin_size <= max_bin_size:

            print(f"Loading bins of size: {bin_size} for {self.name} ...")
            print(f"Number of bins: {math.ceil(self.size / bin_size)}")

            # use the previous bin list to speed up the process
            bins_list, bins_coverage_list = load_bins(prev_bins_list,
                                                      prev_bins_coverage_list)
            prev_bins_list = bins_list
            prev_bins_coverage_list = bins_coverage_list

            self.bins_list.append(bins_list)
            self.bins_list_coverages.append(bins_coverage_list)

            bin_size *= 2

        self.bin_list_numb = len(self.bins_list)
        self.loaded_bins = True

        # test if bins were made correctly
        '''bin_size = self.min_bin_size

        for bin_list_index in range(self.bin_list_numb):
            bin_list = self.bins_list[bin_list_index]
            print(bin_list_index)
            for bin_index in range(len(bin_list)):
                start = bin_index * bin_size
                end = start + bin_size
                test_avg = get_bin_value(self.value_map, self.index_list,
                                          self.intervals[0], self.intervals[1],
                                          start, end)
                if (bin_list[bin_index] == -1 and test_avg != -1) or abs(bin_list[bin_index] - test_avg) > 0.0000001:
                    print(bin_index, bin_list[bin_index], test_avg)
                    exit(-1)
            bin_size *= 2'''

    # get the specific stat to search for
    def get_method(self, stat):
        if stat == "mean":
            return self.get_exact_mean
        elif stat == "approx_mean":
            if self.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return self.get_approx_mean
        elif stat == "median":
            print("Median has not been implemented yet")
            return None
            # return self.get_median
        elif stat == "max":
            return self.get_max
        elif stat == "min":
            return self.get_min
        elif stat == "coverage":
            return self.get_coverage
        elif stat == "std":
            return self.get_std
        else:
            print(f"{stat} is not a valid statistic to search for")
            return None

    def get_approx_mean(self, start_list, end_list):
        return get_approx_means(self.bins_list[0], self.bins_list_coverages[0],
                                self.min_bin_size, start_list, end_list)

    def get_exact_mean(self, start_list, end_list):
        return get_exact_means(self.value_map, self.index_list, self.intervals[0],
                               self.intervals[1], start_list, end_list)

    # slower for now because of usage of numpy instead of implementing O(n)
    # algorithm in Cython
    def get_median(self, start_list, end_list):
        assert len(start_list) == len(end_list)
        num_tests = len(start_list)
        results = np.full(num_tests, -1, dtype=np.float64)
        for i in range(num_tests):
            start = start_list[i]
            end = end_list[i]
            wanted_range = self.value_list[start:end]
            cleaned_range = wanted_range[wanted_range > -1]
            results[i] = np.median(cleaned_range)
        return results

    def get_coverage(self, start_list, end_list):
        return get_coverages(self.index_list, self.intervals[0],
                             self.intervals[1], start_list, end_list)

    def get_max(self, start_list, end_list):
        return get_maximums(self.value_map, self.index_list, self.intervals[0],
                            self.intervals[1], start_list, end_list)

    def get_min(self, start_list, end_list):
        return get_minimums(self.value_map, self.index_list, self.intervals[0],
                            self.intervals[1], start_list, end_list)

    def get_std(self, start_list, end_list):
        return get_stds(self.value_map, self.index_list, self.intervals[0],
                        self.intervals[1], start_list, end_list)
