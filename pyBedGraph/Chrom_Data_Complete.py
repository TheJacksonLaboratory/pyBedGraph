from .Chrom_Data import Chrom_Data
from .include_missing_bp import *
import math

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

MAX_NUMB_BIN_LIST = 1
MIN_BIN_SIZE = 2


class Chrom_Data_Complete(Chrom_Data):

    def __init__(self, name, size):
        super().__init__(name, size)

        # coverage is not used in this class
        self.bins_list_coverages = None

    # assume that a value of 0 in bedGraph file is equivalent to missing space
    def initialize_value_array(self):
        self.value_list = np.zeros(self.size, dtype=np.float64)

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
        print(f"Number of bins: {math.ceil(self.value_list.size / bin_size)}")
        prev_bins_list = load_smallest_bins(self.value_list, bin_size)
        self.bins_list.append(prev_bins_list)
        bin_size *= 2

        # Load larger bins
        while bin_size <= max_bin_size:

            print(f"Loading bins of size: {bin_size} for {self.name} ...")
            print(f"Number of bins: {math.ceil(self.value_list.size / bin_size)}")

            # use the previous bin list to speed up the process
            bins_list = load_bins(prev_bins_list)
            prev_bins_list = bins_list

            self.bins_list.append(bins_list)

            bin_size *= 2

        self.bin_list_numb = len(self.bins_list)
        self.loaded_bins = True

        # test if bins were made correctly
        '''for i in range(self.bin_list_numb):
            bin_size /= 2

        for bin_list_index in range(self.bin_list_numb):
            bin_list = self.bins_list[bin_list_index]
            print(bin_list_index)
            for bin_index in range(len(bin_list)):
                start = bin_index * bin_size
                end = start + bin_size
                test_avg = get_bin_mean(self.value_list, start, end)
                if abs(bin_list[bin_index] - test_avg) > 0.000001:
                    print(bin_index, bin_list[bin_index], test_avg)
                    exit(-1)
            bin_size *= 2'''

    def get_exact_mean(self, start_list, end_list):
        return get_exact_means(self.value_list,
                               self.bins_list[self.bin_list_numb - 1],
                               self.max_bin_size, start_list, end_list)

    def get_approx_mean(self, start_list, end_list):
        return get_approx_means(self.bins_list[self.bin_list_numb - 1],
                                self.max_bin_size, start_list, end_list)

    # slower for now because of usage of numpy instead of implementing O(n)
    # algorithm in Cython
    def get_median(self, start_list, end_list):
        """my_results = get_medians(self.value_list, start_list, end_list)
        return my_results"""

        results = np.zeros(len(start_list), dtype=np.float64)
        for i in range(len(start_list)):
            results[i] = np.median(self.value_list[start_list[i]:end_list[i]])
        return results

    def get_coverage(self, start_list, end_list):
        return get_coverages(self.value_list, start_list, end_list)

    def get_max(self, start_list, end_list):
        return get_maximums(self.value_list, start_list, end_list)

    def get_min(self, start_list, end_list):
        return get_minimums(self.value_list, start_list, end_list)

    def get_std(self, start_list, end_list):
        return get_stds(self.value_list, self.bins_list[self.bin_list_numb - 1],
                        self.max_bin_size, start_list, end_list)
