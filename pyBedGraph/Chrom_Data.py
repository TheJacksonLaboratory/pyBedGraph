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

        # don't use until user loads this chromosome for searching
        self.value_list = None
        self.test_value_list = None

        # length is # of intervals in bedGraph file for the chromosome
        self.value_map = np.zeros(size, dtype=np.float64)
        self.intervals = [
            np.zeros(size, dtype=np.uint32),
            np.zeros(size, dtype=np.uint32)
        ]
        self.current_index = 0
        self.max_index = 0
        self.loaded_value_list = False
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

        self.value_map[self.current_index] = value
        self.current_index += 1

        self.max_index = end

    def trim_extra_space(self):
        self.value_map = self.value_map[:self.current_index]
        self.intervals[0] = self.intervals[0][:self.current_index]
        self.intervals[1] = self.intervals[1][:self.current_index]

    def initialize_value_array(self):
        self.value_list = np.full(self.size, -1, dtype=np.float64)

    def load_value_array(self):

        print(f"Loading {self.name} ...")

        if self.loaded_value_list:
            print(f"{self.name} is already loaded")
            return

        self.initialize_value_array()
        fill_value_array(self.intervals[0], self.intervals[1], self.value_map,
                         self.value_list)

        self.loaded_value_list = True

    def free_value_array(self):
        self.value_list = None
        self.loaded_value_list = False
        print(f"Freed memory for {self.name}'s value_list")

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
        prev_bins_list, prev_bins_coverage_list = load_smallest_bins(self.value_list, bin_size)
        self.bins_list.append(prev_bins_list)
        self.bins_list_coverages.append(prev_bins_coverage_list)
        bin_size *= 2

        # Load larger bins
        while bin_size <= max_bin_size:

            print(f"Loading bins of size: {bin_size} for {self.name} ...")
            print(f"Number of bins: {math.ceil(self.value_list.size / bin_size)}")

            bins_list, bins_coverage_list = load_bins(prev_bins_list,
                                                      prev_bins_coverage_list)
            prev_bins_list = bins_list
            prev_bins_coverage_list = bins_coverage_list

            self.bins_list.append(bins_list)
            self.bins_list_coverages.append(bins_coverage_list)

            bin_size *= 2

        self.bin_list_numb = len(self.bins_list)
        self.loaded_bins = True

        '''for i in range(self.bin_list_numb):
            bin_size /= 2

        # test if bins were made correctly
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

    '''def get_mod_approx_mean(self, start_list, end_list):
        bin_start = int(start / self.max_bin_size)
        bin_end = int(end / self.max_bin_size)

        max_size_bin = self.bins_list[self.bin_list_numb - 1]

        max_bin_index = bin_start

        # special case where interval is within a single bin
        if bin_start == bin_end:
            if max_size_bin[max_bin_index] == -1:
                return None
            return max_size_bin[max_bin_index]

        mean_value = 0
        numb_value = 0

        # first bin
        if max_size_bin[max_bin_index] != -1:
            weight = (max_bin_index + 1) * self.max_bin_size - start

            current_bin_size = self.max_bin_size
            bin_size_index = self.bin_list_numb - 1
            bin_index = max_bin_index
            while weight < current_bin_size / 2 and bin_size_index > 0:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2 + 1

            if self.bins_list[bin_size_index][bin_index] != -1:
                mean_value += weight * self.bins_list[bin_size_index][bin_index]
                numb_value += weight

        # middle bins
        weight = self.max_bin_size
        max_bin_index += 1
        while max_bin_index < bin_end:
            if max_size_bin[max_bin_index] != -1:
                mean_value += weight * max_size_bin[max_bin_index]
                numb_value += weight

            max_bin_index += 1

        # last bin
        if max_size_bin[max_bin_index] != -1:
            weight = end - max_bin_index * self.max_bin_size

            current_bin_size = self.max_bin_size
            bin_size_index = self.bin_list_numb - 1
            bin_index = max_bin_index
            while weight < current_bin_size / 2 and bin_size_index > 0:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2

            if self.bins_list[bin_size_index][bin_index] != -1:
                mean_value += weight * self.bins_list[bin_size_index][bin_index]
                numb_value += weight

        if mean_value == 0:
            return None

        mean_value /= numb_value
        return mean_value

    def get_mod2_approx_mean(self, start_list, end_list):
        bin_start = int(start / self.max_bin_size)
        bin_end = int(end / self.max_bin_size)

        max_size_bin = self.bins_list[self.bin_list_numb - 1]

        max_bin_index = bin_start

        # special case where interval is within a single bin
        if bin_start == bin_end:
            if max_size_bin[max_bin_index] == -1:
                return None
            return max_size_bin[max_bin_index]

        mean_value = 0
        numb_value = 0

        # first bin
        if max_size_bin[max_bin_index] != -1:
            weight = (max_bin_index + 1) * self.max_bin_size - start

            current_bin_size = self.max_bin_size
            bin_size_index = self.bin_list_numb - 1
            bin_index = max_bin_index
            while weight < current_bin_size / 2 and bin_size_index > 0:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2 + 1

            if self.bins_list[bin_size_index][bin_index] != -1:
                weight *= self.bins_list_coverages[bin_size_index][bin_index]
                mean_value += weight * self.bins_list[bin_size_index][bin_index]
                numb_value += weight

        # middle bins
        weight = self.max_bin_size
        max_bin_index += 1
        while max_bin_index < bin_end:
            if max_size_bin[max_bin_index] != -1:
                bin_weight = weight * self.bins_list_coverages[self.bin_list_numb - 1][max_bin_index]
                #print(self.bins_list_coverages[self.bin_list_numb - 1][max_bin_index])
                bin_weight = weight
                mean_value += bin_weight * max_size_bin[max_bin_index]
                numb_value += bin_weight

            max_bin_index += 1

        # last bin
        if max_size_bin[max_bin_index] != -1:
            weight = end - max_bin_index * self.max_bin_size

            current_bin_size = self.max_bin_size
            bin_size_index = self.bin_list_numb - 1
            bin_index = max_bin_index
            while weight < current_bin_size / 2 and bin_size_index > 0:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2

            if self.bins_list[bin_size_index][bin_index] != -1:
                weight *= self.bins_list_coverages[bin_size_index][bin_index]
                mean_value += weight * self.bins_list[bin_size_index][bin_index]
                numb_value += weight

        if mean_value == 0:
            return None

        mean_value /= numb_value
        return mean_value

    '''

    def get_method(self, stat):
        if stat == "mean":
            if self.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return self.get_exact_mean
        elif stat == "approx_mean":
            if self.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return self.get_approx_mean
        elif stat == "mod_approx_mean":
            if self.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return self.get_mod_approx_mean
        elif stat == "median":
            return self.get_median
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

    def get_mod_approx_mean(self, start_list, end_list):
        return get_mod_approx_means(self.bins_list, self.max_bin_size,
                                    self.bin_list_numb, start_list, end_list)

    def get_exact_mean(self, start_list, end_list):
        return get_exact_means(self.value_list,
                               self.bins_list[self.bin_list_numb - 1],
                               self.max_bin_size,
                               self.bins_list_coverages[self.bin_list_numb - 1],
                               start_list, end_list)

    # TODO
    def get_median(self, start_list, end_list):
        return None

    def get_coverage(self, start_list, end_list):
        return get_coverages(self.value_list, start_list, end_list)

    def get_max(self, start_list, end_list):
        return get_maximums(self.value_list, start_list, end_list)

    def get_min(self, start_list, end_list):
        return get_minimums(self.value_list, start_list, end_list)

    # TODO
    def get_std(self, start_list, end_list):
        # my_range = self.value_list[start:end][self.value_list[start:end] > -1]
        # if my_range.size == 0:
        #    return None
        # return np.std(my_range)
        return None
