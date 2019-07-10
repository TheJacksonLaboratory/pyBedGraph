import numpy as np
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
        self.bins_list_coverages = None

    def initialize_value_array(self):
        self.value_list = np.zeros(self.size, dtype=np.float64)

    def load_bins(self, max_bin_size):

        if max_bin_size is None:
            print("Did not specify max_bin_size")
            return

        if max_bin_size == self.max_bin_size and self.loaded_bins is True:
            print(f"Already loaded bins for: {max_bin_size}")
            return

        self.max_bin_size = max_bin_size

        self.bins_list.clear()

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

            bins_list = load_bins(prev_bins_list)
            prev_bins_list = bins_list

            self.bins_list.append(bins_list)

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

    def get_exact_mean(self, start_list, end_list):
        return get_exact_means(self.value_list,
                               self.bins_list[self.bin_list_numb - 1],
                               self.max_bin_size, start_list, end_list)

    def get_median(self, start_list, end_list):
        '''results = np.zeros(len(start_list), dtype=np.float64)
        for i in range(len(start_list)):
            results[i] = np.median(self.value_list[start_list[i]:end_list[i]])

        return results'''

        return get_medians(self.value_list, start_list, end_list)

    def get_coverage(self, start_list, end_list):
        return get_coverages(self.value_list, start_list, end_list)

    def get_max(self, start_list, end_list):
        return get_maximums(self.value_list, start_list, end_list)

    def get_min(self, start_list, end_list):
        return get_minimums(self.value_list, start_list, end_list)

    # TODO
    def get_std(self, start_list, end_list):
        return None

    # maybe 3x faster to vectorize not including time to create matrix
    '''def get_all_mean(self, test_cases):
        length_arr = []

        for i in range(len(test_cases)):
            start = start_list[i]
            end = end_list[i]
            wanted_arr = self.value_list[start:end][self.value_list[start:end] > -1]
            length_arr[wanted_arr.size].append(wanted_arr)

        count = 0
        for i in range(len(length_arr)):
            if len(length_arr[i]) == 0:
                continue
            if i == 0:
                continue
            count += 1
            np.std(length_arr[i], axis=1)

        print(count)'''
