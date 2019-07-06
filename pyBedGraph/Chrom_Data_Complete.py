import numpy as np
from .Chrom_Data import Chrom_Data
from .complete_stats import *

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

MAX_NUMB_BIN_LIST = 1
MIN_BIN_SIZE = 3


class Chrom_Data_Complete(Chrom_Data):

    def __init__(self, name, size):
        super().__init__(name, size)
        self.bins_list_coverages = None
        self.bins_list = None
        self.loaded_bins = True

    def load_value_array(self):

        print(f"Loading {self.name}...")

        if self.loaded_value_list:
            print(f"{self.name} is already loaded")
            return

        # create and populate the value_list
        self.value_list = np.zeros(self.size, dtype=np.float64)

        for i in range(self.intervals[0].size):
            start = self.intervals[0][i]
            end = self.intervals[1][i]
            self.value_list[start:end] = self.value_map[i]

        self.loaded_value_list = True

        print(f"Done with loading {self.name}")

    def get_exact_mean(self, start_list, end_list):
        return get_exact_mean(self.value_list, start_list, end_list)

    # TODO
    def get_median(self, start_list, end_list):
        return None

    def get_coverage(self, start_list, end_list):
        return get_coverage(self.value_list, start_list, end_list)

    def get_max(self, start_list, end_list):
        return get_max(self.value_list, start_list, end_list)

    def get_min(self, start_list, end_list):
        return get_min(self.value_list, start_list, end_list)

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
