import numpy as np
import time
import math

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

MAX_NUMB_BIN_LIST = 6
MIN_BIN_SIZE = 3


class Chrom_Data:

    def __init__(self, name, size, max_bin_size=None):
        self.name = name
        self.size = size

        # don't use until user loads this chromosome for searching
        self.value_list = None
        self.test_value_list = None

        # length is # of intervals in bedgraph file for the chromosome
        self.value_map = np.zeros(size, dtype=np.float64)
        self.intervals = [
            np.zeros(size, dtype=np.uint32),
            np.zeros(size, dtype=np.uint32)
        ]
        self.current_index = 0
        self.max_index = 0
        self.loaded_value_list = False
        self.loaded_bins = False

        self.max_bin_size = max_bin_size
        self.bins_list = []
        self.bins_list_coverages = []
        self.bin_list_numb = 0

        print(f"Initialized {name}.")

    def add_data(self, data):
        start = int(data[START_INDEX])
        end = int(data[END_INDEX])
        value = float(data[VALUE_INDEX].strip())

        self.intervals[0][self.current_index] = start
        self.intervals[1][self.current_index] = end

        self.value_map[self.current_index] = value
        self.current_index += 1

        self.max_index = end

    def trim_extra_space(self):
        self.value_map = self.value_map[:self.current_index]
        self.intervals[0] = self.intervals[0][:self.current_index]
        self.intervals[1] = self.intervals[1][:self.current_index]

    def load_value_array(self):

        print(f"Loading {self.name}...")

        if self.loaded_value_list:
            print(f"{self.name} is already loaded")
            return

        # create and populate the value_list
        self.value_list = np.full(self.size, -1, dtype=np.float64)

        for i in range(self.intervals[0].size):
            start = self.intervals[0][i]
            end = self.intervals[1][i]
            self.value_list[start:end] = self.value_map[i]

        if not self.loaded_bins and self.max_bin_size is not None:
            self.split_bins()

        self.loaded_value_list = True

        print(f"Done with loading {self.name}")

    # approximate
    def free_value_array(self):
        self.value_list = None
        self.loaded_value_list = False
        print(f"Freed memory for {self.name}'s value_list")

    def split_bins(self, max_bin_size=None):

        if max_bin_size is not None:
            if max_bin_size == self.max_bin_size:
                print(f"Already has split bins for: {max_bin_size}")
                return

            self.max_bin_size = max_bin_size

        print(f"Splitting bins for {self.name} into size: {self.max_bin_size}...")

        self.bins_list.clear()
        self.bins_list_coverages.clear()

        bin_size_divisor = 2 ** (MAX_NUMB_BIN_LIST - 1)

        self.bin_list_numb = 0
        for i in range(MAX_NUMB_BIN_LIST):

            # should always be an integer, if not, make sure max_bin_size can be divided
            bin_size = self.max_bin_size / bin_size_divisor
            bin_size_divisor /= 2

            if bin_size < MIN_BIN_SIZE or bin_size != int(bin_size):
                continue

            bin_size = int(bin_size)

            numb_bins = math.ceil(self.size / bin_size)

            # normal python array is twice as fast
            # bins = np.full(numb_bins, -1, dtype=np.float)
            # bins_coverage = np.full(numb_bins, -1, dtype=np.float)
            bins = [-1.0] * numb_bins
            bins_coverage = [-1.0] * numb_bins

            print(f"Bin size: {bin_size}")
            print(f"# of bins: {numb_bins}")

            for current_bin_index in range(numb_bins):
                mean = None
                coverage = None

                # first (smallest) bin list
                if len(self.bins_list) == 0:
                    start = current_bin_index * bin_size
                    end = start + bin_size
                    mean = self.get_exact_mean(start, end)
                    coverage = self.get_coverage(start, end)
                else:
                    smaller_bin_index = current_bin_index * 2

                    bin1_value = self.bins_list[self.bin_list_numb - 1][smaller_bin_index]
                    bin1_coverage = self.bins_list_coverages[self.bin_list_numb - 1][smaller_bin_index]

                    # might not have a bin2 for the last bigger bin to make
                    bin2_value = -1
                    bin2_coverage = 0
                    if smaller_bin_index + 1 < len(self.bins_list[self.bin_list_numb - 1]):
                        bin2_value = self.bins_list[self.bin_list_numb - 1][smaller_bin_index + 1]
                        bin2_coverage = self.bins_list_coverages[self.bin_list_numb - 1][smaller_bin_index + 1]

                    if bin1_value != -1 and bin2_value != -1:
                        mean = (bin2_value * bin2_coverage + bin1_value * bin1_coverage) /\
                               (bin2_coverage + bin1_coverage)
                    elif bin2_value == -1:
                        mean = bin1_value
                    elif bin1_value == -1:
                        mean = bin2_value

                    coverage = (bin1_coverage + bin2_coverage) / 2

                if mean is None:
                    mean = -1

                bins[current_bin_index] = float(mean)
                bins_coverage[current_bin_index] = coverage

            self.bins_list.append(bins)
            self.bins_list_coverages.append(bins_coverage)
            self.bin_list_numb += 1

        self.loaded_bins = True
        print("Done with splitting bins")

    def get_exact_mean(self, start, end):
        wanted_range = self.value_list[start:end]
        cleaned_range = wanted_range[wanted_range > -1]
        if cleaned_range.size == 0:
            return None

        return np.mean(cleaned_range)

    def get_approx_mean(self, start, end):
        bin_start = int(start / self.max_bin_size)
        bin_end = int(end / self.max_bin_size)

        max_size_bin = self.bins_list[self.bin_list_numb - 1]

        bin_index = bin_start

        # special case where interval is within a single bin
        if bin_start == bin_end:
            if max_size_bin[bin_index] == -1:
                return None
            return max_size_bin[bin_index]

        mean_value = 0
        numb_value = 0

        # first bin
        if max_size_bin[bin_index] != -1:
            weight = (bin_index + 1) * self.max_bin_size - start
            mean_value += weight * max_size_bin[bin_index]
            numb_value += weight

        # middle bins
        weight = self.max_bin_size
        bin_index += 1
        while bin_index < bin_end:
            if max_size_bin[bin_index] != -1:
                mean_value += weight * max_size_bin[bin_index]
                numb_value += weight

            bin_index += 1

        # last bin
        if max_size_bin[bin_end] != -1:
            weight = end - bin_end * self.max_bin_size
            mean_value += weight * max_size_bin[bin_end]
            numb_value += weight

        if mean_value == 0:
            return None

        mean_value /= numb_value
        return mean_value

    def get_mod_approx_mean(self, start, end):
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

    def get_mod2_approx_mean(self, start, end):
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
                #bin_weight = weight * self.bins_list_coverages[self.bin_list_numb - 1][max_bin_index]
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

    def get_median(self, start, end):
        my_range = self.value_list[start:end][self.value_list[start:end] > -1]
        if my_range.size == 0:
            return None
        return np.median(my_range)

    def get_coverage(self, start, end):
        my_range = self.value_list[start:end][self.value_list[start:end] > -1]
        orig_range = end - start
        return 1.0 - (orig_range - my_range.size) / orig_range

    def get_max(self, start, end):
        my_range = self.value_list[start:end][self.value_list[start:end] > -1]
        if my_range.size == 0:
            return None
        return np.amax(my_range)

    def get_min(self, start, end):
        my_range = self.value_list[start:end][self.value_list[start:end] > -1]
        if my_range.size == 0:
            return None

        return np.amin(my_range)

    def get_std(self, start, end):
        my_range = self.value_list[start:end][self.value_list[start:end] > -1]
        if my_range.size == 0:
            return None
        return np.std(my_range)

    '''
    def get_all_mean(self, test_cases):
        length_arr = []
        start_list = [i[0] for i in test_cases]
        end_list = [i[1] for i in test_cases]
        for i in range(len(test_cases) + 1):
            length_arr.append([])

        start_time = time.time()
        for i in range(len(test_cases)):
            start = start_list[i]
            end = end_list[i]
            wanted_arr = self.value_list[start:end][self.value_list[start:end] > -1]
            length_arr[wanted_arr.size].append(wanted_arr)

        print(time.time() - start_time)
        count = 0
        for i in range(len(length_arr)):
            if len(length_arr[i]) == 0:
                continue
            if i == 0:
                continue
            count += 1
            np.std(length_arr[i], axis=1)

        print(count)'''
