import numpy as np
import math

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

CHROM_MAX_SIZE = 250000000
MAX_NUMB_BIN_LIST = 6
MIN_BIN_SIZE = 64


class Chrom_Data:

    def __init__(self, name, size):
        self.name = name
        self.size = size

        self.values = np.full(size, -1, dtype=np.float64)
        self.value_indexes = np.full(size, -1, dtype=np.uint32)
        self.size = 0

        self.max_bin_size = None
        self.bins = []
        self.bin_list_numb = MAX_NUMB_BIN_LIST

        print(f"Created {name}")

    def add_data(self, data):
        start = int(data[START_INDEX])
        end = int(data[END_INDEX])
        value = float(data[VALUE_INDEX].strip())

        self.values[start:end] = value

        self.size = end

    def split_bins(self, max_bin_size):

        print(f"Splitting bins for {self.name}...")

        self.max_bin_size = max_bin_size
        bin_size_divisor = 2 ** (MAX_NUMB_BIN_LIST - 1)

        prev_bin_list = None
        for i in range(MAX_NUMB_BIN_LIST):

            # should always be an integer, if not, make sure max_bin_size can be divided
            bin_size = max_bin_size / bin_size_divisor
            bin_size_divisor /= 2

            if bin_size != int(bin_size):
                print(f"Bin size: {bin_size} must be an integer")
                exit(-1)
            bin_size = int(bin_size)

            if bin_size < MIN_BIN_SIZE:
                self.bin_list_numb -= 1
                continue

            numb_bins = math.ceil(self.size / bin_size)

            # normal python array seems faster
            bins = np.full(numb_bins, -1, dtype=np.float64)

            print(f"Bin size: {bin_size}")
            print(f"# of bins: {numb_bins}")

            for bin_index in range(numb_bins):
                mean = None

                # first bin
                if len(self.bins) == 0:
                    mean = self.get_exact_mean(bin_index * bin_size,
                                               (bin_index + 1) * bin_size)
                else:
                    new_index = bin_index * 2

                    bin1_value = prev_bin_list[new_index]
                    bin1_start = bin_index * bin_size
                    bin1_end = bin1_start + bin_size // 2
                    bin1_coverage = self.get_coverage(bin1_start, bin1_end)

                    bin2_value = -1
                    bin2_coverage = -1
                    if new_index + 1 < prev_bin_list.size:
                        bin2_value = prev_bin_list[new_index + 1]
                        bin2_start = bin1_end
                        bin2_end = bin2_start + bin_size // 2
                        bin2_coverage = self.get_coverage(bin2_start, bin2_end)

                    if bin1_value != -1 and bin2_value != -1:
                        mean = (bin2_value * bin2_coverage + bin1_value * bin1_coverage) /\
                               (bin2_coverage + bin1_coverage)
                    elif bin2_value == -1:
                        mean = bin1_value
                    elif bin1_value == -1:
                        mean = bin2_value

                if mean is None:
                    mean = -1
                bins[bin_index] = mean

            self.bins.append(bins)
            prev_bin_list = bins

        print("Done.\n")

    def get_exact_mean(self, start, end):
        wanted_range = self.values[start:end]
        cleaned_range = wanted_range[wanted_range > -1]
        if cleaned_range.size == 0:
            return None

        return np.mean(cleaned_range)

    def get_approx_mean(self, start, end):
        bin_start = int(start / self.max_bin_size)
        bin_end = int(end / self.max_bin_size)

        max_size_bin = self.bins[self.bin_list_numb - 1]

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
        tries = 0
        while bin_index < bin_end:
            if max_size_bin[bin_index] != -1:
                mean_value += weight * max_size_bin[bin_index]
                numb_value += weight

            bin_index += 1
            tries += 1

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

        max_size_bin = self.bins[self.bin_list_numb - 1]

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
            while weight < current_bin_size / 2:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2 + 1

                if bin_size_index == 0:
                    break

            if self.bins[bin_size_index][bin_index] != -1:
                mean_value += weight * self.bins[bin_size_index][bin_index]
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
            while weight < current_bin_size / 2:
                current_bin_size /= 2
                bin_size_index -= 1
                bin_index = bin_index * 2

                if bin_size_index == 0:
                    break

            if self.bins[bin_size_index][bin_index] != -1:
                mean_value += weight * self.bins[bin_size_index][bin_index]
                numb_value += weight

        if mean_value == 0:
            return None

        mean_value /= numb_value
        return mean_value

    def get_coverage(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        orig_range = end - start
        return 1.0 - (orig_range - my_range.size) / orig_range

    def get_max(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        if my_range.size == 0:
            return None
        return np.amax(my_range)

    def get_min(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        if my_range.size == 0:
            return None

        return np.amin(my_range)

    def get_std(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        if my_range.size == 0:
            return None
        return np.std(my_range)
