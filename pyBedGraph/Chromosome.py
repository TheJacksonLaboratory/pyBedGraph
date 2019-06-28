import numpy as np
import math

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

CHROM_MAX_SIZE = 250000000
NUMB_BIN_LIST = 2


class Chromosome:

    def __init__(self, name):
        self.name = name

        self.values = np.full(CHROM_MAX_SIZE, -1, dtype=np.float64)
        self.size = 0

        self.max_bin_size = None
        self.bins = []
        for i in range(NUMB_BIN_LIST):
            self.bins.append([])

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

        for i in range(NUMB_BIN_LIST):

            # should always be an integer, if not, make sure max_bin_size can be divided
            bin_size = max_bin_size / (NUMB_BIN_LIST - i)

            if bin_size != int(bin_size):
                print(f"Bin size: {bin_size} must be an integer")
                exit(-1)
            bin_size = int(bin_size)

            numb_bins = math.ceil(self.size / bin_size)

            # normal python array seems faster
            self.bins[i] = [-1.0] * numb_bins
            # self.bins = np.full(int(CHROM_SIZE / bin_size), -1, dtype=np.float64)

            print(f"Bin size: {bin_size}")
            print(f"# of bins: {numb_bins}")

            for bin_index in range(numb_bins):
                mean = self.get_exact_mean(bin_index * bin_size,
                                           (bin_index + 1) * bin_size)
                if mean is None:
                    mean = -1
                self.bins[i][bin_index] = mean

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

        max_size_bin = self.bins[NUMB_BIN_LIST - 1]

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

        # middle self.bins
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

        max_size_bin = self.bins[NUMB_BIN_LIST - 1]

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
            if weight < self.max_bin_size / 2:
                small_bin_index = int(start / (self.max_bin_size / 2))
                if (small_bin_index + 1) * self.max_bin_size / 2 - start != weight:
                    print(weight, (small_bin_index + 1) * self.max_bin_size / 2 - start)
                    exit(-1)
                if self.bins[0][small_bin_index] != -1:
                    mean_value += weight * self.bins[0][small_bin_index]
                    numb_value += weight
            else:
                mean_value += weight * max_size_bin[bin_index]
                numb_value += weight

        # middle self.bins
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
            if weight < self.max_bin_size / 2:
                small_bin_index = bin_end * 2
                if self.bins[0][small_bin_index] != -1:
                    mean_value += weight * self.bins[0][small_bin_index]
                    numb_value += weight
            else:
                mean_value += weight * max_size_bin[bin_end]
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
