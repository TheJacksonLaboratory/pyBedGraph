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

    def get_average(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        if my_range.size == 0:
            return -1

        return np.mean(my_range)

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
                self.bins[i][bin_index] = self.get_average(bin_index * bin_size,
                                                           (bin_index + 1) * bin_size)

        print("Finished\n")

    def get_exact_average(self, start, end):
        wanted_range = self.values[start:end]
        cleaned_range = wanted_range[wanted_range > -1]
        if cleaned_range.size == 0:
            return None

        return np.mean(cleaned_range)

    def get_approx_average(self, start, end):
        bin_start = int(start / self.max_bin_size)
        bin_end = int(end / self.max_bin_size)

        max_size_bin = self.bins[NUMB_BIN_LIST - 1]

        bin_index = bin_start

        # special case where interval is within a single bin
        if bin_start == bin_end:
            return max_size_bin[bin_index]

        average_value = 0
        numb_value = 0

        # first bin
        if max_size_bin[bin_index] != -1:
            weight = (bin_index + 1) * self.max_bin_size - start
            average_value += weight * max_size_bin[bin_index]
            numb_value += weight

        # middle self.bins
        weight = self.max_bin_size
        bin_index += 1
        tries = 0
        while bin_index < bin_end:
            if max_size_bin[bin_index] != -1:
                average_value += weight * max_size_bin[bin_index]
                numb_value += weight

            bin_index += 1
            tries += 1

        # last bin
        if max_size_bin[bin_end] != -1:
            weight = end - bin_end * self.max_bin_size
            average_value += weight * max_size_bin[bin_end]
            numb_value += weight

        if average_value == 0:
            return None

        average_value /= numb_value
        return average_value

    def get_mod_approx_average(self, start, end):
        bin_start = int(start / self.max_bin_size)
        bin_end = int(end / self.max_bin_size)

        max_size_bin = self.bins[NUMB_BIN_LIST - 1]

        bin_index = bin_start

        # special case where interval is within a single bin
        if bin_start == bin_end:
            return max_size_bin[bin_index]

        average_value = 0
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
                    average_value += weight * self.bins[0][small_bin_index]
                    numb_value += weight
            else:
                average_value += weight * max_size_bin[bin_index]
                numb_value += weight

        # middle self.bins
        weight = self.max_bin_size
        bin_index += 1
        tries = 0
        while bin_index < bin_end:
            if max_size_bin[bin_index] != -1:
                average_value += weight * max_size_bin[bin_index]
                numb_value += weight

            bin_index += 1
            tries += 1

        # last bin
        if max_size_bin[bin_end] != -1:
            weight = end - bin_end * self.max_bin_size
            average_value += weight * max_size_bin[bin_end]
            numb_value += weight

        if average_value == 0:
            return None

        average_value /= numb_value
        return average_value

    def get_coverage(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        orig_range = end - start
        return 1.0 - (orig_range - my_range.size) / orig_range
