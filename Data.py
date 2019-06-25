from Chromosome import Chromosome
import math

CHROM_NAME_INDEX = 0


class Data:

    def __init__(self, data_file_name, bin_size):

        current = None

        # could make chromosomes a map if there are many
        self.chromosomes = {}

        with open(data_file_name, 'r') as data_file:

            for line in data_file:
                data = line.split("\t")

                if current is None or data[CHROM_NAME_INDEX] != current.name:
                    current = Chromosome(data[CHROM_NAME_INDEX])
                    self.chromosomes[current.name] = current

                current.add_data(data)

        print(f"Finished reading in {data_file_name}")

        # contains floats
        self.bins = []
        self.bin_size = 0

    def split_bins(self, bin_size):
        self.bins.clear()

        chromosome = self.chromosomes["chr2L"]

        total_size = chromosome.windows[-1].end
        numb_bins = math.ceil(total_size / bin_size)
        self.bin_size = bin_size
        print(f"Bin size: {self.bin_size}")
        print(f"# of bins: {numb_bins}")

        # keeps track of current window
        data_index = 0

        done = False

        current_window = chromosome.windows[data_index]
        for bin_index in range(numb_bins):

            average_value = 0
            numb_values = 0
            i = 0

            while i < self.bin_size:
                allele_index = i + bin_index * self.bin_size

                # nothing to do before, so go to the start of the window
                if allele_index < current_window.start:
                    i = current_window.start - bin_index * self.bin_size
                    continue

                # go to the next window
                if allele_index >= current_window.end:
                    data_index += 1

                    # no more windows
                    if len(chromosome.windows) == data_index:
                        done = True
                        break

                    current_window = chromosome.windows[data_index]

                    if allele_index > current_window.end:
                        print("ERROR in finding current_window")

                    continue

                # if bin_index < 52:
                # print(f"{i} | {allele_index}: {current_window.value}")

                average_value += current_window.value
                numb_values += 1

                i += 1

            if numb_values > 0:
                average_value /= numb_values
                self.bins.append(average_value)

            else:
                self.bins.append(-1)

            #if bin_index > 5120:
            #    print(f"From: average{bin_index * bin_size} -> {bin_index * bin_size + bin_size} | Bin value: {self.bins[-1]}")

            if done:

                # pad the rest of the list with empty bins
                while len(self.bins) < numb_bins:
                    self.bins.append(-1)

                break

    def get_average(self, chromosome_name, start, end):
        chromosome = self.chromosomes[chromosome_name]

        start_bin_index = int(start / self.bin_size)
        end_bin_index = int(end / self.bin_size)

        bin_index = start_bin_index

        # special case where interval is within a bin
        if start_bin_index == end_bin_index:
            return self.bins[bin_index]

        average_value = 0
        numb_value = 0

        # first bin
        if self.bins[bin_index] != -1:
            weight = (bin_index + 1) * self.bin_size - start
            average_value += weight * self.bins[bin_index]
            numb_value += weight

        # middle bins
        weight = self.bin_size
        bin_index += 1
        while bin_index < end_bin_index:
            if self.bins[bin_index] != -1:
                average_value += weight * self.bins[bin_index]
                numb_value += weight

            bin_index += 1

        # last bin
        if self.bins[end_bin_index] != -1:
            weight = end - end_bin_index * self.bin_size
            average_value += weight * self.bins[end_bin_index]
            numb_value += weight

        if average_value == 0:
            return None

        average_value /= numb_value
        return average_value
