from ChromosomeNumpy import ChromosomeNumpy
import numpy as np
import pyBigWig

CHROM_NAME_INDEX = 0


class Data:

    def __init__(self, data_file_name, chromosome_name=None):

        current = None

        # could make chromosomes a map if there are many
        self.chromosomes = {}
        self.bins = None

        print(f"Reading in {data_file_name}...")

        with open(data_file_name, 'r') as data_file:

            for line in data_file:
                data = line.split("\t")

                if chromosome_name is not None and chromosome_name != data[CHROM_NAME_INDEX]:

                    # have not created chromosome and added data to it yet
                    if current is None:
                        continue

                    # done with adding data, stop now
                    else:
                        break

                # update current to be the chromosome to add data to
                if current is None or data[CHROM_NAME_INDEX] != current.name:
                    current = ChromosomeNumpy(data[CHROM_NAME_INDEX])
                    self.chromosomes[current.name] = current

                current.add_data(data)

        print(f"Finished")

        """
        # Tests if chromosome is inputted correctly for bin size 1
        bw = pyBigWig.open('data/P2MC7N8HCE3K.bw')
        ranges = bw.intervals("chr2L", 0, 1000001)
        range_length = len(ranges)
        for i in range(range_length):
            for j in range(ranges[i][0], ranges[i][1]):
                if current.bins[j] != ranges[i][2]:
                    print(ranges[i], current.bins[j])
        exit(-1)
        """

    def get_coverage(self, chromosome_name, start, end):
        chrom = self.chromosomes[chromosome_name]

        my_range = chrom.values[start:end][chrom.values[start:end] > -1]
        orig_range = end - start
        return 1.0 - (orig_range - my_range.size) / orig_range

    def get_exact_average(self, chromosome_name, start, end):
        chrom = self.chromosomes[chromosome_name]

        wanted_range = chrom.values[start:end]
        cleaned_range = wanted_range[wanted_range > -1]
        if cleaned_range.size == 0:
            return None

        return np.mean(cleaned_range)

    def get_approx_average(self, chromosome_name, start, end):
        chrom = self.chromosomes[chromosome_name]

        self.bins = chrom.bins
        bins = self.bins
        bin_size = chrom.bin_size

        bin_start = int(start / bin_size)
        bin_end = int(end / bin_size)

        """
        # only works for numpy arrays
        my_range = (bins[bin_start:bin_end])[bins[bin_start:bin_end] > -1]
        if my_range.size == 0:
            return None

        return np.mean(my_range)
        """

        bin_index = bin_start

        # special case where interval is within a bin
        if bin_start == bin_end:
            return bins[bin_index]

        average_value = 0
        numb_value = 0

        # first bin
        if bins[bin_index] != -1:
            weight = (bin_index + 1) * bin_size - start
            average_value += weight * bins[bin_index]
            numb_value += weight

        # middle bins
        weight = bin_size
        bin_index += 1
        tries = 0
        while bin_index < bin_end:
            if bins[bin_index] != -1:
                average_value += weight * bins[bin_index]
                numb_value += weight

            bin_index += 1
            tries += 1

        # last bin
        if bins[bin_end] != -1:
            weight = end - bin_end * bin_size
            average_value += weight * bins[bin_end]
            numb_value += weight

        if average_value == 0:
            return None

        average_value /= numb_value
        return average_value
