import numpy as np
import math

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3

CHROM_SIZE = 250000000


class ChromosomeNumpy:

    def __init__(self, name):
        self.name = name
        self.data_index = 0

        self.values = np.full(CHROM_SIZE, -1, dtype=np.float64)
        self.current_index = 0

        self.bin_size = None
        self.bins = None

        print(f"Created {name}")

    def add_data(self, data):
        start = int(data[START_INDEX])
        end = int(data[END_INDEX])
        value = float(data[VALUE_INDEX].strip())

        self.values[start:end] = value

        self.current_index = end

    def get_average(self, start, end):
        my_range = self.values[start:end][self.values[start:end] > -1]
        if my_range.size == 0:
            return -1

        return np.mean(my_range)

    def split_bins(self, bin_size):

        print(f"Splitting bins for {self.name}")

        total_size = self.current_index
        numb_bins = math.ceil(total_size / bin_size)

        self.bin_size = bin_size
        #self.bins = [-1.0] * numb_bins
        self.bins = np.full(int(CHROM_SIZE / bin_size), -1, dtype=np.float64)

        print(f"Bin size: {self.bin_size}")
        print(f"# of bins: {numb_bins}")

        for bin_index in range(numb_bins):
            self.bins[bin_index] = self.get_average(bin_index * bin_size,
                                                    (bin_index + 1) * bin_size)

        with open('test.txt', 'w') as test:
            output = ""
            for i in range(numb_bins):
                output += f'{self.bins[i]}\n'

            test.write(output)
        """
            for value_index in range(bin_size):
                allele_index = value_index + bin_index * self.bin_size

                # nothing to do before, so go to the start of the window
                if allele_index < current_window.start:
                    value_index = current_window.start - bin_index * self.bin_size
                    continue

                # go to the next window
                if allele_index >= current_window.end:
                    data_index += 1

                    # no more windows
                    if len(self.windows) == data_index:
                        done = True
                        break

                    current_window = self.windows[data_index]

                    if allele_index > current_window.end:
                        print("ERROR in finding current_window")

                    continue

                # if bin_index < 52:
                # print(f"{i} | {allele_index}: {current_window.value}")

                average_value += current_window.value
                numb_values += 1

                value_index += 1

            if numb_values > 0:
                average_value /= numb_values
                self.values.append(average_value)

            else:
                self.values.append(-1)

            #if bin_index > 5120:
            #    print(f"From: average{bin_index * bin_size} -> {bin_index * bin_size + bin_size} | Bin value: {self.bins[-1]}")

            if done:

                # pad the rest of the list with empty bins
                while len(self.values) < numb_bins:
                    self.values.append(-1)

                break

        print()
        """
