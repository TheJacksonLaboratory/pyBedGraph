# .Dependency for uploading to pypi
# from .Chromosome import Chromosome

# normal import for local use
from .Chrom_Data import Chrom_Data
from .Chrom_Data_Complete import Chrom_Data_Complete
import numpy as np
import time

CHROM_NAME_INDEX = 0
BUFFER_COUNTER = 10000  # 10^4


class BedGraph:

    def __init__(self, chrom_size_file_name, data_file_name, chrom_wanted=None,
                 like_pyBigWig=True):

        self.chromosome_map = {}
        self.chrom_sizes = {}
        self.like_pyBigWig = like_pyBigWig

        print(f"Reading in {chrom_size_file_name}...", end=" ")
        with open(chrom_size_file_name) as chrom_size_file:
            for line in chrom_size_file:
                data = line.split()
                if len(data) != 2:
                    print(f"\n{chrom_size_file} has an incorrect format."
                          f"It must be in the format:\n"
                          "chr1 100000\n"
                          "chr2 50000")
                    break

                self.chrom_sizes[data[0]] = int(data[1])
        print("Done")

        current_chrom = None
        print(f"Reading in {data_file_name}...")
        with open(data_file_name) as data_file:
            for line in data_file:
                data = line.split()
                chrom_name = data[CHROM_NAME_INDEX]

                if chrom_name not in self.chrom_sizes:
                    print(f"{chrom_name} was not included in {chrom_size_file}")
                    continue

                if chrom_wanted is not None and chrom_wanted != chrom_name:
                    # have not yet created specified chromosome and added data
                    if current_chrom is None:
                        continue

                    # done with adding data, stop now
                    else:
                        break

                if current_chrom is None or chrom_name != current_chrom.name:
                    # clean up the current chromosome before moving on
                    if current_chrom is not None:
                        current_chrom.trim_extra_space()
                        print("Done")

                    if like_pyBigWig:
                        current_chrom =\
                            Chrom_Data(chrom_name, self.chrom_sizes[chrom_name])
                    else:
                        current_chrom =\
                            Chrom_Data_Complete(chrom_name,
                                                self.chrom_sizes[chrom_name])

                    self.chromosome_map[chrom_name] = current_chrom

                current_chrom.add_data(data)

            # trim for the last chromosome found in the bedGraph file
            current_chrom.trim_extra_space()
            print("Done")

            if current_chrom is None:
                print(f"{chrom_wanted} was not found in {data_file_name}")
                exit(-1)

        print(f"Done with reading {data_file_name}")

    def load_chrom_data(self, chrom_name):
        self.chromosome_map[chrom_name].load_value_array()

    def load_chrom_bins(self, chrom_name, max_bins_size):
        self.chromosome_map[chrom_name].load_bins(max_bins_size)

    def free_chrom_data(self, chrom_name):
        self.chromosome_map[chrom_name].free_value_array()

    def get_method(self, chrom_name, stat):

        if chrom_name not in self.chromosome_map:
            print(f"{chrom_name} is not a valid chromosome")
            return None

        chromosome = self.chromosome_map[chrom_name]
        if stat == "mean":
            if chromosome.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return chromosome.get_exact_mean
        elif stat == "approx_mean":
            if chromosome.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return chromosome.get_approx_mean
        elif stat == "mod_approx_mean":
            if chromosome.loaded_bins is False:
                print(f'Bins were not loaded')
                return None
            return chromosome.get_mod_approx_mean
        elif stat == "median":
            return chromosome.get_median
        elif stat == "max":
            return chromosome.get_max
        elif stat == "min":
            return chromosome.get_min
        elif stat == "coverage":
            return chromosome.get_coverage
        elif stat == "std":
            return chromosome.get_std
        else:
            print(f"{stat} is not a valid statistic to search for")
            return None

    @staticmethod
    def change_shape(intervals):
        num_tests = len(intervals)
        start_list = np.zeros(num_tests, dtype=np.int32)
        end_list = np.zeros(num_tests, dtype=np.int32)
        for i in range(num_tests):
            if len(intervals[i]) != 3:
                print(f"List given has incorrect formatting. It must be in the format:\n"
                      "[[chr1, 1, 100], [chr1, 101, 200], ...]")
                print(intervals[i])
                continue

            start_list[i] = intervals[i][1]
            end_list[i] = intervals[i][2]

        return start_list, end_list

    # can only deal with one chromosome at a time
    def stats(self, stat="mean", intervals=None, start_list=None,
              end_list=None, chrom_name=None):

        if intervals is not None:
            method_to_call = self.get_method(intervals[0][0], stat)
            current_chrom = self.chromosome_map[intervals[0][0]]

            start_list, end_list = self.change_shape(intervals)
        elif start_list is None or end_list is None or chrom_name is None:
            print("Must either have intervals or start_list, end_list, chrom_name")
            return None
        else:
            method_to_call = self.get_method(chrom_name, stat)
            current_chrom = self.chromosome_map[chrom_name]

            assert end_list.size == start_list.size

        if not current_chrom.loaded_value_list:
            print(f"{current_chrom.name} needs to be loaded before it can be searched.")
            return

        if method_to_call is None:
            return

        start_time = time.time()
        result = method_to_call(start_list, end_list)
        print(f"Time for {stat}:", time.time() - start_time)
        return result

    def stats_from_file(self, interval_file, output_file=None, stat="mean"):
        output = ""

        out_file = None
        if output_file:
            out_file = open(output_file, 'w')

        result_counter = 0

        with open(interval_file) as in_file:

            current_chrom = None
            method_to_call = None

            for line in in_file:
                interval = line.split()

                if len(interval) != 3:
                    print(f"{interval_file} has incorrect formatting. It must be in the format:\n"
                          "chr1\t100\t401\n"
                          "chr1\t600\t1000\n"
                          "...")
                    print(interval)
                    break

                if current_chrom is None or current_chrom.name != interval[0]:
                    method_to_call = self.get_method(interval[0], stat)
                    current_chrom = self.chromosome_map[interval[0]]

                if not current_chrom.loaded_value_list:
                    print(f"{current_chrom.name} needs to be loaded before it can be searched.")
                    return

                if method_to_call is None:
                    return

                result = method_to_call(int(interval[1]), int(interval[2]))

                if out_file is None:
                    print(interval[0], interval[1], interval[2], result)
                else:
                    output += f"{result}\n"
                    result_counter += 1

                    if result_counter == BUFFER_COUNTER:
                        result_counter = 0
                        out_file.write(output)
                        output = ""

        if out_file and result_counter > 0:
            out_file.write(output)

        if out_file:
            out_file.close()
