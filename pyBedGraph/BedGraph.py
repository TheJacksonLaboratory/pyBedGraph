# .Dependency for uploading to pypi
#from .Chromosome import Chromosome

# normal import for local use
from Chrom_Data import Chrom_Data

CHROM_NAME_INDEX = 0


class BedGraph:

    def __init__(self, chrom_size_file_name, data_file_name, chrom_name_wanted=None, bin_size=64):

        self.chromosome = None
        self.chrom_sizes = {}

        print(f"Reading in {chrom_size_file_name}")
        with open(chrom_size_file_name) as chrom_size_file:
            for line in chrom_size_file:
                data = line.split()
                if len(data) != 2:
                    print(f"{chrom_size_file} has an incorrect format. It must be in the format:\n"
                          "chr1 100000\n"
                          "chr2 50000")
                    break

                self.chrom_sizes[data[0]] = int(data[1])

        print(f"Reading in {data_file_name}...")
        with open(data_file_name) as data_file:
            for line in data_file:
                data = line.split()
                chrom_name = data[CHROM_NAME_INDEX]

                if chrom_name not in self.chrom_sizes:
                    print(f"{chrom_name} was not included in chromosome size file")
                    break

                if chrom_name_wanted is not None and chrom_name_wanted != chrom_name:

                    # have not yet created specified chromosome and added data to it
                    if self.chromosome is None:
                        continue

                    # done with adding data, stop now
                    else:
                        break

                # update self.chromosome to be the chromosome to add data to
                if self.chromosome is None or chrom_name != self.chromosome.name:
                    self.chromosome = Chrom_Data(chrom_name, self.chrom_sizes[chrom_name])

                self.chromosome.add_data(data)

            if self.chromosome is None:
                print(f"{chrom_name_wanted} was not found in {data_file_name}")
                exit(-1)

        print("Done\n")

        self.split_bins(bin_size)

    def split_bins(self, bin_size):
        self.chromosome.split_bins(bin_size)

    def get_method(self, stat):
        if stat == "mean":
            return self.chromosome.get_exact_mean
        elif stat == "approx_mean":
            return self.chromosome.get_approx_mean
        elif stat == "mod_approx_mean":
            return self.chromosome.get_mod_approx_mean
        elif stat == "max":
            return self.chromosome.get_max
        elif stat == "min":
            return self.chromosome.get_min
        elif stat == "coverage":
            return self.chromosome.get_coverage
        elif stat == "std":
            return self.chromosome.get_std
        else:
            print(f"{stat} is not a valid statistic to search for")
            return None

    def stats(self, intervals, stat="mean"):

        if len(intervals[0]) != 3:
            print(f"List given has incorrect formatting. It must be in the format:\n"
                  "[[chr1, 1, 100], [chr1, 101, 200], ...]")
            return None

        method = self.get_method(stat)
        results = []

        for interval in intervals:
            results.append(method(interval[1], interval[2]))

        return results

    def stats_from_file(self, interval_file, output_file=None, stat="mean"):
        output = ""
        method = self.get_method(stat)
        with open(interval_file) as in_file:
            for line in in_file:
                interval = line.split()

                if len(interval) != 3:
                    print(f"{interval_file} has incorrect formatting. It must be in the format:\n"
                          "chr1\t100\t401\n"
                          "chr1\t600\t1000\n"
                          "...")
                    break

                result = method(int(interval[1]), int(interval[2]))

                if output_file is None:
                    print(interval[0], interval[1], interval[2], result)
                else:
                    output += f"{result}\n"

        if output_file:
            with open(output_file, 'w') as out_file:
                out_file.write(output)
