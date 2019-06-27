from Chromosome import Chromosome

CHROM_NAME_INDEX = 0


class Genome:

    def __init__(self, data_file_name, chromosome_name=None, bin_size=64):

        self.chromosome = None

        print(f"Reading in {data_file_name}...")

        with open(data_file_name, 'r') as data_file:

            for line in data_file:
                data = line.split()

                if chromosome_name is not None and chromosome_name != data[CHROM_NAME_INDEX]:

                    # have not yet created specified chromosome and added data to it
                    if self.chromosome is None:
                        continue

                    # done with adding data, stop now
                    else:
                        break

                # update self.chromosome to be the chromosome to add data to
                if self.chromosome is None or data[CHROM_NAME_INDEX] != self.chromosome.name:
                    self.chromosome = Chromosome(data[CHROM_NAME_INDEX])

                self.chromosome.add_data(data)

            if self.chromosome is None:
                print(f"{chromosome_name} was not found in {data_file_name}")
                exit(-1)

        print("Finished\n")

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

        method = self.get_method(stat)
        results = []

        for interval in intervals:
            results.append(method(interval[0], interval[1]))

        return results

    def stats_from_file(self, interval_file, output_file=None, stat="mean"):
        output = ""
        method = self.get_method(stat)
        with open(interval_file) as in_file:
            for line in in_file:
                interval = line.split()

                if len(interval) != 2:
                    print(f"{interval_file} has incorrect formatting. It must be in the format:\n"
                          "start_of_interval1\tend_of_interval1\n"
                          "start_of_interval2\tend_of_interval2\n"
                          "...")
                    break

                result = method(int(interval[0]), int(interval[1]))

                if output_file is None:
                    print(interval[0], interval[1], result)
                else:
                    output += f"{result}\n"

        if output_file:
            with open(output_file, 'w') as out_file:
                out_file.write(output)
