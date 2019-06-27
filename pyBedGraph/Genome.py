from Chromosome import Chromosome

CHROM_NAME_INDEX = 0


class Genome:

    def __init__(self, data_file_name, chromosome_name=None):

        current = None

        self.chromosomes = {}

        print(f"Reading in {data_file_name}...")

        with open(data_file_name, 'r') as data_file:

            for line in data_file:
                data = line.split("\t")

                if chromosome_name is not None and chromosome_name != data[CHROM_NAME_INDEX]:

                    # have not yet created specified chromosome and added data to it
                    if current is None:
                        continue

                    # done with adding data, stop now
                    else:
                        break

                # update current to be the chromosome to add data to
                if current is None or data[CHROM_NAME_INDEX] != current.name:
                    current = Chromosome(data[CHROM_NAME_INDEX])
                    self.chromosomes[current.name] = current

                current.add_data(data)

            if current is None:
                print(f"{chromosome_name} was not found in {data_file_name}")
                exit(-1)

        print("Finished\n")

