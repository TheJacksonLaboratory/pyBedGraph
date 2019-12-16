from .include_missing_bp import *
from .ignore_missing_bp import *
from .Chrom_Data import Chrom_Data
from .Chrom_Data_Complete import Chrom_Data_Complete
from .BedGraph import BedGraph

'''
NUM_TESTS = 10000

MIN_INTERVAL_SIZE = 100
MAX_INTERVAL_SIZE = 1000
INTERVAL_STEP = 300

MIN_BIN_SIZE = 50
MAX_BIN_SIZE = 500
BIN_STEP = 100
# List of arguments (required)
# arg1 - bedgraph file
# arg2 - name of chromosome
# arg3 - 'b' for benchmark, 'r' to run

# arg4 - file containing intervals to search for
# or
# arg4 - bigWig file to benchmark pyBigWig
if len(sys.argv) != 3:
    print("List of required arguments:\n"
          "arg1 - bedgraph file\n"
          "arg2 - name of chromosome to search in\n"
          "arg4 - bigWig file to benchmark pyBigWig | file containing intervals to search for\n")
    exit(-1)

genome = Genome(sys.argv[1], sys.argv[2])

bin_size = 512
while True:

    action = input("Enter 'b' for benchmarking pyBigWig | 's' to search for intervals: ")

    if action == 'b':
        bigWig_file = input("Enter a file to benchmark pyBigWig: ")
        Benchmark(genome.chromosomes[sys.argv[2]], NUM_TESTS, bin_size, 500, sys.argv[4], 'exact')

    elif action == 's':
        interval_file = input("Enter a file that contains intervals to search for: ")
        stat = input("Enter the statistic to search for")
        with open(interval_file) as in_file:
            for line in in_file:
                interval = line.split("")

                if len(interval) != 2:
                    print(f"{sys.argv[4]} has incorrect formatting. It must be in the format:\n"
                          "start_of_interval1 end_of_interval1\n"
                          "start_of_interval2 end_of_interval2\n"
                          "...")
                    break

                genome.stats()

    else:
        print(f"{action} is not a valid action")

    if input('enter y to run again: ') != 'y':
        break

    bin_size = int(input('enter maximum bin size: '))
'''
