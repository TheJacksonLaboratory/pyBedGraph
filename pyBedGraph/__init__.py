from Genome import Genome
from Benchmark import Benchmark
import sys

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

if len(sys.argv) != 5:
    print("List of required arguments:\n"
          "arg1 - bedgraph file\n"
          "arg2 - name of chromosome\n"
          "arg3 - 'b' for benchmarking pyBigWig | 's' to search for intervals\n"
          "arg4 - bigWig file to benchmark pyBigWig | file containing intervals to search for\n")
    exit(-1)

if sys.argv[3] != 'b' and sys.argv[3] != 's':
    print(f"{sys.argv[3]} is not a valid input. It must be either 'b' or 's'")

genome = Genome(sys.argv[1], sys.argv[2])

bin_size = 512
while True:
    if sys.argv[3] == 'b':
        Benchmark(genome.chromosomes[sys.argv[2]], NUM_TESTS, bin_size, 500, sys.argv[4], 'exact')
    elif sys.argv[3] == 's':
        print("Have not implemented searching for intervals yet")

    if input('enter y to run again: ') != 'y':
        break

    bin_size = int(input('enter maximum bin size: '))
