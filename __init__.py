from Data import Data
from Benchmark import benchmark
import pyBigWig
import time
import random

NUM_TESTS = 5000
NUM_INDEX = 23011544 # largest allele_index in CHROM_NAME for test
MAX_WIDTH = 100000
CHROM_NAME = "chr2L"

MIN_INTERVAL_SIZE = 100
MAX_INTERVAL_SIZE = 1000
INTERVAL_STEP = 100

MIN_BIN_SIZE = 50
MAX_BIN_SIZE = 1000
BIN_STEP = 50

test_cases = []
actual = []

# create test cases
# test[0]: start of interval
# test[1]: end of interval
for i in range(NUM_TESTS):
    test_case = {}
    test_case[0] = random.randint(1, NUM_INDEX - MAX_WIDTH)
    test_case[1] = test_case[0] + MAX_WIDTH
    # test[1] = random.randint(start + 1, start + MAX_WIDTH)

    test_cases.append(test_case)

################################################################################

# baseline
bw = pyBigWig.open('data/P2MC7N8HCE3K.bw')
start_time = time.time()
for i in range(NUM_TESTS):
    value = bw.stats(CHROM_NAME, test_cases[i][0], test_cases[i][1])
    actual.append(value[0])
bigWig_time = time.time() - start_time
print(f"pyBigWig: {bigWig_time} seconds for {NUM_TESTS} trials")

################################################################################

obj = Data('data/P2MC7N8HCE3K.bedgraph', 500)
print()

for interval_size in range(MIN_INTERVAL_SIZE, MAX_INTERVAL_SIZE + 1, INTERVAL_STEP):
    for bin_size in range(MIN_BIN_SIZE, MAX_BIN_SIZE + 1, BIN_STEP):
        benchmark(obj, NUM_TESTS, bin_size, interval_size, bigWig_time, test_cases, actual)
