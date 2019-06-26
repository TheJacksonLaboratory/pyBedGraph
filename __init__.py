from Data import Data
from Benchmark import benchmark

NUM_TESTS = 50000

MIN_INTERVAL_SIZE = 100
MAX_INTERVAL_SIZE = 1000
INTERVAL_STEP = 300

MIN_BIN_SIZE = 50
MAX_BIN_SIZE = 500
BIN_STEP = 100
"""
obj = Data('data/human.bedgraph', 1, 'chr1')
print()

"""
obj = Data('data/P2MC7N8HCE3K_chr2L.bedgraph', "chr2L")
print()

benchmark(obj, NUM_TESTS, 100, 500)
"""
for interval_size in range(MIN_INTERVAL_SIZE, MAX_INTERVAL_SIZE + 1, INTERVAL_STEP):
    for bin_size in range(MIN_BIN_SIZE, MAX_BIN_SIZE + 1, BIN_STEP):
        if bin_size > interval_size:
            continue
        benchmark(obj, NUM_TESTS, bin_size, interval_size)
        """
