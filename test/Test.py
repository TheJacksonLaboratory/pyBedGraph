import time
import sys
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
from pyBedGraph.Benchmark import ALL_STATS, Benchmark
import pyBigWig

TEST_COUNTER = 1000000

if len(sys.argv) != 6:
    print("Needs 5 arguments:\n"
          "arg 1 - chrom_sizes_file\n"
          "arg 2 - bedgraph_file\n"
          "arg 3 - interval file\n"
          "arg 4 - out file\n",
          "arg 5 - bigWig_file")
    exit(-1)

bin_size = 200
test_intervals = []
num_tests = 1000000
interval_size = 500
chrom_name = 'chr1'
stats = ['mean', 'approx_mean', 'mod_approx_mean']

bedGraph = BedGraph(sys.argv[1], sys.argv[2], chrom_name)
bench = Benchmark(bedGraph, sys.argv[5])
result = bench.benchmark(num_tests, interval_size, chrom_name, bin_size, stats,
                         bench_pyBigWig=False)
for key in result:
    print(key, result[key])

exit()

complete_bedGraph = BedGraph(sys.argv[1], sys.argv[2], chrom_name,
                             like_pyBigWig=False)
complete_bedGraph.load_chrom_data(chrom_name)


with open(sys.argv[3]) as interval_file:
    count = 0
    for line in interval_file:
        if count > num_tests:
            break
        test_intervals.append(line.split())
        count += 1

while True:
    stat = input('Enter a stat name: ')

    print("Testing", stat)

    start_time = time.time()
    bedGraph.stats(test_intervals, stat)
    print("Time for bedGraph stats:", time.time() - start_time)

    start_time = time.time()
    complete_bedGraph.stats(test_intervals, stat)
    print("Time for complete_bedGraph stats:", time.time() - start_time)

    bw = pyBigWig.open(sys.argv[5])
    start_time = time.time()
    for test_interval in test_intervals:
        bw.stats(chrom_name, int(test_interval[1]), int(test_interval[2]),
                 type=stat, exact=False)
    print("Time for pyBigWig stats:", time.time() - start_time)
