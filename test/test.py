import time
import sys
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
from pyBedGraph.Benchmark import Benchmark

if len(sys.argv) != 5:
    print("Needs 3 arguments:\n"
          "arg 1 - chrom_sizes_file\n"
          "arg 2 - bedgraph_file\n"
          "arg 3 - bigWig file\n"
          "arg 4 - chromosome name")
    exit(-1)

start_time = time.time()
bedGraph = BedGraph(sys.argv[1], sys.argv[2], sys.argv[4])
print("Time for loading bedGraph file: ", time.time() - start_time)

bench = Benchmark(bedGraph, sys.argv[3])

num_tests = 2500
interval_size = 1000
chrom_name = sys.argv[4]
bin_size = 500
stats = ['mean', 'approx_mean', 'mod_approx_mean', 'mod2_approx_mean']
results = bench.benchmark(num_tests, interval_size, chrom_name, bin_size, stats)

for stat in results:
    print(stat, results[stat])
    if 'run_time' in results[stat]:
        print(results[stat]['run_time'] / results['pyBigWig_mean']['exact_run_time'])