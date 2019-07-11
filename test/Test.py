import time
import sys
import math
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
from pyBedGraph.Benchmark import Benchmark

if len(sys.argv) < 5:
    print("Needs 4 arguments:\n"
          "arg 1 - chrom_sizes_file\n"
          "arg 2 - bedGraph_file\n"
          "arg 3 - bigWig_file\n"
          "arg 4 - interval file")
    exit(-1)

MAX_NUM_TESTS = 1000000
test_intervals = []
average_interval_size = 500
num_tests = 100000
chrom_name = 'chr1'
bin_size = int(math.sqrt(average_interval_size))
stats = ['mean']

bedGraph = BedGraph(sys.argv[1], sys.argv[2], chrom_name)
bedGraph.load_chrom_data(chrom_name)
bins = [x for x in range(10, 101, 1)]
with open('out/values_indexed.txt', 'w') as out_file:
    out_file.write(' '.join([str(x) for x in bins]) + '\n')
    for bin_size in bins:
        bedGraph.load_chrom_bins(chrom_name, bin_size)
        bench = Benchmark(bedGraph, sys.argv[3])
        result, bins = bench.benchmark(num_tests, average_interval_size, chrom_name, bin_size, stats,
                                 bench_pyBigWig=False, pyBigWig_baseline=False, only_runtime=True)
        out_file.write(str(bins) + " ")
        for key in result:
            print(key, result[key])

exit()
'''

complete_bedGraph = BedGraph(sys.argv[1], sys.argv[2], chrom_name,
                             ignore_missing_bp=False)
complete_bedGraph.load_chrom_data(chrom_name)
complete_bedGraph.load_chrom_bins(chrom_name, bin_size)


with open(sys.argv[4]) as interval_file:
    count = 0
    for line in interval_file:
        if count > MAX_NUM_TESTS:
            break
        test_intervals.append(line.split())
        count += 1

for stat in stats:
    start_time = time.time()
    stat_values = bedGraph.stats(stat, test_intervals)
    print("Time for bedGraph stats:", time.time() - start_time)

    start_time = time.time()
    bedGraph.stats_from_file(sys.argv[4], stat=stat)
    print("Time for bedGraph stats_from_file to output file:", time.time() - start_time)

    start_time = time.time()
    file_values = bedGraph.stats_from_file(sys.argv[4], output_to_file=False, stat=stat)
    print("Time for bedGraph stats_from_file to return:", time.time() - start_time)

    file_values = file_values[chrom_name]
    try:
        assert len(file_values) == len(stat_values)
    except AssertionError:
        print(len(file_values), len(stat_values))
        exit(-1)
    size = len(file_values)
    with open('chr1_out.txt') as in_file:
        for i in range(size):
            file_value = float(in_file.readline())
            try:
                assert stat_values[i] == file_values[i]
                assert file_value == file_values[i]
            except AssertionError:
                print(stat_values[i], file_values[i])
                print(file_value, file_values[i])
                exit(-1)

    start_time = time.time()
    stat_values = complete_bedGraph.stats(stat, test_intervals)
    print("Time for complete_bedGraph stats:", time.time() - start_time)

    start_time = time.time()
    complete_bedGraph.stats_from_file(sys.argv[4], stat=stat)
    print("Time for complete_bedGraph stats_from_file to output file:", time.time() - start_time)

    start_time = time.time()
    file_values = complete_bedGraph.stats_from_file(sys.argv[4], output_to_file=False, stat=stat)
    print("Time for complete_bedGraph stats_from_file to return:", time.time() - start_time)

    file_values = file_values[chrom_name]
    try:
        assert len(file_values) == len(stat_values)
    except AssertionError:
        print(len(file_values), len(stat_values))
        exit(-1)
    size = len(file_values)
    with open('chr1_out.txt') as in_file:
        for i in range(size):
            file_value = float(in_file.readline())
            try:
                assert stat_values[i] == file_values[i]
                assert file_value == file_values[i]
            except AssertionError:
                print(stat_values[i], file_values[i])
                print(file_value, file_values[i])
                exit(-1)

a = input('end program')
complete_bedGraph.free_chrom_data(chrom_name)
bedGraph.free_chrom_data(chrom_name)
a = input('end program')'''
