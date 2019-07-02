import sys
import os
import time
from pathlib import Path
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
from pyBedGraph.Benchmark import Benchmark
import generate_images

MIN_NUM_TEST = 125
MAX_NUM_TEST = 10000


def interval_size_error_benchmark():
    num_tests = 5000

    interval_error_results = {}
    for name in generate_images.INTERVAL_ERROR_NAMES:
        interval_error_results[name] = []

    for interval_size in interval_test_list:
        result = bench.benchmark(num_tests, interval_size, chrom_name, None, ['mean'])
        interval_error_results['pyBigWig_approx'].append(result['pyBigWig_mean']['error'])

    stats_to_bench = ['approx_mean', 'mod_approx_mean']
    for bin_size in bin_size_test_list:
        interval_error_results['pyBedGraph_approx_' + str(bin_size)] = []
        interval_error_results['pyBedGraph_mod_approx_' + str(bin_size)] = []
        for interval_size in interval_test_list:
            result = bench.benchmark(num_tests, interval_size, chrom_name, bin_size, stats_to_bench)
            interval_error_results['pyBedGraph_approx_' + str(bin_size)].append(result['approx_mean']['error'])
            interval_error_results['pyBedGraph_mod_approx_' + str(bin_size)].append(result['mod_approx_mean']['error'])

    print(interval_error_results)

    generate_images.create_error_interval(data_name, interval_test_list, interval_error_results)


def runtime_benchmark():
    # create list of num_tests
    num_test_list = []
    num_test = MIN_NUM_TEST
    while num_test <= MAX_NUM_TEST:
        num_test_list.append(num_test)
        diff = 1000
        if num_test > diff:
            num_test += diff
        else:
            num_test *= 2

    stats_to_bench = ['mean']
    run_time_results = {}
    for name in generate_images.RUN_TIME_NAMES:
        run_time_results[name] = []

    for num_test in num_test_list:
        result = bench.benchmark(num_test, 500, chrom_name, None, stats_to_bench)
        run_time_results['pyBedGraph_exact'].append(result['mean']['run_time'])
        run_time_results['pyBigWig_exact'].append(result['pyBigWig_mean']['exact_run_time'])
        run_time_results['pyBigWig_approx'].append(result['pyBigWig_mean']['approx_run_time'])

    stats_to_bench = ['approx_mean', 'mod_approx_mean']
    for bin_size in bin_size_test_list:
        run_time_results['pyBedGraph_approx_' + str(bin_size)] = []
        run_time_results['pyBedGraph_mod_approx_' + str(bin_size)] = []
        for num_test in num_test_list:
            result = bench.benchmark(num_test, 500, chrom_name, bin_size, stats_to_bench)
            run_time_results['pyBedGraph_approx_' + str(bin_size)].append(result['approx_mean']['run_time'])
            run_time_results['pyBedGraph_mod_approx_' + str(bin_size)].append(result['mod_approx_mean']['run_time'])

    print(run_time_results)
    generate_images.create_runtime_num_test(data_name, num_test_list, run_time_results)


if len(sys.argv) != 5:
    print("Needs 3 arguments:\n"
          "arg 1 - chrom_sizes_file\n"
          "arg 2 - bedgraph_file\n"
          "arg 3 - bigWig file\n"
          "arg 4 - chrom_name")
    exit(-1)

chrom_name = sys.argv[4]

start_time = time.time()
bedGraph = BedGraph(sys.argv[1], sys.argv[2], chrom_name)
print("Time for loading bedGraph file: ", time.time() - start_time)

start_time = time.time()
print(f"Time for loading {sys.argv[4]}: ", time.time() - start_time, '\n')

bench = Benchmark(bedGraph, sys.argv[3])

data_name = Path(sys.argv[2]).stem
if not os.path.isdir(f'graphs'):
    os.mkdir(f'graphs')
if not os.path.isdir(f'graphs/{data_name}'):
    os.mkdir(f'graphs/{data_name}')

test_numbers = {}
with open('test_numbers.txt') as in_file:
    for i in range(3):
        line = in_file.readline()
        chrom_size_file = line.strip()
        test_numbers[chrom_size_file] = {}

        intervals = in_file.readline()
        bin_sizes = in_file.readline()
        intervals = [int(x) for x in intervals.split()]
        bin_sizes = [int(x) for x in bin_sizes.split()]

        test_numbers[chrom_size_file]['intervals'] = intervals
        test_numbers[chrom_size_file]['bin_sizes'] = bin_sizes

bin_size_test_list = test_numbers[Path(sys.argv[1]).stem + '.sizes']['bin_sizes']
interval_test_list = test_numbers[Path(sys.argv[1]).stem + '.sizes']['intervals']

interval_size_error_benchmark()
# runtime_benchmark()
