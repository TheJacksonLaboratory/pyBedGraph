import sys
import os
import time
from pathlib import Path
import math
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
from pyBedGraph.Benchmark import Benchmark
import generate_images

MIN_NUM_TEST = 100000
MAX_NUM_TEST = 1000000
DEFAULT_INTERVAL_SIZE = 500
DEFAULT_NUM_TESTS = 10000

total_start_time = time.time()


def interval_size_error_benchmark():
    interval_error_results = {}
    for name in generate_images.INTERVAL_ERROR_NAMES:
        interval_error_results[name] = []

    stats_to_bench = ['mean']
    for interval_size in interval_test_list:
        bin_size = int(math.sqrt(interval_size))
        result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size, chrom_name,
                                 bin_size, stats_to_bench)
        interval_error_results['pyBigWig_approx'].append(result['pyBigWig_mean']['error'])

        print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    # stats_to_bench = ['approx_mean', 'mod_approx_mean']
    stats_to_bench = ['approx_mean']
    for bin_size_divide in bin_size_test_list:
        interval_error_results['pyBedGraph_approx_' + str(bin_size_divide)] = []
        # interval_error_results['pyBedGraph_mod_approx_' + str(bin_size)] = []
        for interval_size in interval_test_list:
            bin_size = int(interval_size / bin_size_divide)
            result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size,
                                     chrom_name, bin_size, stats_to_bench,
                                     bench_pyBigWig=False)
            interval_error_results['pyBedGraph_approx_' + str(bin_size_divide)].append(result['approx_mean']['error'])
            # interval_error_results['pyBedGraph_mod_approx_' + str(bin_size)].append(result['mod_approx_mean']['error'])

            print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    print(interval_error_results)
    with open(f'graphs/{data_name}/interval_error_results.txt', 'w') as out:
        out.write(" ".join([str(x) for x in interval_test_list]) + '\n')
        for key in interval_error_results:
            output = key + " " + " ".join([str(x) for x in interval_error_results[key]]) + '\n'
            out.write(output)

    # generate_images.create_error_interval(data_name, interval_test_list, interval_error_results)


def interval_size_runtime_benchmark():
    interval_runtime_results = {}
    for name in generate_images.INTERVAL_RUNTIME_NAMES:
        interval_runtime_results[name] = []

    stats_to_bench = ['mean']
    for interval_size in interval_test_list:
        bin_size = int(math.sqrt(interval_size))
        result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size, chrom_name,
                                 bin_size, stats_to_bench, True)
        interval_runtime_results['pyBigWig_approx'].append(result['pyBigWig_mean']['approx_run_time'])
        interval_runtime_results['pyBigWig_exact'].append(result['pyBigWig_mean']['exact_run_time'])
        interval_runtime_results['pyBedGraph_exact'].append(result['mean']['run_time'])

        print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    # stats_to_bench = ['approx_mean', 'mod_approx_mean']
    stats_to_bench = ['approx_mean']
    for bin_size_divide in bin_size_test_list:
        interval_runtime_results['pyBedGraph_approx_' + str(bin_size_divide)] = []
        # interval_error_results['pyBedGraph_mod_approx_' + str(bin_size)] = []
        for interval_size in interval_test_list:
            bin_size = int(interval_size / bin_size_divide)
            result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size,
                                     chrom_name, bin_size, stats_to_bench,
                                     True, False)
            interval_runtime_results['pyBedGraph_approx_' + str(bin_size_divide)].append(result['approx_mean']['run_time'])
            # interval_error_results['pyBedGraph_mod_approx_' + str(bin_size)].append(result['mod_approx_mean']['error'])

            print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    print(interval_runtime_results)
    with open(f'graphs/{data_name}/interval_runtime_results.txt', 'w') as out:
        out.write(" ".join([str(x) for x in interval_test_list]) + '\n')
        for key in interval_runtime_results:
            output = key + " " + " ".join([str(x) for x in interval_runtime_results[key]]) + '\n'
            out.write(output)


def runtime_benchmark():
    # create list of num_tests
    num_test_list = [x for x in range(MIN_NUM_TEST, MAX_NUM_TEST + 1, MIN_NUM_TEST)]

    stats_to_bench = ['mean']
    run_time_results = {}
    for name in generate_images.RUN_TIME_NAMES:
        run_time_results[name] = []

    bin_size = int(math.sqrt(DEFAULT_INTERVAL_SIZE))
    for num_test in num_test_list:
        result = bench.benchmark(num_test, DEFAULT_INTERVAL_SIZE, chrom_name,
                                 bin_size, stats_to_bench, True, True)
        run_time_results['pyBedGraph_exact'].append(result['mean']['run_time'])
        run_time_results['pyBigWig_approx'].append(result['pyBigWig_mean']['approx_run_time'])
        run_time_results['pyBigWig_exact'].append(result['pyBigWig_mean']['exact_run_time'])

        print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    stats_to_bench = ['approx_mean']
    for bin_size_divide in bin_size_test_list:
        run_time_results['pyBedGraph_approx_' + str(bin_size_divide)] = []
        bin_size = int(DEFAULT_INTERVAL_SIZE / bin_size_divide)
        # run_time_results['pyBedGraph_mod_approx_' + str(bin_size)] = []
        for num_test in num_test_list:
            result = bench.benchmark(num_test, DEFAULT_INTERVAL_SIZE,
                                     chrom_name, bin_size, stats_to_bench, True,
                                     False, False)
            run_time_results['pyBedGraph_approx_' + str(bin_size_divide)].append(result['approx_mean']['run_time'])
            # run_time_results['pyBedGraph_mod_approx_' + str(bin_size)].append(result['mod_approx_mean']['run_time'])

            print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    print(run_time_results)
    with open(f'graphs/{data_name}/run_time_results.txt', 'w') as out:
        out.write(" ".join([str(x) for x in num_test_list]) + '\n')
        for key in run_time_results:
            output = key + " " + " ".join([str(x) for x in run_time_results[key]]) + '\n'
            out.write(output)

    # generate_images.create_runtime_num_test(data_name, num_test_list, run_time_results)


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

# bin_size_test_list = test_numbers[Path(sys.argv[1]).stem + '.sizes']['bin_sizes']
interval_test_list = test_numbers[Path(sys.argv[1]).stem + '.sizes']['intervals']
bin_size_test_list = [2, 10, 20]

runtime_benchmark()
#interval_size_error_benchmark()
#interval_size_runtime_benchmark()
