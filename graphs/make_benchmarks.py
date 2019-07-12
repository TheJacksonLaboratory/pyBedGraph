import sys
import os
import time
from pathlib import Path
import math
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
from pyBedGraph.Benchmark import Benchmark
import generate_graphs

MIN_NUM_TEST = 100000
MAX_NUM_TEST = 1000000
DEFAULT_INTERVAL_SIZE = 500
DEFAULT_NUM_TESTS = 10000

total_start_time = time.time()

interval_test_list = [100, 250, 500, 750, 1000, 2000, 3000, 4000, 5000]
bin_size_test_list = [5, 10, 20]


def interval_size_error_benchmark():
    interval_error_results = {}
    for name in generate_graphs.INTERVAL_ERROR_NAMES:
        interval_error_results[name] = []

    stats_to_bench = ['mean']
    for interval_size in interval_test_list:
        result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size, chrom_name,
                                 None, stats_to_bench)
        interval_error_results['pyBW app.'].append(result['pyBigWig_mean']['error'])

        print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    stats_to_bench = ['approx_mean']
    for bin_size_divide in bin_size_test_list:
        interval_error_results['pyBG app. bin=int_size/' + str(bin_size_divide)] = []
        for interval_size in interval_test_list:
            bin_size = int(interval_size / bin_size_divide)
            result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size,
                                     chrom_name, bin_size, stats_to_bench,
                                     bench_pyBigWig_approx=False)
            interval_error_results['pyBG app. bin=int_size/' + str(bin_size_divide)].append(result['approx_mean']['error'])

            print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    print(interval_error_results)
    error_types = ['percent_error', 'ms_error', 'abs_error', 'num_actual_0']
    with open(f'graphs/{data_name}/interval_error_results.txt', 'a+') as out:
        out.write(" ".join([str(x) for x in interval_test_list]) + '\n')
        for key in interval_error_results:
            out.write(key + '\n')
            error_list = interval_error_results[key]
            for error_type in error_types:
                out.write(error_type + " ")
                for error_dict in error_list:
                    out.write(str(error_dict[error_type]) + " ")
                out.write('\n')

    # generate_images.create_error_interval(data_name, interval_test_list, interval_error_results)


def interval_size_runtime_benchmark():
    interval_runtime_results = {}
    for name in generate_graphs.INTERVAL_RUNTIME_NAMES:
        interval_runtime_results[name] = []

    stats_to_bench = ['mean']
    for interval_size in interval_test_list:
        result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size, chrom_name,
                                 None, stats_to_bench, True, True)
        interval_runtime_results['pyBW app.'].append(result['pyBigWig_mean']['approx_run_time'])
        interval_runtime_results['pyBW exact'].append(result['pyBigWig_mean']['exact_run_time'])
        interval_runtime_results['pyBG exact'].append(result['mean']['run_time'])

        print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    stats_to_bench = ['approx_mean']
    for bin_size_divide in bin_size_test_list:
        interval_runtime_results['pyBG app. bin=int_size/' + str(bin_size_divide)] = []
        for interval_size in interval_test_list:
            bin_size = int(interval_size / bin_size_divide)
            result = bench.benchmark(DEFAULT_NUM_TESTS, interval_size,
                                     chrom_name, bin_size, stats_to_bench,
                                     True, False, False)
            interval_runtime_results['pyBG app. bin=int_size/' + str(bin_size_divide)].append(result['approx_mean']['run_time'])

            print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    print(interval_runtime_results)
    with open(f'graphs/{data_name}/interval_runtime_results.txt', 'a+') as out:
        out.write(" ".join([str(x) for x in interval_test_list]) + '\n')
        for key in interval_runtime_results:
            output = key + "\n" + " ".join([str(x) for x in interval_runtime_results[key]]) + '\n'
            out.write(output)


def runtime_benchmark():
    # create list of num_tests
    num_test_list = [x for x in range(MIN_NUM_TEST, MAX_NUM_TEST + 1, MIN_NUM_TEST)]

    run_time_results = {}

    stats_to_bench = ['mean']
    for name in generate_graphs.RUN_TIME_NAMES:
        run_time_results[name] = []

    for num_test in num_test_list:
        result = bench.benchmark(num_test, DEFAULT_INTERVAL_SIZE, chrom_name,
                                 None, stats_to_bench, True, True)
        run_time_results['pyBG exact'].append(result['mean']['run_time'])
        run_time_results['pyBW app.'].append(result['pyBigWig_mean']['approx_run_time'])
        run_time_results['pyBW exact'].append(result['pyBigWig_mean']['exact_run_time'])

        print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    stats_to_bench = ['approx_mean']
    for bin_size_divide in bin_size_test_list:
        bin_size = int(DEFAULT_INTERVAL_SIZE / bin_size_divide)
        run_time_results['pyBG app. bin=' + str(bin_size)] = []
        for num_test in num_test_list:
            result = bench.benchmark(num_test, DEFAULT_INTERVAL_SIZE,
                                     chrom_name, bin_size, stats_to_bench, True,
                                     False, False)
            run_time_results['pyBG app. bin=' + str(bin_size)].append(result['approx_mean']['run_time'])

            print(f"Total time taken so far (min): {(time.time() - total_start_time) / 60}")

    print(run_time_results)
    with open(f'graphs/{data_name}/run_time_results.txt', 'a') as out:
        out.write(" ".join([str(x) for x in num_test_list]) + '\n')
        for key in run_time_results:
            output = key + "\n" + " ".join([str(x) for x in run_time_results[key]]) + '\n'
            out.write(output)

    # generate_images.create_runtime_num_test(data_name, num_test_list, run_time_results)


if len(sys.argv) != 4:
    print("Needs 3 arguments:\n"
          "arg 1 - chrom_sizes_file\n"
          "arg 2 - bedgraph_file\n"
          "arg 3 - bigWig file")
    exit(-1)

chrom_name = 'chr1'

start_time = time.time()
bedGraph = BedGraph(sys.argv[1], sys.argv[2], chrom_name)
print("Time for loading bedGraph file: ", time.time() - start_time)

start_time = time.time()
print(f"Time for loading {chrom_name}: ", time.time() - start_time, '\n')

bench = Benchmark(bedGraph, sys.argv[3])

data_name = Path(sys.argv[2]).stem
if not os.path.isdir(f'graphs'):
    os.mkdir(f'graphs')
if not os.path.isdir(f'graphs/{data_name}'):
    os.mkdir(f'graphs/{data_name}')

runtime_benchmark()
interval_size_error_benchmark()
interval_size_runtime_benchmark()
