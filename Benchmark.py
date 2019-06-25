from Data import Data
import time

CHROM_NAME = "chr2L"


def benchmark(obj, num_tests, bin_size, interval_size, bigWig_time, test_cases, actual):

    print(f"Benchmark for interval size: {interval_size} | bin size: {bin_size}")

    obj.split_bins(bin_size)

    predicted = []

    start_time = time.time()
    for i in range(num_tests):
        value = obj.get_average(CHROM_NAME, test_cases[i][0], test_cases[i][1])
        predicted.append(value)
    num_sec = time.time() - start_time

    print(f"myTest: {num_sec} seconds for {num_tests} trials")

    ################################################################################

    error = 0
    num = 0
    for i in range(num_tests):
        if actual[i] is not None and predicted[i] is not None:
            error += (actual[i] - predicted[i]) * (actual[i] - predicted[i])
            num += 1

    print(f"My test takes {num_sec / bigWig_time} of bigWig's time")
    print(f"Error: {error / num}")
    print(f"Num: {num}")
    print()
