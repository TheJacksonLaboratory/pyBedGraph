import random
import pyBigWig
import time

NUM_INDEX = 23011544  # largest allele_index in CHROM_NAME for test
CHROM_NAME = "chr2L"


# Make sure # of bins is small enough to avoid calling np.mean too many times
def benchmark(obj, num_tests, bin_size, interval_size):
    test_cases = create_test_cases(num_tests, interval_size)
    bigWig_time, actual = find_baseline(test_cases)

    print(f"Benchmark for interval size: {interval_size} | bin size: {bin_size}")

    chrom = obj.chromosomes[CHROM_NAME]
    chrom.split_bins(bin_size)

    predicted = []

    start_time = time.time()
    for i in range(num_tests):
        value = obj.get_approx_average(CHROM_NAME, test_cases[i][0], test_cases[i][1])
        #value = obj.get_exact_average(CHROM_NAME, test_cases[i][0], test_cases[i][1])
        predicted.append(value)
    num_sec = time.time() - start_time

    print(f"myTest: {num_sec} seconds for {num_tests} trials")

    ################################################################################

    error = 0
    num = 0
    for i in range(num_tests):
        if actual[i] is None:
            continue

        if predicted[i] is None or not predicted[i] >= 0:
            # print("Error with predicted or actual", test_cases[i][0], test_cases[i][1],
            #     predicted[i], actual[i])
            continue
        # if (actual[i] - predicted[i]) * (actual[i] - predicted[i]) > 0.0001:
        #    print(test_cases[i][0], test_cases[i][1], actual[i], predicted[i])

        error += (actual[i] - predicted[i]) * (actual[i] - predicted[i])
        num += 1

    print(f"My test takes {num_sec / bigWig_time} of bigWig's time")
    print(f"Error: {error / num}")
    print(f"Num: {num}")
    print()


def create_test_cases(num_tests, interval_size):
    test_cases = []

    # test[0]: start of interval
    # test[1]: end of interval
    for i in range(num_tests):
        test_case = {}
        test_case[0] = random.randint(1, NUM_INDEX - interval_size)
        # test_case[0] = 8591040
        test_case[1] = test_case[0] + interval_size
        # test[1] = random.randint(start + 1, start + MAX_WIDTH)

        test_cases.append(test_case)

    return test_cases


def find_baseline(test_cases):
    num_tests = len(test_cases)
    actual = []

    bw = pyBigWig.open('data/P2MC7N8HCE3K.bw')
    start_time = time.time()
    for i in range(num_tests):
        value = bw.stats(CHROM_NAME, test_cases[i][0], test_cases[i][1])
        actual.append(value[0])
    bigWig_time = time.time() - start_time
    print(f"pyBigWig: {bigWig_time} seconds for {num_tests} trials")

    return bigWig_time, actual
