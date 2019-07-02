import random
import pyBigWig
import time
import numpy as np

EXACT_MEAN_INDEX = 0
APPROX_MEAN_INDEX = 1
MOD_APPROX_MEAN_INDEX = 2
RANDOM_SEED = 1

mean_names = [
    'pyBigWig_exact',
    'pyBigWig_approx',
    'pyBedGraph_exact',
    'pyBedGraph_approx',
    'pyBedGraph_mod_approx'
]

ALL_STATS = [
    "mean",
    "approx_mean",
    "mod_approx_mean",
    "max",
    "min",
    "coverage",
    "std"
]



class Benchmark:

    def __init__(self, bedGraph, bigWig_file):

        self.bw = pyBigWig.open(bigWig_file)
        self.bedGraph = bedGraph

        self.intervals_list = []

        self.num_tests = 0
        self.test_cases = []

        self.chromosome = None

    def benchmark(self, num_tests, interval_size, chrom_name, bin_size, stats=None):

        print("Benchmarking:\n"
              f"Number of tests: {num_tests}\n"
              f"Interval size: {interval_size}\n"
              f"Chromosome name: {chrom_name}\n"
              f"Bin size: {bin_size}\n")

        # self.find_intervals()

        self.bedGraph.load_chrom_data(chrom_name)
        self.chromosome = self.bedGraph.chromosome_map[chrom_name]
        if bin_size is not None:
            self.chromosome.split_bins(bin_size)

        self.create_test_cases(num_tests, interval_size)

        # benchmark all stats if none given
        if stats is None:
            stats = ALL_STATS

        results = {}
        actual = {}
        predictions = {}

        for stat_name in stats:

            results[stat_name] = {}

            # actual value for approx_mean/mod_approx_mean is mean
            actual_stat_name = stat_name
            if 'mean' in stat_name and stat_name != 'mean':
                actual_stat_name = 'mean'
            pyBigWig_name = 'pyBigWig_' + actual_stat_name

            if pyBigWig_name not in results:
                results[pyBigWig_name] = {}

            # get actual value of the stat
            if actual_stat_name not in actual:
                results[pyBigWig_name]['exact_run_time'], actual[actual_stat_name] = self.benchmark_pyBigWig(actual_stat_name)

            # get corresponding pyBigWig non-exact stat
            if pyBigWig_name not in predictions:
                results[pyBigWig_name]['approx_run_time'], predictions[pyBigWig_name]\
                    = self.benchmark_pyBigWig(actual_stat_name, False)

            # get stat from pyBedGraph
            results[stat_name]['run_time'], predictions[stat_name] =\
                self.benchmark_self(stat_name)

        # find error
        for stat_name in predictions:
            actual_stat_name = stat_name
            if 'mean' in stat_name and stat_name != 'mean':
                actual_stat_name = 'mean'

            if 'pyBigWig_' in stat_name:
                actual_stat_name = stat_name[9:]  # get rid of the pyBigWig tag

            results[stat_name]['error'] = self.get_error(predictions[stat_name],
                                                         actual[actual_stat_name])

        return results

    # do not call
    def output_results(self):

        for stat in range(len(self.stats)):

            baseline_values = self.actual_values[stat]
            baseline_time = self.pyBigWig_times[stat]
            baseline_name = 'pyBigWig'
            my_values = self.self_values[stat]
            my_time = self.self_times[stat]

            if baseline_values is None:
                baseline_values = self.actual_values[EXACT_MEAN_INDEX]
                baseline_time = self.pyBigWig_times[EXACT_MEAN_INDEX]

            print(f"Results for {self.stats[stat]}:")

            mse_values = []
            absolute_error_values = []
            percent_error_values = []
            not_0_values = []
            for i in range(self.num_tests):
                actual = baseline_values[i]
                predicted = my_values[i]

                if actual is None:
                    if predicted is None:
                        continue
                    if stat != APPROX_MEAN_INDEX and stat != MOD_APPROX_MEAN_INDEX:
                        if self.stats[stat] != "coverage":
                            print(predicted, f"is supposed to be None")
                    actual = 0

                if predicted is None:
                    print(f"Correct value is ({actual}) for the range"
                          f" ({self.test_cases[i][0]} - {self.test_cases[i][1]}),"
                          f" but {self.stats[stat]} found ({predicted})")
                    if len(self.intervals_list) > 0:
                        print(self.intervals_list[i])
                    predicted = 0

                mean_squared_error = (actual - predicted) * (actual - predicted)
                mse_values.append(mean_squared_error)

                absolute_error = abs(actual - predicted)
                absolute_error_values.append(absolute_error)

                if actual == 0 and predicted != 0:
                    not_0_values.append(abs(predicted))
                elif predicted == 0:
                    percent_error_values.append(0)
                else:
                    percent_error = abs(actual - predicted) / actual
                    percent_error_values.append(percent_error)

            print(f"Getting {self.stats[stat]} values takes"
                  f" {round(my_time / baseline_time * 100, 2)}% of {baseline_name}'s time")
            print(f"Mean Squared Error: {np.mean(mse_values)}")
            print(f"Mean Squared Error Standard Deviation: {np.std(mse_values)}")
            print(f"Absolute Error: {np.mean(absolute_error_values)}")
            print(f"Absolute Error Standard Deviation: {np.std(absolute_error_values)}")
            print(f"Percent Error: {round(np.mean(percent_error_values) * 100, 2)}%")
            print(f"Percent Error Standard Deviation: {round(np.std(percent_error_values) * 100, 2)}")

            if len(not_0_values) != 0:
                print(f"Number of actual values that were 0 and predicted was not: {len(not_0_values)}")
                print(f"Mean: {np.mean(not_0_values)}")
                print(f"Standard Deviation: {np.std(not_0_values)}")
            print()

    def create_test_cases(self, num_tests, interval_size):

        # random.seed(RANDOM_SEED)

        self.test_cases = []

        # test[0]: start of interval
        # test[1]: end of interval
        for i in range(num_tests):
            test_case = {
                0: random.randint(1, self.chromosome.max_index - interval_size)
            }
            test_case[1] = test_case[0] + interval_size

            self.test_cases.append(test_case)

        self.num_tests = num_tests

    def find_intervals(self):

        self.intervals_list.clear()

        print("Finding intervals using pyBigWig's interval function...")
        start_time = time.time()

        for i in range(self.num_tests):
            intervals = self.bw.intervals(self.chrom_name, self.test_cases[i][0], self.test_cases[i][1])
            self.intervals_list.append(intervals)

        time_taken = time.time() - start_time
        print(f"Time taken to get intervals: {time_taken} seconds")

    def benchmark_pyBigWig(self, stat, want_exact=True):

        if want_exact is False:
            print(f"Finding benchmark for pyBigWig's approximate {stat}...")
        else:
            print(f"Finding benchmark for pyBigWig's exact {stat}...")

        values = []
        start_time = time.time()
        for i in range(self.num_tests):
            value = self.bw.stats(self.chromosome.name, self.test_cases[i][0], self.test_cases[i][1],
                                   type=stat, exact=want_exact)
            values.append(value[0])
        time_taken = time.time() - start_time

        print(f"Time for {stat}: {time_taken} seconds for {self.num_tests} trials\n")
        return time_taken, values

    def benchmark_self(self, stat):

        print(f"Finding pyBedGraph benchmark for {stat}...")

        method = self.bedGraph.get_method(self.chromosome.name, stat)
        values = []

        start_time = time.time()
        for i in range(self.num_tests):
            value = method(self.test_cases[i][0], self.test_cases[i][1])
            values.append(value)
        time_taken = time.time() - start_time

        print(f"Time for {stat}: {time_taken} seconds for {self.num_tests} trials\n")
        return time_taken, values

    @staticmethod
    def get_error(predicted_values, actual_values):

        if len(predicted_values) != len(actual_values):
            print(f"Length of predicted values: {len(predicted_values)}, does"
                  f"not equal length of actual values: {len(actual_values)}")
            return

        percent_error_values = []
        for i in range(len(predicted_values)):
            actual = actual_values[i]
            predicted = predicted_values[i]

            if actual is None:
                actual = 0

            if predicted is None:
                predicted = 0

            if actual == 0 and predicted != 0:
                continue
            elif predicted == 0:
                percent_error_values.append(0)
            else:
                percent_error = abs(actual - predicted) / actual
                percent_error_values.append(percent_error)

        return np.mean(percent_error_values)
