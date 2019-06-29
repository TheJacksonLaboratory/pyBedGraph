import random
import pyBigWig
import time
import numpy as np

EXACT_MEAN_INDEX = 0
APPROX_MEAN_INDEX = 1
MOD_APPROX_MEAN_INDEX = 2


class Benchmark:

    def __init__(self, genome, num_tests, interval_size, bigWig_file, test='full'):
        self.stats = [
            "mean",
            "approx_mean",
            "mod_approx_mean",
            "max",
            "min",
            "coverage",
            "std"
        ]

        if test == 'approx':
            self.stats = [
                "mean",
                "approx_mean",
                "mod_approx_mean"
            ]

        self.bw = pyBigWig.open(bigWig_file)
        self.pyBigWig_times = []

        self.actual_values = []
        self.self_values = []
        self.self_times = []

        self.test_cases = []
        self.intervals_list = []

        self.num_tests = num_tests
        self.genome = genome

        print(f"Benchmark for interval size: {interval_size}")

        self.create_test_cases(interval_size)
        self.find_intervals()
        self.benchmark_pyBigWig()
        self.benchmark_self()

        self.output_results()

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

    def create_test_cases(self, interval_size):
        # test[0]: start of interval
        # test[1]: end of interval
        for i in range(self.num_tests):
            test_case = {
                0: random.randint(1, self.genome.chromosome.size - interval_size)
            }
            test_case[1] = test_case[0] + interval_size

            self.test_cases.append(test_case)

    def find_intervals(self):

        print("Finding actual values using bigWig's intervals function...")
        start_time = time.time()

        for i in range(self.num_tests):
            intervals = self.bw.intervals(self.genome.chromosome.name, self.test_cases[i][0], self.test_cases[i][1])
            self.intervals_list.append(intervals)

        time_taken = time.time() - start_time
        print(f"Time taken to get intervals: {time_taken} seconds")

    def benchmark_pyBigWig(self):

        for stat_name in self.stats:
            if stat_name == "approx_mean" or stat_name == "mod_approx_mean":
                self.pyBigWig_times.append(None)
                self.actual_values.append(None)
                continue

            print(f"Finding benchmark for pyBigWig's {stat_name}...")

            start_time = time.time()
            values = []
            for i in range(self.num_tests):
                value = self.bw.stats(self.genome.chromosome.name, self.test_cases[i][0], self.test_cases[i][1],
                              type=stat_name, exact=True)
                values.append(value[0])
            time_taken = time.time() - start_time

            self.pyBigWig_times.append(time_taken)
            self.actual_values.append(values)

            print(f"Time for {stat_name}: {time_taken} seconds for {self.num_tests} trials\n")

    def benchmark_self(self):
        for stat in self.stats:
            print(f"Finding self benchmark for {stat}...")

            method = self.genome.get_method(stat)
            values = []

            start_time = time.time()
            for i in range(self.num_tests):
                value = method(self.test_cases[i][0], self.test_cases[i][1])
                values.append(value)
            time_taken = time.time() - start_time

            self.self_values.append(values)
            self.self_times.append(time_taken)

            print(f"Time for {stat}: {time_taken} seconds for {self.num_tests} trials\n")
