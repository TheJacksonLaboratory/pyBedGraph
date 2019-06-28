import random
import pyBigWig
import time
import numpy as np
import sys

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
        self.find_actual_values()
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

            actual_error_values = []
            percent_error_values = []
            for i in range(self.num_tests):
                actual = baseline_values[i]
                predicted = my_values[i]

                if actual is None:
                    if predicted is None:
                        continue
                    if stat != APPROX_MEAN_INDEX and stat != MOD_APPROX_MEAN_INDEX:
                        print(predicted, f"is supposed to be None")
                    actual = 0

                if predicted is None:
                    print(f"Correct value is ({actual}) for the range"
                          f" ({self.test_cases[i][0]} - {self.test_cases[i][1]}),"
                          f" but {self.stats[stat]} found ({predicted})")
                    predicted = 0

                actual_error = (actual - predicted) * (actual - predicted)
                actual_error_values.append(actual_error)

                # equation from stackoverflow, calculating relative error when true value is 0
                if actual - predicted == 0:
                    percent_error = 0
                else:
                    percent_error = abs((actual - predicted) / (actual + predicted) * 2)
                percent_error_values.append(percent_error)

            print(f"Getting {self.stats[stat]} values takes"
                  f" {round(my_time / baseline_time * 100, 2)}% of {baseline_name}'s time")
            print(f"Squared Error: {np.mean(actual_error_values)}")
            print(f"Squared Error Standard Deviation: {np.std(actual_error_values)}")
            print(f"Percent Error: {round(np.mean(percent_error_values) * 100, 2)}%")
            print(f"Percent Error Standard Deviation: {round(np.std(percent_error_values) * 100, 2)}\n")

    def get_mean_from_intervals(self):

        mean_values = []
        for i in range(self.num_tests):
            if self.intervals_list[i] is None:
                mean_values.append(None)
                continue

            test_case = self.test_cases[i]
            numb = 0
            mean = 0
            for interval in self.intervals_list[i]:
                start = interval[0]
                end = interval[1]
                if interval[0] < test_case[0]:
                    start = test_case[0]

                if interval[1] > test_case[1]:
                    end = test_case[1]

                numb += end - start
                mean += (end - start) * interval[2]

            mean /= numb
            mean_values.append(mean)

        self.actual_values.append(mean_values)

    def get_max_from_intervals(self):
        max_values = []
        for i in range(self.num_tests):
            if self.intervals_list[i] is None:
                max_values.append(None)
                continue

            max_value = -1
            for interval in self.intervals_list[i]:
                if interval[2] > max_value:
                    max_value = interval[2]

            if max_value == -1:
                max_value = None

            max_values.append(max_value)

        self.actual_values.append(max_values)

    def get_min_from_intervals(self):
        min_values = []
        for i in range(self.num_tests):
            if self.intervals_list[i] is None:
                min_values.append(None)
                continue

            min_value = sys.maxsize
            for interval in self.intervals_list[i]:
                if interval[2] < min_value:
                    min_value = interval[2]

            if min_value == sys.maxsize:
                min_value = None

            min_values.append(min_value)

        self.actual_values.append(min_values)

    def get_coverage_from_intervals(self):
        coverage_values = []
        for i in range(self.num_tests):
            intervals = self.intervals_list[i]
            if intervals is None:
                coverage_values.append(0)
                continue

            test_case = self.test_cases[i]
            total_size = test_case[1] - test_case[0]
            covered = 0
            for interval in intervals:
                start = interval[0]
                end = interval[1]
                if interval[0] < test_case[0]:
                    start = test_case[0]

                if interval[1] > test_case[1]:
                    end = test_case[1]

                covered += end - start

            coverage_values.append(covered / total_size)

        self.actual_values.append(coverage_values)

    def get_std_from_intervals(self):
        std_values = []
        for i in range(self.num_tests):
            if self.intervals_list[i] is None:
                std_values.append(None)
                continue

            values = []
            test_case = self.test_cases[i]
            for interval in self.intervals_list[i]:
                start = interval[0]
                end = interval[1]
                if interval[0] < test_case[0]:
                    start = test_case[0]

                if interval[1] > test_case[1]:
                    end = test_case[1]

                interval_length = end - start
                to_append = [interval[2]] * interval_length
                values = values + to_append

            std = np.std(values)
            std_values.append(std)

        self.actual_values.append(std_values)

    def create_test_cases(self, interval_size):
        # test[0]: start of interval
        # test[1]: end of interval
        for i in range(self.num_tests):
            test_case = {
                0: random.randint(1, self.genome.chromosome.size - interval_size)
            }
            test_case[1] = test_case[0] + interval_size

            self.test_cases.append(test_case)

    def find_actual_values(self):

        print("Finding actual values using bigWig's intervals function...")
        start_time = time.time()

        for i in range(self.num_tests):
            intervals = self.bw.intervals(self.genome.chromosome.name, self.test_cases[i][0], self.test_cases[i][1])
            self.intervals_list.append(intervals)

        time_taken = time.time() - start_time
        print(f"Time taken to get intervals: {time_taken} seconds")

        self.get_mean_from_intervals()
        self.actual_values.append(None)
        self.actual_values.append(None)
        self.get_max_from_intervals()
        self.get_min_from_intervals()
        self.get_coverage_from_intervals()
        self.get_std_from_intervals()

        print("Done.\n")

    def benchmark_pyBigWig(self):

        for stat_name in self.stats:
            if stat_name == "approx_mean" or stat_name == "mod_approx_mean":
                self.pyBigWig_times.append(None)
                continue

            print(f"Finding benchmark for pyBigWig's {stat_name}...")

            start_time = time.time()
            for i in range(self.num_tests):
                self.bw.stats(self.genome.chromosome.name, self.test_cases[i][0], self.test_cases[i][1],
                         type=stat_name)
            time_taken = time.time() - start_time

            self.pyBigWig_times.append(time_taken)

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
