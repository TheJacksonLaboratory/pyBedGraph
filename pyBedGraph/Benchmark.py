import random
import pyBigWig
import time


class Benchmark:

    def __init__(self, chrom, num_tests, bin_size, interval_size, bigWig_file, baseline="pyBigWig"):
        self.pyBigWig_time = None
        self.exact_time = None
        self.approx_time = None
        self.mod_approx_time = None
        self.pyBigWig_values = []
        self.exact_values = []
        self.approx_values = []
        self.mod_approx_values = []

        self.test_cases = []

        self.num_tests = num_tests
        self.chrom = chrom

        print(f"Benchmark for interval size: {interval_size} | bin size: {bin_size}")
        chrom.split_bins(bin_size)

        self.create_test_cases(interval_size)
        self.find_pyBigWig(bigWig_file)
        self.find_approx()
        self.find_mod_approx()
        self.find_exact()

        self.output_results(baseline, bigWig_file)

    def output_results(self, baseline, bigWig_file):

        # baseline can only be exact or pyBigWig
        baseline_values = self.pyBigWig_values
        if baseline == "exact":
            baseline_values = self.exact_values

        found = {}
        found["pyBigWig"] = (self.pyBigWig_values, self.pyBigWig_time)
        found["exact"] = (self.exact_values, self.exact_time)
        found["approx"] = (self.approx_values, self.approx_time)
        found["mod_approx"] = (self.mod_approx_values, self.mod_approx_time)

        bw = pyBigWig.open(bigWig_file)

        print("Results:")
        for data_set in found:

            if data_set == baseline:
                continue

            values = found[data_set][0]

            squared_error = 0
            percent_error = 0
            num = 0
            for i in range(self.num_tests):
                if baseline_values[i] is None or not baseline_values[i] >= 0:
                    continue

                if values[i] is None or not values[i] >= 0:
                    print(f"Baseline found a value ({baseline_values[i]}) for the range"
                          f" ({self.test_cases[i][0]} - {self.test_cases[i][1]}),"
                          f" but {data_set} found ({values[i]})")
                    continue

                squared_error += (baseline_values[i] - values[i]) * (baseline_values[i] - values[i])
                percent_error += abs(baseline_values[i] - values[i]) / baseline_values[i]
                num += 1

                # if percent_error > 0.001 and data_set == 'exact':
                #   intervals = bw.intervals(self.chrom.name, self.test_cases[i][0], self.test_cases[i][1])
                #   self.test(intervals, self.test_cases[i], baseline_values[i], values[i])

            print(f"Getting {data_set} values takes"
                  f" {round(found[data_set][1] / found[baseline][1] * 100, 2)}% of {baseline}'s time")
            print(f"Squared Error: {squared_error / num}")
            print(f"Percent Error: {round(percent_error / num * 100, 2)}%\n")

    @staticmethod
    def test(intervals, testcase, baseline, predicted):
        average = 0
        numb = 0
        for interval in intervals:
            start = interval[0]
            end = interval[1]
            if interval[0] < testcase[0]:
                start = testcase[0]

            if interval[1] > testcase[1]:
                end = testcase[1]

            numb += end - start
            average += (end - start) * interval[2]
        print(average / numb, baseline, predicted)
        print()


    def create_test_cases(self, interval_size):
        # test[0]: start of interval
        # test[1]: end of interval
        for i in range(self.num_tests):
            test_case = {}
            test_case[0] = random.randint(1, self.chrom.size - interval_size)
            test_case[1] = test_case[0] + interval_size

            self.test_cases.append(test_case)

    def find_pyBigWig(self, input_file):
        print("Finding benchmark for pyBigWig...")
        bw = pyBigWig.open(input_file)

        start_time = time.time()
        for i in range(self.num_tests):
            value = bw.stats(self.chrom.name, self.test_cases[i][0], self.test_cases[i][1])
            self.pyBigWig_values.append(value[0])
        self.pyBigWig_time = time.time() - start_time

        print(f"pyBigWig: {self.pyBigWig_time} seconds for {self.num_tests} trials\n")

    def find_approx(self):
        print("Finding benchmark for approx...")
        start_time = time.time()
        for i in range(self.num_tests):
            value = self.chrom.get_approx_average(self.test_cases[i][0], self.test_cases[i][1])
            self.approx_values.append(value)
        self.approx_time = time.time() - start_time

        print(f"Time for approx_average: {self.approx_time} seconds for {self.num_tests} trials\n")

    def find_mod_approx(self):
        print("Finding benchmark for modified approx...")
        start_time = time.time()
        for i in range(self.num_tests):
            value = self.chrom.get_mod_approx_average(self.test_cases[i][0], self.test_cases[i][1])
            self.mod_approx_values.append(value)
        self.mod_approx_time = time.time() - start_time

        print(f"Time for mod_approx_average: {self.mod_approx_time} seconds for {self.num_tests} trials\n")

    def find_exact(self):
        print("Finding benchmark for exact...")
        start_time = time.time()
        for i in range(self.num_tests):
            value = self.chrom.get_exact_average(self.test_cases[i][0], self.test_cases[i][1])
            self.exact_values.append(value)
        self.exact_time = time.time() - start_time

        print(f"Time for approx_average: {self.exact_time} seconds for {self.num_tests} trials\n")
