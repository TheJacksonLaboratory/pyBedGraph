import pyBigWig
import time
import numpy as np

EXACT_MEAN_INDEX = 0
APPROX_MEAN_INDEX = 1
MOD_APPROX_MEAN_INDEX = 2
RANDOM_SEED = 1

ERROR = 0.00000001

mean_names = [
    'pyBigWig_exact',
    'pyBigWig_approx',
    'pyBedGraph_exact',
    'pyBedGraph_approx'
]

ALL_STATS = [
    "mean",
    "approx_mean",
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
        self.test_cases = None

        self.chromosome = None

    def benchmark(self, num_tests, interval_size, chrom_name, bin_size=None,
                  stats=None, just_runtime=False, bench_pyBigWig_approx=True,
                  make_pyBigWig_baseline=True):

        # benchmark all stats if none given
        if stats is None:
            stats = ALL_STATS

        print("Benchmarking:\n"
              f"Number of tests: {num_tests}\n"
              f"Interval size: {interval_size}\n"
              f"Chromosome name: {chrom_name}\n"
              f"Bin size: {bin_size}\n"
              f"Stats to bench: {stats}\n"
              f"Just bench run time: {just_runtime}\n"
              f"Bench pyBigWig approx: {bench_pyBigWig_approx}\n"
              f"Baseline is pyBigWig exact: {make_pyBigWig_baseline}\n")

        # self.find_intervals()

        self.chromosome = self.bedGraph.chromosome_map[chrom_name]
        self.bedGraph.load_chrom_data(chrom_name)
        if 'approx_mean' in stats:
            if bin_size is None:
                print("Must give a bin_size if benchmarking approx_mean")
                return
            self.bedGraph.load_chrom_bins(chrom_name, bin_size)

        self.create_test_cases(num_tests, interval_size)

        results = {}
        actual = {}
        predictions = {}

        for stat_name in stats:

            results[stat_name] = {}

            actual_stat_name = stat_name
            if 'mean' in stat_name and stat_name != 'mean':
                actual_stat_name = 'mean'

            # get actual value of the stat from pyBedGraph if wanted
            if actual_stat_name not in actual:
                if not make_pyBigWig_baseline:
                    actual[actual_stat_name] = \
                        self.bedGraph.stats(actual_stat_name,
                                            start_list=self.test_cases[0],
                                            end_list=self.test_cases[1],
                                            chrom_name=chrom_name)

            pyBigWig_name = 'pyBigWig_' + actual_stat_name
            if bench_pyBigWig_approx:
                if pyBigWig_name not in results:
                    results[pyBigWig_name] = {}

                # get corresponding pyBigWig non-exact stat
                if pyBigWig_name not in predictions:
                    results[pyBigWig_name]['approx_run_time'], predictions[pyBigWig_name] =\
                        self.benchmark_pyBigWig(actual_stat_name, False)

            # get actual stat from pyBigWig if haven't gotten it yet
            if actual_stat_name not in actual:
                if pyBigWig_name not in results:
                    results[pyBigWig_name] = {}

                results[pyBigWig_name]['exact_run_time'], actual[actual_stat_name] = \
                    self.benchmark_pyBigWig(actual_stat_name)

            # get stat from pyBedGraph
            start_time = time.time()
            predictions[stat_name] = self.bedGraph.stats(stat_name,
                                                         start_list=self.test_cases[0],
                                                         end_list=self.test_cases[1],
                                                         chrom_name=chrom_name)
            results[stat_name]['run_time'] = time.time() - start_time

        if just_runtime:
            return results

        # find error
        for stat_name in predictions:
            actual_stat_name = stat_name
            if 'mean' in stat_name and stat_name != 'mean':
                actual_stat_name = 'mean'

            care_about_high_error = True
            if 'pyBigWig_' in stat_name:
                actual_stat_name = stat_name[9:]  # get rid of the pyBigWig tag
                care_about_high_error = False

            percent_error, ms_error, abs_error, not_included =\
                self.get_error(predictions[stat_name], actual[actual_stat_name],
                               care_about_high_error)

            results[stat_name]['error'] = {}
            r = results[stat_name]['error']
            r['percent_error'] = percent_error
            r['ms_error'] = ms_error
            r['abs_error'] = abs_error
            r['not_included'] = not_included

        return results

    def create_test_cases(self, num_tests, interval_size):

        np.random.seed(RANDOM_SEED)

        self.num_tests = num_tests
        # change maximum to chrom size file
        test_cases = np.random.randint(0, self.chromosome.max_index - interval_size,
                                       num_tests, dtype=np.int32)
        self.test_cases = np.vstack((test_cases, test_cases + interval_size))

    def find_intervals(self):

        self.intervals_list.clear()

        print("Finding intervals using pyBigWig's interval function...")
        start_time = time.time()

        for i in range(self.num_tests):
            intervals = self.bw.intervals(self.chromosome.name, self.test_cases[0][i], self.test_cases[1][i])
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
            value = self.bw.stats(self.chromosome.name, self.test_cases[0][i], self.test_cases[1][i],
                                   type=stat, exact=want_exact)
            values.append(value[0])

            '''start = self.test_cases[0][i]
            end = self.test_cases[1][i]

            intervals = self.bw.intervals(self.chromosome.name, start, end)
            current_interval = intervals[0]
            interval_index = 0
            total = 0
            count = 0
            for i in range(start, end, 1):
                if i == current_interval[1]:
                    interval_index += 1

                    if interval_index == len(intervals):
                        break

                    current_interval = intervals[interval_index]

                if current_interval[0] <= i < current_interval[1]:
                    total += current_interval[2]
                    count += 1'''

        time_taken = time.time() - start_time

        print(f"Time for {stat}: {time_taken} seconds for {self.num_tests} trials\n")
        return time_taken, values

    def benchmark_self(self, stat):

        print(f"Finding pyBedGraph benchmark for {stat}...")

        method = self.bedGraph.get_method(self.chromosome.name, stat)
        values = []

        start_time = time.time()
        for i in range(self.num_tests):
            value = method(self.test_cases[0][i], self.test_cases[1][i])
            values.append(value)
        time_taken = time.time() - start_time

        print(f"Time for {stat}: {time_taken} seconds for {self.num_tests} trials\n")
        return time_taken, values

    def get_error(self, predicted_values, actual_values, care_about_high_error):

        if len(predicted_values) != len(actual_values):
            print(f"Length of predicted values: {len(predicted_values)}, does"
                  f"not equal length of actual values: {len(actual_values)}")
            return

        percent_error_values = []
        mse_error_values = []
        abs_error_values = []
        not_included = 0
        for i in range(len(predicted_values)):
            actual = actual_values[i]
            predicted = predicted_values[i]

            if actual is None or actual == -1:
                actual = 0

            if predicted is None or predicted == -1 or predicted <= ERROR:
                predicted = 0

            abs_error_values.append(abs(actual - predicted))
            mse_error_values.append((actual - predicted) * (actual - predicted))

            if actual == 0 and predicted > ERROR:
                not_included += 1
                #start = self.test_cases[0][i]
                #end = self.test_cases[1][i]
                #bin_start = int(start / 1000)
                #bin_end = int(end / 1000)
                #bin_list = self.bedGraph.chromosome_map['chr1'].bins_list[0]
                #bin_coverage_list = self.bedGraph.chromosome_map['chr1'].bins_list_coverages[0]
                #print(i, bin_list[bin_start:bin_end], bin_coverage_list[bin_start:bin_end], predicted,
                #      self.bedGraph.stats('approx_mean', intervals=[['chr1', start, end]]))
                continue
            elif actual == 0 and predicted <= ERROR:
                percent_error_values.append(0)
            else:
                percent_error = abs(actual - predicted) / actual

                # if care_about_high_error and percent_error > 0.0001:
                #   print(percent_error)

                percent_error_values.append(percent_error)

        return np.mean(percent_error_values), np.mean(mse_error_values),\
            np.mean(abs_error_values), not_included
