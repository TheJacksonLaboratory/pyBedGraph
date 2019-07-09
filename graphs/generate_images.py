import matplotlib.pyplot as plt
import numpy as np
import os

colors = [
    'blue',
    'red',
    'green',
    'cyan',
    'magenta',
    'yellow',
    'orange',
    'black',
    'darkgreen'
]

FONT_SIZE = 16
LEGEND_FONT_SIZE = 8

RUN_TIME_NAMES = [
    'pyBigWig_exact',
    'pyBedGraph_exact',
    'pyBigWig_approx'
]

INTERVAL_ERROR_NAMES = [
    'pyBigWig_approx'
]

INTERVAL_RUNTIME_NAMES = [
    'pyBigWig_exact',
    'pyBedGraph_exact',
    'pyBigWig_approx'
]

GRAPH_ROOT_LOCATION = 'graphs'
NUM_ERROR_TYPES = 4


# Interval size: 500
def create_runtime_num_test(infile, data_name):
    line = infile.readline()
    num_tests = [int(x) for x in line.split()]

    run_times = {}
    for line in infile:
        results = line.split()
        name = results[0]
        run_times[name] = [float(results[x]) for x in range(1, len(results), 1)]

    i = 0
    for name in run_times:
        plt.plot([np.log10(x) for x in num_tests], [np.log10(x) for x in run_times[name]],
                 color=colors[i], label=name)
        i += 1

    plt.title(f"Run Time Test for {data_name}")
    plt.xlabel("# of Tests (log10)", fontsize=FONT_SIZE)
    plt.ylabel("Runtime (log10 sec)", fontsize=FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/run_time.png')
    plt.close()


# Number of tests: 10,000
def create_interval_error(in_file, data_name):
    line = in_file.readline()
    intervals = [int(x) for x in line.split()]

    errors = {}
    while(True):
        name = in_file.readline()

        if name == "":
            break

        name = name.strip()
        error_list = []
        for i in range(NUM_ERROR_TYPES):
            error_list.append(in_file.readline().split())

        errors[name] = {}
        for line in error_list:
            error_name = line[0]
            errors[name][error_name] = [float(line[x]) for x in range(1, len(line), 1)]

    i = 0
    for name in errors:
        # only include error of less than 100%
        '''interval_list = []
        error_list = []
        for test_index in range(len(data[name])):
            error = data[name][test_index]
            if error > 1:
                if len(interval_list) == 0:
                    continue
                break
            interval_list.append(interval_test_list[test_index])
            error_list.append(error)'''

        # plt.plot(interval_test_list, [np.log10(x) for x in data[name]], color=colors[i], label=names[i])
        plt.plot(intervals, errors[name]['percent_error'], color=colors[i], label=name)
        i += 1

    plt.title(f"Error vs. Interval Size for {data_name}")
    plt.xlabel("Interval Size", fontsize=FONT_SIZE)
    plt.ylabel("Percentage Error Rate", fontsize=FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/interval_error.png')
    plt.close()


# Number of tests: 10,000
def create_interval_runtime(in_file, data_name):
    line = in_file.readline()
    intervals = [int(x) for x in line.split()]

    run_times = {}
    for line in in_file:
        results = line.split()
        name = results[0]
        run_times[name] = [float(results[x]) for x in range(1, len(results), 1)]

    i = 0
    for name in run_times:
        plt.plot(intervals, [np.log10(x) for x in run_times[name]],
                 color=colors[i], label=name)
        i += 1

    plt.title(f"Interval Run Time Test for {data_name}")
    plt.xlabel("Interval Size", fontsize=FONT_SIZE)
    plt.ylabel("Runtime (log10 sec)", fontsize=FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/interval_run_time.png')
    plt.close()


def main():
    for subdir, dirs, files in os.walk(GRAPH_ROOT_LOCATION):
        data_name = subdir[7:]
        for file_name in files:
            file_path = subdir + '/' + file_name

            with open(file_path) as in_file:
                if file_name == 'run_time_results.txt':
                    create_runtime_num_test(in_file, data_name)
                elif file_name == 'interval_error_results.txt':
                    create_interval_error(in_file, data_name)
                elif file_name == 'interval_runtime_results.txt':
                    create_interval_runtime(in_file, data_name)
                else:
                    print(f"Unknown file: {file_name}")


if __name__ == '__main__':
    main()
