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

def create_plots():
    with open('graphs/ENCFF487MJL/run_time_data.txt') as input_file:
        line = input_file.readline()
        num_tests = line.split()
        num_tests = [np.log10(int(x)) for x in num_tests]

        run_times = []
        for line in input_file:
            values = line.split()
            run_times.append([float(x) for x in values])

        create_runtime_num_test('ENCFF487MJL', num_tests, run_times)


# Interval size: 500
def create_runtime_num_test(data_name, num_tests, run_times):
    names = list(run_times.keys())
    for i in range(len(names)):
        name = names[i]
        plt.plot([np.log10(x) for x in num_tests], [np.log10(x) for x in run_times[name]],
                 color=colors[i], label=names[i])
    plt.title(f"Runtime Test for {data_name}")
    plt.xlabel("# of tests (log10)", fontsize=FONT_SIZE)
    plt.ylabel("Runtime (log10 sec)", fontsize=FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/run_time.png')
    plt.close()


# Number of tests: 100,000
def create_error_interval(data_name, interval_test_list, data):

    names = list(data.keys())
    for i in range(len(names)):
        name = names[i]

        # only include error of less than 100%
        interval_list = []
        error_list = []
        for test_index in range(len(data[name])):
            error = data[name][test_index]
            if error > 1:
                if len(interval_list) == 0:
                    continue
                break
            interval_list.append(interval_test_list[test_index])
            error_list.append(error)

        # plt.plot(interval_test_list, [np.log10(x) for x in data[name]], color=colors[i], label=names[i])
        plt.plot(interval_list, error_list, color=colors[i], label=names[i])
    plt.title(f"Error vs. Interval Size for {data_name}")
    plt.xlabel("Interval Size", fontsize=FONT_SIZE)
    plt.ylabel("Percentage Error Rate", fontsize=FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/error_interval.png')
    plt.close()


#create_plots()

def main():
    pass


if __name__ == '__main__':
    main()
