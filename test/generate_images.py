import matplotlib.pyplot as plt
import numpy as np

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
    'pyBigWig_approx',
]

INTERVAL_ERROR_NAMES = [
    'pyBigWig_approx',
]

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


def sample_method():
    # runtime
    rt_x = np.arange(1, 6)
    rt = [[1, 5, 25, 125, 625], [1, 2, 4, 8, 16], [0.5, 1, 2, 4, 8]]
    rt_col = ['blue', 'red', 'green']
    rt_lab = ['pyBW', 'pyBG exact', 'pyBG approx.']
    # error rate
    er_x = [64, 128, 256, 512, 1024]
    er = [[1, 1, 1.2, 1.3, 1.4], [0.1, 0.4, 1.6, 3, 5], [0.1, 0.3, 1.2, 2, 4]]
    er_col = ['blue', 'red', 'green']
    er_lab = ['pyBW approx.', 'pyBG approx.', 'pyBG mod approx.']

    fig = plt.figure(figsize=(8, 12))
    fig.subplots_adjust(hspace=.3)
    plt.subplot(2, 1, 1)

    # plot runtime
    for i in range(len(rt_lab)):
        plt.plot(rt_x, [np.log10(y) for y in rt[i]], color=rt_col[i], marker='o', ms=9, linewidth=2.5, label=rt_lab[i])
    plt.title("Test runtime", fontsize=16)
    plt.xlabel("log10(Sample size)", fontsize=16)
    plt.ylabel("log10(Runtime in seconds)", fontsize=16)
    plt.xticks(np.arange(1, 6, 1), fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc='best', fontsize=16)

    # plot error rate
    plt.subplot(2, 1, 2)
    for i in range(len(er_lab)):
        plt.plot(er_x, er[i], color=er_col[i], marker='o', ms=9, linewidth=2.5, label=er_lab[i])
    plt.title("Test error rate for binsize = " + str(128), fontsize=16)
    plt.xlabel("Interval size (bps)", fontsize=16)
    plt.ylabel("Error rate (%)", fontsize=16)
    plt.xticks(er_x, fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc='best', fontsize=16)

    # plt.show()
    # plt.savefig('pyBedGraph_example_figures.pdf', dpi=300, bbox_inches="tight")
    plt.savefig('graphs/pyBedGraph_example_figures_binsize_' + str(128) + '.pdf', dpi=300)
    plt.close()
