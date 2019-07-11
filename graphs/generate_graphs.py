import matplotlib.pyplot as plt
import numpy as np
import os
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.offline

colors = {
    'pyBW exact': 'orange',
    'pyBW app.': 'red',
    'pyBG exact': 'green',
    'pyBG app. bin=int_size/5': 'cyan',
    'pyBG app. bin=int_size/10': 'magenta',
    'pyBG app. bin=int_size/20': 'blue',
    'pyBG app. bin=100': 'cyan',
    'pyBG app. bin=50': 'magenta',
    'pyBG app. bin=25': 'blue'
}

TABLE_HEADERS = [
    'Interval Size (bPS)',
    'Error Rate (%)',
    'Mean Squared Error',
    'Absolute Error',
    '# Actual is 0'
]

TITLE_FONT_SIZE = 16
AXIS_FONT_SIZE = 12
LEGEND_FONT_SIZE = 10

RUN_TIME_NAMES = [
    'pyBW exact',
    'pyBG exact',
    'pyBW app.'
]

INTERVAL_ERROR_NAMES = [
    'pyBW app.'
]

INTERVAL_RUNTIME_NAMES = [
    'pyBW exact',
    'pyBG exact',
    'pyBW app.'
]

GRAPH_ROOT_LOCATION = 'graphs'
NUM_ERROR_TYPES = 4


# Interval size: 500
def create_runtime_num_test(infile, data_name):
    line = infile.readline()
    num_tests = [int(x) for x in line.split()]

    run_times = {}
    while True:
        name = infile.readline().strip().strip()

        if name == "":
            break

        results = infile.readline().split()
        run_times[name] = [float(x) for x in results]

    i = 0
    for name in run_times:
        plt.plot([np.log10(x) for x in num_tests], [np.log10(x) for x in run_times[name]],
                 color=colors[name], label=name)
        i += 1

    plt.title(f"Run Time for {data_name}", fontsize=TITLE_FONT_SIZE)
    plt.xlabel("log10(# of tests)", fontsize=AXIS_FONT_SIZE)
    plt.ylabel("log10(runtime (seconds))", fontsize=AXIS_FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/run_time.png')
    plt.close()


# Number of tests: 10,000
def create_interval_error(in_file, data_name):
    line = in_file.readline()
    intervals = [int(x) for x in line.split()]

    errors = {}
    table_cells = {}
    while True:
        name = in_file.readline().strip()

        if name == "":
            break

        table_cells[name] = [[x] for x in intervals]

        error_list = []
        for i in range(NUM_ERROR_TYPES):
            error_list.append(in_file.readline().split())

        errors[name] = {}
        for line in error_list:
            error_name = line[0]
            errors[name][error_name] = [float(line[x]) for x in range(1, len(line), 1)]
            for x in range(len(intervals)):
                error = errors[name][error_name][x]

                if error_name == 'num_actual_0':
                    error = int(error)

                if error_name == 'percent_error':
                    error *= 100
                error = round(error, 5)

                table_cells[name][x].append(error)

    i = 0
    for name in errors:
        plt.plot(intervals, [x * 100 for x in errors[name]['percent_error']], color=colors[name], label=name)
        i += 1

    plt.title(f"Error vs. Interval Size for {data_name}", fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Interval Size (basePairs)", fontsize=AXIS_FONT_SIZE)
    plt.ylabel("Percentage Error Rate (%)", fontsize=AXIS_FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/interval_error.png')
    plt.close()

    for name in table_cells:
        plt.figure()
        plt.title(f"{data_name} --- {name}", fontsize=AXIS_FONT_SIZE)
        table = plt.table(
            cellText=table_cells[name],
            colWidths=[0.027, 0.023, 0.03, 0.022, 0.02],
            colLabels=TABLE_HEADERS,
            loc='center'
        )
        table.auto_set_font_size(False)
        table.set_fontsize(LEGEND_FONT_SIZE)
        table.scale(11, 2)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)

        for pos in ['right', 'top', 'bottom', 'left']:
            plt.gca().spines[pos].set_visible(False)

        save_name = name.replace(" ", "_").replace('/', '')
        plt.savefig(f'graphs/{data_name}/{save_name}_table.png',
                    bbox_inches='tight', pad_inches=0.05)
        plt.close()


# Number of tests: 10,000
def create_interval_runtime(in_file, data_name):
    line = in_file.readline()
    intervals = [int(x) for x in line.split()]

    run_times = {}
    while True:
        name = in_file.readline().strip()

        if name == '':
            break

        results = in_file.readline().split()
        run_times[name] = [float(x) for x in results]

    i = 0
    for name in run_times:
        plt.plot(intervals, [np.log10(x) for x in run_times[name]],
                 color=colors[name], label=name)
        i += 1

    plt.title(f"Interval vs. Run Time for {data_name}", fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Interval Size", fontsize=AXIS_FONT_SIZE)
    plt.ylabel("log10(runtime (seconds))", fontsize=AXIS_FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/interval_run_time.png')
    plt.close()


def create_values_indexed(in_file, data_name):

    bin_sizes = in_file.readline().split()
    bin_sizes = [int(x) for x in bin_sizes]

    values_indexed = in_file.readline().split()
    values_indexed = [int(x) for x in values_indexed]

    plt.plot(bin_sizes, values_indexed)

    plt.title(f"Values Indexed vs. Bin Size for Exact Mean Calculation", fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Bin Size", fontsize=AXIS_FONT_SIZE)
    plt.ylabel("Values Indexed", fontsize=AXIS_FONT_SIZE)
    plt.savefig(f'graphs/{data_name}/values_indexed.png')
    plt.close()


def main():
    create_values_indexed(open('graphs/ENCFF376VCU/values_indexed.txt'), 'ENCFF376VCU')
    exit()

    for subdir, dirs, files in os.walk(GRAPH_ROOT_LOCATION):
        data_name = subdir[7:]
        print(data_name)
        for file_name in files:
            file_path = subdir + '/' + file_name

            with open(file_path) as in_file:
                if file_name == 'run_time_results.txt':
                    create_runtime_num_test(in_file, data_name)
                elif file_name == 'interval_error_results.txt':
                    create_interval_error(in_file, data_name)
                elif file_name == 'interval_runtime_results.txt':
                    create_interval_runtime(in_file, data_name)
                elif file_name[-4:] == '.png':
                    continue
                else:
                    print(f"Unknown file: {file_name}")


if __name__ == '__main__':
    main()
