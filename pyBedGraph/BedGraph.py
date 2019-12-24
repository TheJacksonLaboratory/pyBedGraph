# .Dependency for uploading to pypi
# from .Chromosome import Chromosome

# normal import for local use
from .Chrom_Data import Chrom_Data
from .Chrom_Data_Complete import Chrom_Data_Complete
import numpy as np
import time
import os
import logging

CHROM_NAME_INDEX = 0
BUFFER_COUNTER = 10000  # 10^4
BUFFER_COUNTER = 100000  # 10^4

log = logging.getLogger()


class BedGraph:

    """
    min_value only works when loading bedgraphs
    """
    def __init__(self, chrom_size_file_name, data_file_name, chroms_to_load=None,
                 ignore_missing_bp=True, min_value=-1, debug=False):

        file_parts = os.path.basename(data_file_name).split('.')
        using_bigwig = (file_parts[1].lower() == 'bigwig')

        self.name = file_parts[0]
        self.chromosome_map = {}
        self.chrom_sizes = {}
        self.ignore_missing_bp = ignore_missing_bp

        if chroms_to_load:
            chroms_to_load = set(chroms_to_load)

        log.info(f"Reading in {chrom_size_file_name} ...")
        with open(chrom_size_file_name) as chrom_size_file:
            for line in chrom_size_file:
                data = line.split()
                if len(data) != 2:
                    error_msg = f"\n{chrom_size_file} has an incorrect format."\
                                f"It must be in the format:\n"\
                                "chr1 100000\n"\
                                "chr2 50000"
                    log.critical(error_msg)

                    raise RuntimeError(error_msg)

                chrom_name = data[0]
                self.chrom_sizes[chrom_name] = int(data[1])

        if not using_bigwig:
            current_chrom = None
            log.info(f"Reading in {data_file_name} ...")
            with open(data_file_name) as data_file:
                for line in data_file:
                    data = line.split()
                    chrom_name = data[CHROM_NAME_INDEX]

                    if chrom_name not in self.chrom_sizes:
                        log.warning(
                            f"{chrom_name} was not included in {chrom_size_file}")
                        continue

                    if chroms_to_load is not None and chrom_name not in chroms_to_load:
                        continue

                    # Create chromosome object here to trim right after adding
                    if chrom_name not in self.chromosome_map:
                        if ignore_missing_bp:
                            self.chromosome_map[chrom_name] = \
                                Chrom_Data(chrom_name,
                                           self.chrom_sizes[chrom_name],
                                           min_value, debug=debug)
                        else:
                            self.chromosome_map[chrom_name] = \
                                Chrom_Data_Complete(chrom_name,
                                                    self.chrom_sizes[
                                                        chrom_name],
                                                    min_value, debug=debug)

                    if current_chrom is None:
                        current_chrom = self.chromosome_map[chrom_name]

                    if current_chrom.name != chrom_name:
                        current_chrom.trim_extra_space()
                        current_chrom = self.chromosome_map[chrom_name]

                    current_chrom.add_data(data)

                # clean up the last chromosome found in the bedGraph file
                if current_chrom:
                    current_chrom.trim_extra_space()

        else:
            # start_time = time.time()
            import pyBigWig
            bw = pyBigWig.open(data_file_name)

            for chrom_name in self.chrom_sizes:

                if chroms_to_load is not None and \
                        chrom_name not in chroms_to_load:
                    continue

                try:
                    chrom_intervals = bw.intervals(chrom_name)
                except RuntimeError:
                    continue

                if ignore_missing_bp:
                    self.chromosome_map[chrom_name] = \
                        Chrom_Data(chrom_name,
                                   self.chrom_sizes[chrom_name],
                                   min_value, debug=debug)
                else:
                    self.chromosome_map[chrom_name] = \
                        Chrom_Data_Complete(chrom_name,
                                            self.chrom_sizes[chrom_name],
                                            min_value, debug=debug)

                current_chrom = self.chromosome_map[chrom_name]
                if min_value > -1:
                    for interval in chrom_intervals:
                        current_chrom.add_data(interval)
                else:
                    current_chrom.add_bigwig_data(chrom_intervals)
                current_chrom.trim_extra_space()

            # print(time.time() - start_time)

        if chroms_to_load:
            for chrom_name in chroms_to_load:
                if chrom_name not in self.chromosome_map:
                    error_msg = f"{chrom_name} was not found in " \
                                f"{data_file_name}"
                    log.critical(error_msg)

                    raise RuntimeError(error_msg)

    def get_chrom(self, chrom_name):
        return self.chromosome_map[chrom_name]

    def has_chrom(self, chrom_name):
        return chrom_name in self.chromosome_map

    def load_chrom_data(self, chrom_name):
        self.chromosome_map[chrom_name].load_index_array()

    def load_chrom_bins(self, chrom_name, max_bins_size):
        self.chromosome_map[chrom_name].load_bins(max_bins_size)

    def free_chrom_data(self, chrom_name):
        self.chromosome_map[chrom_name].free_index_list()

    def get_method(self, chrom_name, stat):

        if chrom_name not in self.chromosome_map:
            log.error(f"{chrom_name} is not a valid chromosome")
            return None

        chrom = self.chromosome_map[chrom_name]

        if not chrom.loaded_chrom:
            log.error(
                f"{chrom.name} needs to be loaded before it can be searched.")
            return None

        return chrom.get_method(stat)

    # change the shape of intervals to be two lists: start_list, end_list
    @staticmethod
    def change_shape(intervals):
        num_tests = len(intervals)
        start_list = np.zeros(num_tests, dtype=np.int32)
        end_list = np.zeros(num_tests, dtype=np.int32)
        for i in range(num_tests):
            if len(intervals[i]) != 3:
                log.error(
                    f"List given has incorrect formatting. It must be in the format:\n"
                    "[[chr1, 1, 100], [chr1, 101, 200], ...]")
                log.error(intervals[i])
                continue

            start_list[i] = intervals[i][1]
            end_list[i] = intervals[i][2]

        return start_list, end_list

    # can only search one chromosome at a time
    def stats(self, stat="mean", intervals=None, start_list=None,
              end_list=None, chrom_name=None):

        # convert intervals to start_list, end_list
        if intervals is not None:
            chrom_name = intervals[0][0]
            start_list, end_list = self.change_shape(intervals)

        if start_list is None or end_list is None or chrom_name is None:
            log.error(
                "Must either have intervals or start_list, end_list, chrom_name")
            return None

        if not type(start_list) is np.ndarray:
            assert type(start_list) is list
            start_list = np.asarray(start_list, dtype=np.int32)
        if not type(end_list) is np.ndarray:
            assert type(end_list) is list
            end_list = np.asarray(end_list, dtype=np.int32)

        assert end_list.size == start_list.size
        method_to_call = self.get_method(chrom_name, stat)

        if method_to_call is None:
            return

        start_time = time.time()
        result = method_to_call(start_list, end_list)
        # log.info(f"Time for {stat}:", time.time() - start_time)
        return result

    # output to output_file if given, otherwise return a list of results
    def stats_from_file(self, interval_file, output_to_file=True, stat="mean"):
        results = {}
        test_intervals = {}

        with open(interval_file) as in_file:

            current_chrom_name = None
            current_interval = None

            for line in in_file:
                interval = line.split()

                if len(interval) != 3:
                    log.error(
                        f"{interval_file} has incorrect formatting. It must be in the format:\n"
                        "chr1\t100\t401\n"
                        "chr1\t600\t1000\n"
                        "chr2\t0\t1000\n"
                        "...")
                    log.error(interval)
                    break

                # change current_chrom if interval wants a different chrom
                if current_chrom_name is None or current_chrom_name != interval[
                    0]:
                    current_chrom_name = interval[0]
                    if current_chrom_name not in test_intervals:
                        test_intervals[current_chrom_name] = {}
                        current_interval = test_intervals[current_chrom_name]
                        current_interval['start_list'] = []
                        current_interval['end_list'] = []

                current_interval['start_list'].append(int(interval[1]))
                current_interval['end_list'].append(int(interval[2]))

        for chrom_name in test_intervals:
            method_to_call = self.get_method(chrom_name, stat)

            if method_to_call is None:
                return

            start_time = time.time()
            result = method_to_call(
                np.array(test_intervals[chrom_name]['start_list'],
                         dtype=np.int32),
                np.array(test_intervals[chrom_name]['end_list'],
                         dtype=np.int32))
            log.info(f"Time for {stat}: {time.time() - start_time}")

            results[chrom_name] = result
            if output_to_file:
                with open(chrom_name + '_out.txt', 'w') as out_file:
                    result_counter = 0

                    output = ''
                    for value_index in range(result.size):
                        output += f"{result[value_index]}\n"
                        result_counter += 1

                        if result_counter == BUFFER_COUNTER:
                            result_counter = 0
                            out_file.write(output)
                            output = ""

                    if result_counter > 0:
                        out_file.write(output)

        return results
