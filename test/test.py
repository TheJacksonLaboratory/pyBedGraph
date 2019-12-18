import sys
import time
import numpy as np

sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph

# arg1 - chromosome sizes file
# arg2 - bedgraph file
# arg3 - (optional) chromosome_name
# Just load chromosome 'chr1' (uses less memory and takes less time)
bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bedGraph', 'chr1')

# Load the whole bedGraph file
bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bedGraph')

# Option to not ignore missing basePairs when calculating statistics
# Used the exact same way but produces slightly different results
inclusive_bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bedGraph', ignore_missing_bp=False)

bedGraph.load_chrom_data('chr1')
inclusive_bedGraph.load_chrom_data('chr1')

bedGraph.load_chrom_bins('chr1', 3)
inclusive_bedGraph.load_chrom_bins('chr1', 3)

# Option 1
test_intervals = [
    ['chr1', 24, 26],
    ['chr1', 12, 15],
    ['chr1', 8, 12],
    ['chr1', 9, 10],
    ['chr1', 0, 5]
]
values = bedGraph.stats(intervals=test_intervals)

# Option 2
start_list = np.array([24, 12, 8, 9, 0], dtype=np.int32)
end_list = np.array([26, 15, 12, 10, 5], dtype=np.int32)
chrom_name = 'chr1'

# arg1 - (optional) stat (default is 'mean')
# arg2 - intervals
# arg3 - start_list
# arg4 - end_list
# arg5 - chrom_name
# must have either intervals or start_list, end_list, chrom_name
# returns a numpy array of values
result = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)

correct = [-1, 0.9, 0.1, -1, 0.82]
assert len(result) == len(values)
for i in range(len(result)):
    assert result[i] == values[i]
    assert result[i] == correct[i]

# arg1 - interval file
# arg2 - (optional) output_to_file (default is True and outputs to 'chr1_out.txt'
# arg3 - (optional) stat (default is 'mean')
# returns a dictionary; keys are chromosome names, values are numpy arrays
result = bedGraph.stats_from_file('test_files/test_intervals.txt', output_to_file=False, stat='mean')

assert len(result) == 1
assert 'chr1' in result
arr = result['chr1']
for i in range(len(arr)):
    assert arr[i] == correct[i]

correct = [-1, 0.9, 0.1, -1, 0.8076923076923077]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 1/3, 0.25, 0, 1]
result = bedGraph.stats('coverage', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [-1, 0.9, 0.1, -1, 0.7]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [-1, 0.9, 0.1, -1, 0.9]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [-1, 0, 0, -1, 0.09797959]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.3, 0.025, 0, 0.82]
result = inclusive_bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.3, 0.00833333, 0, 0.7]
result = inclusive_bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.33333333, 0.25, 0, 1]
result = inclusive_bedGraph.stats('coverage', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0, 0.1, 0, 0.7]
result = inclusive_bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.9, 0.1, 0, 0.9]
result = inclusive_bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.42426407, 0.04330127, 0, 0.09797959]
result = inclusive_bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

print("Passed all simple tests!")

try:
    bedGraph = BedGraph('test_files/hg38.chrom.sizes', 'test_files/ENCFF376VCU.bigWig', 'chr14')
    assert False
except RuntimeError:
    print("Passed giving wrong chrom size test!")

start_time = time.time()
bedGraph = BedGraph('test_files/mm10.chrom.sizes', 'test_files/ENCFF376VCU.bedGraph')
print(f"Loading ENCFF376VCU.bedgraph took {time.time() - start_time}")

bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 100)

bedGraph.load_chrom_data('chr4')
bedGraph.load_chrom_bins('chr4', 100)

total_num_intervals = 0
avg_interval_sizes = {
    'chr1': 26.447609,
    'chr10': 25.53135
}

for chrom in avg_interval_sizes:
    assert abs(bedGraph.chromosome_map[chrom].avg_interval_size - avg_interval_sizes[chrom]) < 0.00001

test_intervals = [
    ['chr1', 0, 3000000],
    ['chr1', 3000000, 3100000],
    ['chr1', 3100000, 3200000],
    ['chr1', 3000297, 3000387]
]
correct = [0.0, 0.2783213414900005, 0.37037220081612465, 0.733430027961731]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.27832134, 0.3703722, 0.61754806]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 4.400169849395752, 5.133600234985352, 0.733430027961731]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.0, 0.0, 0.733430027961731]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.5458894746933545, 0.645887192684785, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

test_intervals = [
    ['chr4', 0, 50000],
    ['chr4', 3000000, 3100000],
    ['chr4', 3100000, 3200000],
    ['chr4', 3049939, 3049957]
]
correct = [0.0, 0.16491732227414846, 0.028100075026154518, 0.4027400016784668]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.16491732, 0.02810007, 0.213423]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 2.8125500679016113, 2.200079917907715, 0.4027400016784668]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.0, 0.0, 0.4027400016784668]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.3551471309944084, 0.15098161797799015, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

print("Passed all bedgraph tests!")

start_time = time.time()
bedGraph = BedGraph('test_files/mm10.chrom.sizes', 'test_files/ENCFF376VCU.bigWig')
print(f"Loading ENCFF376VCU.bigWig took {time.time() - start_time}")

bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 100)

bedGraph.load_chrom_data('chr4')
bedGraph.load_chrom_bins('chr4', 100)

total_num_intervals = 0
avg_interval_sizes = {
    'chr1': 26.447609,
    'chr10': 25.53135
}

for chrom in avg_interval_sizes:
    assert abs(bedGraph.chromosome_map[chrom].avg_interval_size - avg_interval_sizes[chrom]) < 0.00001

test_intervals = [
    ['chr1', 0, 3000000],
    ['chr1', 3000000, 3100000],
    ['chr1', 3100000, 3200000],
    ['chr1', 3000297, 3000387]
]
correct = [0.0, 0.2783213414900005, 0.37037220081612465, 0.733430027961731]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.27832134, 0.3703722, 0.61754806]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 4.400169849395752, 5.133600234985352, 0.733430027961731]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.0, 0.0, 0.733430027961731]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.5458894746933545, 0.645887192684785, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

test_intervals = [
    ['chr4', 0, 50000],
    ['chr4', 3000000, 3100000],
    ['chr4', 3100000, 3200000],
    ['chr4', 3049939, 3049957]
]
correct = [0.0, 0.16491732227414846, 0.028100075026154518, 0.4027400016784668]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.16491732, 0.02810007, 0.213423]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 2.8125500679016113, 2.200079917907715, 0.4027400016784668]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.0, 0.0, 0.4027400016784668]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0.0, 0.3551471309944084, 0.15098161797799015, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

print("Passed all bigwig tests!")
print("Passed all tests!")
