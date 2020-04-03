import numpy as np
import pyBedGraph
from pyBedGraph import BedGraph

print(f'Using {pyBedGraph.__file__}')

try:
    bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bigWig', ['chrasdf'])
    assert False
except RuntimeError:
    print("Passed giving unknown chromosome!")

single_bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bedGraph', ['chr1'])
bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bedGraph')
bigWig = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bigWig')

single_bedGraph.load_chrom_data('chr1')
single_bedGraph.load_chrom_bins('chr1', 3)
bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 3)
bigWig.load_chrom_data('chr1')
bigWig.load_chrom_bins('chr1', 3)

inclusive_bedGraph = BedGraph('test_files/myChrom.sizes', 'test_files/random_test.bedGraph', ignore_missing_bp=False)
inclusive_bedGraph.load_chrom_data('chr1')
inclusive_bedGraph.load_chrom_bins('chr1', 3)

test_intervals = [
    ['chr1', 24, 26],
    ['chr1', 12, 15],
    ['chr1', 8, 12],
    ['chr1', 9, 10],
    ['chr1', 0, 5],
    ['chr1', 0, 30]
]
values = bedGraph.stats(intervals=test_intervals)
single_values = bedGraph.stats(intervals=test_intervals)
bigWig_values = bigWig.stats(intervals=test_intervals)

start_list = np.array([24, 12, 8, 9, 0, 0], dtype=np.int32)
end_list = np.array([26, 15, 12, 10, 5, 30], dtype=np.int32)
chrom_name = 'chr1'

result = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)

correct = [-1, 0.9, 0.1, -1, 0.82, 13/18]
assert len(result) == len(values)
for i in range(len(result)):
    assert result[i] == values[i]
    assert single_values[i] == values[i]
    assert result[i] == correct[i]
    assert abs(bigWig_values[i] - correct[i]) < 0.0000001

result = bedGraph.stats_from_file('test_files/test_intervals.txt', output_to_file=False, stat='mean')

assert len(result) == 1
assert 'chr1' in result
arr = result['chr1']
for i in range(len(arr)):
    assert arr[i] == correct[i]

correct = [-1, 0.9, 0.1, -1, 0.8076923076923077, 13/18]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 1/3, 0.25, 0, 1, 0.3]
result = bedGraph.stats('coverage', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [-1, 0.9, 0.1, -1, 0.7, 0.1]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [-1, 0.9, 0.1, -1, 0.9, 1]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [-1, 0, 0, -1, 0.09797959, 0.27799991]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.9, 0.1, 0, 4.1, 6.5]
result = bedGraph.stats('sum', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.3, 0.025, 0, 0.82, 13/60]
result = inclusive_bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.3, 0.00833333, 0, 0.7, 13/60]
result = inclusive_bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.33333333, 0.25, 0, 1, 0.3]
result = inclusive_bedGraph.stats('coverage', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0, 0.1, 0, 0.7, 0.1]
result = inclusive_bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.9, 0.1, 0, 0.9, 1]
result = inclusive_bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.42426407, 0.04330127, 0, 0.09797959, 0.36431061]
result = inclusive_bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

correct = [0, 0.9, 0.1, 0, 4.1, 6.5]
result = inclusive_bedGraph.stats('sum', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.00001

print("Passed all simple tests!")
