import sys
sys.path.append("../..")
from pyBedGraph import BedGraph

# arg1 - chromosome sizes file
# arg2 - bedgraph file
# arg3 - (optional) chromosome_name
# Just load chromosome 'chr1' (uses less memory and takes less time)
bedGraph = BedGraph('myChrom.sizes', 'random_test.bedGraph', 'chr1')

# Load the whole bedGraph file
bedGraph = BedGraph('myChrom.sizes', 'random_test.bedGraph', 'chr1')

# Option to not ignore missing basePairs when calculating statistics
# Used the exact same way but produces slightly different results
inclusive_bedGraph = BedGraph('myChrom.sizes', 'random_test.bedGraph', ignore_missing_bp=False)

bedGraph.load_chrom_data('chr1')
inclusive_bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 3)
inclusive_bedGraph.load_chrom_bins('chr1', 3)

import numpy as np

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
value = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)
print(value)

# arg1 - interval file
# arg2 - (optional) output_to_file (default is True and outputs to 'chr1_out.txt'
# arg3 - (optional) stat (default is 'mean')
# returns a dictionary; keys are chromosome names, values are numpy arrays
result = bedGraph.stats_from_file('test_intervals.txt', output_to_file=False, stat='mean')
# {'chr1': array([-1.  ,  0.9 ,  0.1 , -1.  ,  0.82])}
print(result)

result = bedGraph.stats('mean', test_intervals)
print(result)
result = bedGraph.stats('approx_mean', test_intervals)
print(result)
result = bedGraph.stats('coverage', test_intervals)
print(result)
result = bedGraph.stats('min', test_intervals)
print(result)
result = bedGraph.stats('max', test_intervals)
print(result)
result = bedGraph.stats('std', test_intervals)
print(result)

result = inclusive_bedGraph.stats('mean', test_intervals)
print(result)
result = inclusive_bedGraph.stats('approx_mean', test_intervals)
print(result)
result = inclusive_bedGraph.stats('coverage', test_intervals)
print(result)
result = inclusive_bedGraph.stats('min', test_intervals)
print(result)
result = inclusive_bedGraph.stats('max', test_intervals)
print(result)
result = inclusive_bedGraph.stats('std', test_intervals)
print(result)
