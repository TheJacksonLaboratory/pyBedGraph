# pyBedGraph
pyBedGraph is an alternative to pyBigWig for bedGraph files.

# Features:
- Finds the mean, approx. mean, max, min, coverage, or standard deviation for a given interval in a bedGraph file

# Improvements over pyBigWig:
- Much faster (>200x) for most exact statistics
- Even faster for approximate statistics

# Downsides:
- Uses much more memory
    - 16 bytes per line in bedGraph file
    - 4 bytes per basePair in every chromosome loaded
- Loading the bedGraph file takes a few minutes if it is large

# Usage:
### Create the object:
```python
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
```

### Choose and load a chromosome to search for:
```python
bedGraph.load_chrom_data('chr1')
inclusive_bedGraph.load_chrom_data('chr1')
```
### Load bins for finding mean:
For approx_mean:
1. Smaller bin size -> more accurate but slower
2. Larger bin size -> less accurate but faster
```python
bedGraph.load_chrom_bins('chr1', 3)
inclusive_bedGraph.load_chrom_bins('chr1', 3)
```
### Choose a specific statistic to search for:
  - `'mean'`
  - `'approx_mean'` - an approximate mean that is slightly faster for a 0-1% error
  - `'max'`
  - `'min'`
  - `'coverage'`
  - `'std'` - (population standard deviation)

### Search from a list of intervals:
```python
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
result = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)

# [-1.    0.9   0.1  -1.    0.82]
print(result)
```

### Search from a file:
```python
# arg1 - interval file
# arg2 - (optional) output_to_file (default is True and outputs to 'chr1_out.txt'
# arg3 - (optional) stat (default is 'mean')
# returns a dictionary; keys are chromosome names, values are numpy arrays
result = bedGraph.stats_from_file('test_intervals.txt', output_to_file=False, stat='mean')

# {'chr1': array([-1.  ,  0.9 ,  0.1 , -1.  ,  0.82])}
print(result)
```

### Sample Tests (from included test files):
```python
# [-1.    0.9   0.1  -1.    0.82]
bedGraph.stats('mean', test_intervals)

# [-1.          0.9        -1.         -1.          0.76666667]
bedGraph.stats('approx_mean', test_intervals)

# [0.         0.33333333 0.25       0.         1.        ]
bedGraph.stats('coverage', test_intervals)

# [-1.   0.9  0.1 -1.   0.7]
bedGraph.stats('min', test_intervals)

# [-1.   0.9  0.1 -1.   0.9]
bedGraph.stats('max', test_intervals)

# [-1.          0.          0.         -1.          0.09797959]
bedGraph.stats('std', test_intervals)
```

```python
# [0.    0.3   0.025 0.    0.82 ]
inclusive_bedGraph.stats('mean', test_intervals)

# [0.         0.3        0.00833333 0.         0.7       ]
inclusive_bedGraph.stats('approx_mean', test_intervals)

# [0.         0.33333333 0.25       0.         1.        ]
inclusive_bedGraph.stats('coverage', test_intervals)

# [0.  0.  0.1 0.  0.7]
inclusive_bedGraph.stats('min', test_intervals)

# [0.  0.9 0.1 0.  0.9]
inclusive_bedGraph.stats('max', test_intervals)

# [0.         0.42426407 0.04330127 0.         0.09797959]
inclusive_bedGraph.stats('std', test_intervals)
```

# Benchmark:
Actual values are found from the `stats` function in pyBigWig with the `exact` argument being `True`. The error for exact stats will be ~1e-8 due to rounding error of conversion of bigWig and bedGraph files.

Alternatively, one can make actual values be from pyBedGraph. 
```python
from pyBedGraph import Benchmark

bedGraph = BedGraph('mm10.chrom.sizes', ENCFF376VCU.bedgraph', 'chr1')
# arg1 - BedGraph object
# arg2 - bigwig file
bench = Benchmark(bedGraph, 'ENCFF376VCU.bigWig')

# arg1 - num_tests
# arg2 - interval_size
# arg3 - chrom_nam
# arg4 - bin_size
# arg5 - stats (optional) (Default is all stats)
# arg6 - just_runtime (optional) (Default is False)
# arg6 - bench_pyBigWig_approx (optional) (Default is True)
# arg6 - make_pyBigWig_baseline (optional) (Default is True)
result = bench.benchmark(10000, 500, 'chr1', 100, stats='mean')

# formatted
# mean {'run_time': 0.003580808639526367, 'error': {'percent_error': 1.1133849453411403e-08, 'ms_error': 1.1558877957200436e-15, 'abs_error': 5.565259658128112e-09, 'num_actual_0': 0}}
# pyBigWig_mean {'approx_run_time': 0.6421082019805908, 'exact_run_time': 0.6379795074462891, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}
# approx_mean {'run_time': 0.0011749267578125, 'error': {'percent_error': 0.13400600725529524, 'ms_error': 0.00964614706312478, 'abs_error': 0.068980199063462, 'num_actual_0': 0}}

# max {'run_time': 0.0027365684509277344, 'error': {'percent_error': 2.1245231544977356e-08, 'ms_error': 9.128975974031677e-13, 'abs_error': 6.218157096711807e-08, 'num_actual_0': 0}}
# pyBigWig_max {'approx_run_time': 0.6533908843994141, 'exact_run_time': 0.6436026096343994, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}

# min {'run_time': 0.002889871597290039, 'error': {'percent_error': 2.3296755440892273e-10, 'ms_error': 9.931400247350677e-19, 'abs_error': 7.883071898306948e-11, 'num_actual_0': 0}}
# pyBigWig_min {'approx_run_time': 0.6556143760681152, 'exact_run_time': 0.6390907764434814, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}

# coverage {'run_time': 0.002706289291381836, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}
# pyBigWig_coverage {'approx_run_time': 0.6487991809844971, 'exact_run_time': 0.6407179832458496, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}

# std {'run_time': 0.008781194686889648, 'error': {'percent_error': 0.0008802452423860437, 'ms_error': 3.5123006260771487e-07, 'abs_error': 0.0004987475752671237, 'num_actual_0': 0}}
# pyBigWig_std {'approx_run_time': 0.6418542861938477, 'exact_run_time': 0.6490097045898438, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}

# Unformatted
# {'mean': {'run_time': 0.003580808639526367, 'error': {'percent_error': 1.1133849453411403e-08, 'ms_error': 1.1558877957200436e-15, 'abs_error': 5.565259658128112e-09, 'num_actual_0': 0}}, 'pyBigWig_mean': {'approx_run_time': 0.6421082019805908, 'exact_run_time': 0.6379795074462891, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}, 'approx_mean': {'run_time': 0.0011749267578125, 'error': {'percent_error': 0.13400600725529524, 'ms_error': 0.00964614706312478, 'abs_error': 0.068980199063462, 'num_actual_0': 0}}, 'max': {'run_time': 0.0027365684509277344, 'error': {'percent_error': 2.1245231544977356e-08, 'ms_error': 9.128975974031677e-13, 'abs_error': 6.218157096711807e-08, 'num_actual_0': 0}}, 'pyBigWig_max': {'approx_run_time': 0.6533908843994141, 'exact_run_time': 0.6436026096343994, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}, 'min': {'run_time': 0.002889871597290039, 'error': {'percent_error': 2.3296755440892273e-10, 'ms_error': 9.931400247350677e-19, 'abs_error': 7.883071898306948e-11, 'num_actual_0': 0}}, 'pyBigWig_min': {'approx_run_time': 0.6556143760681152, 'exact_run_time': 0.6390907764434814, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}, 'coverage': {'run_time': 0.002706289291381836, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}, 'pyBigWig_coverage': {'approx_run_time': 0.6487991809844971, 'exact_run_time': 0.6407179832458496, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}, 'std': {'run_time': 0.008781194686889648, 'error': {'percent_error': 0.0008802452423860437, 'ms_error': 3.5123006260771487e-07, 'abs_error': 0.0004987475752671237, 'num_actual_0': 0}}, 'pyBigWig_std': {'approx_run_time': 0.6418542861938477, 'exact_run_time': 0.6490097045898438, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'num_actual_0': 0}}}
```
