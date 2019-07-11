# pyBedGraph
pyBedGraph is an alternative to pyBigWig for bedGraph files.

# Features:
- Finds the mean, approximate mean, median, max, min, coverage, or standard deviation for an interval

# Improvements:
- Much faster (>100x)

# Downsides:
- Uses much more memory
    - 16 bytes per line in bedGraph file
    - 8 bytes per basePair in chromosome loaded
- Loading the bedGraph file takes a few minutes if it is large

# Usage:
### Create the object:
```python
from pyBedGraph import BedGraph

# arg1 - chromosome sizes file
# arg2 - bedgraph file
# arg3 - (optional) chromosome_name
# Just load chromosome 'chr1' (uses less memory and takes less time)
bedGraph = BedGraph('hg38.chrom.sizes', 'ENCFF321FZQ.bedGraph', 'chr1')

# Load the whole bedGraph file
bedGraph = BedGraph('hg38.chrom.sizes', 'ENCFF321FZQ.bedGraph')

# Option to not ignore missing basePairs when calculating statistics
bedGraph = BedGraph('hg38.chrom.sizes', 'ENCFF321FZQ.bedGraph', ignore_missing_bp=False)
```

### Choose a specific statistic:
  - `'mean'`
  - `'approx_mean'` - an approximate mean that is around 5x faster for a 0-1% error
  - `'median'` - Much slower (>>10x) due to not being implemented in Cython
  - `'max'`
  - `'min'`
  - `'coverage'`
  - `'std'` - 10x slower than the rest of the stats

### Choose and load a chromosome to search for:
```python
bedGraph.load_chrom_data('chr1')
```
### Load bins for finding mean:
For mean:
- sqrt(interval size)

For approx_mean:
1. Smaller bin size -> more accurate but slower
2. Larger bin size -> less accurate but faster
```python
bedGraph.load_chrom_bins('chr1', 100)
```

### Search from a file:
```python
# arg1 - interval file
# arg2 - (optional) output_to_file (default is True and outputs to 'chr1_out.txt'
# arg3 - (optional) stat (default is 'mean')
# returns a dictionary; keys are chromosome names, values are numpy arrays
result = bedGraph.stats_from_file('intervals_to_search_for.txt', output_to_file=False, stat='mean')
```

### Search from a list of intervals:
```python
import numpy as np

# Option 1
intervals = [
    ['chr1', 0, 1000],
    ['chr1', 1001, 1500],
    ['chr1', 2000, 2200],
    ['chr1', 3000, 5000],
    ['chr1', 5001, 10000],
    ['chr1', 100000, 101000]
]

# Option 2
start_list = np.array([0, 1001, 2000, 3000, 5001, 100000], dtype=np.int32)
end_list = np.array([1000, 1500, 2200, 5000, 10000, 101000], dtype=np.int32)
chrom_name = 'chr1'

# arg1 - (optional) stat (default is 'mean')
# arg2 - intervals
# arg3 - start_list
# arg4 - end_list
# arg5 - chrom_name
# must have either intervals or start_list, end_list, chrom_name
# returns a numpy array of values
values = bedGraph.stats(intervals=intervals)

values = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)

# Output is [0.         0.         0.         0.         0.00207475 0.05981362]
print(values)
```

# Benchmark:
Actual values are found from the `stats` function in pyBigWig with the `exact` argument being `True`. The error for exact stats will be ~1e-8 due to rounding error of conversion of bigWig and bedGraph files.

Alternatively, one can make actual values be from pyBedGraph. 
```python
from pyBedGraph import Benchmark

# arg1 - BedGraph object
# arg2 - bigwig file
bedGraph = BedGraph('P2MC7N8HCE3K.bedgraph')
bench = Benchmark(bedGraph, 'P2MC7N8HCE3K.bw')

# arg1 - num_tests
# arg2 - interval_size
# arg3 - chrom_nam
# arg4 - bin_size
# arg5 - stats (optional) (Default is all stats)
# arg6 - only_runtime (optional) (Default is False)
# arg6 - bench_pyBigWig (optional) (Default is True)
# arg6 - pyBigWig_baseline (optional) (Default is True)
result = bench.benchmark(num_test, interval_size, chrom_name, bin_size, stats='mean')
```
