# pyBedGraph

pyBedGraph is an alternative to pyBigWig.

# Features:

  - Finds mean, approximate mean, max, min, coverage, and standard deviation

# Improvements:
  - Much faster

# Usage:
Create the object:
In general, a smaller bin size (relative to the interval size) will be more accurate but slower compared to a larger bin size. Also, bin sizes should be powers of two.
```python
from pyBedGraph import BedGraph

# arg1 - chromosome sizes file
# arg2 - bedgraph file
# arg3 - (optional) chromosome_name
# arg4 - (optional) bin_size (default is 64)
bedGraph = BedGraph('dm6.chrom.sizes','bedgraph_file.bedgraph', 'chr1')

bedGraph = bedGraph('bedgraph_file.bedgraph', chromosome_name='chr1', bin_size=128)
```

Search from a file:
```python
# arg1 - interval file
# arg2 - (optional) output_file (default goes to stdout)
# arg3 - (optional) stat (default is 'mean')
bedGraph.stats_from_file('intervals_to_search_for.txt', output_file='out.txt', stat='mean')
```

Search from a list of intervals:
```python
intervals = [
    ['chr2L', 0, 100],
    ['chr2L', 101, 200],
    ['chr2L', 4, 100],
    ['chr2L', 100000, 999999]
]

# arg1 - interval list
# arg2 - (optional) stat (default is 'mean')
print(bedGraph.stats(intervals))
# output is [value1, value2, value3, value4]
```

Choose a specific statistic:
  - `'mean'`
  - `'approx_mean'` - an approximate mean that is around 5x faster
  - `'mod_approx_mean'` - a slightly slower but almost halves the error in `'approx_mean'`
  - `'max'`
  - `'min'`
  - `'coverage'`
  - `'std'`
```python
bedGraph.stats_from_file('intervals_to_search_for.txt', 'out.txt', 'std')
```

# Benchmark:
Actual values are found from the `intervals` function in pyBigWig. Times are compared to pyBigWig's `stat` function with the `exact` argument being `True`.
```python
from pyBedGraph import Benchmark

# arg1 - BedGraph object
# arg2 - number of tests
# arg3 - interval size
# arg4 - bigwig file
# arg5 - (optional) test (not working)
bedGraph = BedGraph('P2MC7N8HCE3K.bedgraph')
Benchmark(bedGraph, 10000, 128, 'P2MC7N8HCE3K.bw')
```