# pyBedGraph

pyBedGraph is an alternative to pyBigWig.

# Features:

  - Finds mean, approximate mean, median, max, min, coverage, and standard deviation

# Improvements:
  - Much faster

# Usage:
Create the object:
In general, a smaller bin size (relative to the interval size) will be more accurate but slower compared to a larger bin size.
```python
from pyBedGraph import BedGraph

# arg1 - chromosome sizes file
# arg2 - bedgraph file
# arg3 - (optional) chromosome_name
# Just load chromosome 'chr1' (uses less memory and takes less time)
bedGraph = BedGraph('dm6.chrom.sizes', 'bedgraph_file.bedgraph', 'chr1')

# Load the whole bedGraph file
bedGraph = BedGraph('dm6.chrom.sizes', 'bedgraph_file.bedgraph')
```
Choose a specific statistic:
  - `'mean'`
  - `'approx_mean'` - an approximate mean that is around 5x faster for a 0-1% error
  - `'median'`
  - `'max'`
  - `'min'`
  - `'coverage'`
  - `'std'`

Search from a file: (TODO)
```python
# arg1 - interval file
# arg2 - (optional) output_file (default goes to stdout)
# arg3 - (optional) stat (default is 'mean')
bedGraph.stats_from_file('intervals_to_search_for.txt', output_file='out.txt', stat='mean')
```

Search from a list of intervals:
```python
# Option 1
intervals = [
    ['chr2L', 0, 100],
    ['chr2L', 101, 200],
    ['chr2L', 4, 100],
    ['chr2L', 100000, 999999]
]

# Option 2
start_list = [0, 101, 4, 100000]
end_list = [100, 200, 100, 999999]
chrom_name = 'chr2L'

# arg1 - (optional) stat (default is 'mean')
# arg2 - intervals
# arg3 - start_list
# arg4 - end_list
# arg5 - chrom_name
# must have either intervals or start_list, end_list, chrom_name
# returns a list of values
values = bedGraph.stats(intervals=intervals)

values = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)

# Output is [value1, value2, value3, value4]
print(values)
```



# Benchmark:
Actual values are found from the `stats` function in pyBigWig with the `exact` argument being `True`. (TODO)
```python
from pyBedGraph import Benchmark

# arg1 - BedGraph object
# arg2 - bigwig file
bedGraph = BedGraph('P2MC7N8HCE3K.bedgraph')
bench = Benchmark(bedGraph, 'P2MC7N8HCE3K.bw')
result = bench.benchmark(num_test, interval_size, chrom_name, bin_size, stats=None)
```