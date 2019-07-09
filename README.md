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
bedGraph = BedGraph('hg38.chrom.sizes', 'ENCFF321FZQ.bedGraph', 'chr1')

# Load the whole bedGraph file
bedGraph = BedGraph('hg38.chrom.sizes', 'ENCFF321FZQ.bedGraph')
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
    ['chr1', 0, 1000],
    ['chr1', 1001, 1500],
    ['chr1', 2000, 2200],
    ['chr1', 3000, 5000],
    ['chr1', 5001, 10000],
    ['chr1', 1000000, 101000]
]

# Option 2
start_list = [0, 1001, 2000, 3000, 5001, 1000000]
end_list = [1000, 1500, 2200, 5000, 10000, 101000]
chrom_name = 'chr1'

# arg1 - (optional) stat (default is 'mean')
# arg2 - intervals
# arg3 - start_list
# arg4 - end_list
# arg5 - chrom_name
# must have either intervals or start_list, end_list, chrom_name
# returns a list of values
values = bedGraph.stats(intervals=intervals)

values = bedGraph.stats(start_list=start_list, end_list=end_list, chrom_name=chrom_name)

# Output is [0.         0.         0.         0.         0.00207475 0.05981362]
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