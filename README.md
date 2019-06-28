# pyBedGraph

pyBedGraph is an alternative to pyBigWig.

# Features:

  - Finds mean, approximate mean, max, min, coverage, and standard deviation

# Improvements:
  - Much faster
  - Avoids rounding errors found in pyBigWig

# Usage:
Create the object:
In general, a smaller bin size (relative to the interval size) will be more accurate but slower compared to a larger bin size. Also, bin sizes should be powers of two.
```python
# arg1 - bedgraph file
# arg2 - (optional) chromosome_name
# arg3 - (optional) bin_size (default is 64)
genome = Genome('bedgraph_file.bedgraph', 'chr1')

genome = Genome('bedgraph_file.bedgraph', chromosome_name='chr1', bin_size=128)
```

Search from a file:
```python
# arg1 - interval file
# arg2 - (optional) output_file (default goes to stdout)
# arg3 - (optional) stat (default is 'mean')
genome.stats_from_file('intervals_to_search_for.txt', output_file='out.txt', stat='mean')
```

Search from a list of intervals:
```python
intervals = [
    [0, 100],
    [101, 200],
    [4, 100],
    [100000, 999999]
]

# arg1 - interval list
# arg2 - (optional) stat (default is 'mean')
print(genome.stats(intervals))
# output is [value1, value2, value3, value4]
```

Choose a specific statistic:
  - `'mean'`
  - `'approx_mean'` - an approximate mean that is around 2x faster (with default bin size)
  - `'mod_approx_mean'` - a slightly slower but more accurate approximate mean
  - `'max'`
  - `'min'`
  - `'coverage'`
  - `'std'`
```python
genome.stats_from_file('intervals_to_search_for.txt', 'out.txt', 'std')
```

# Benchmark:
Actual values are found from the `intervals` function in pyBigWig. Times are compared to pyBigWig's `stat` function.
```python
# arg1 - genome object
# arg2 - number of tests
# arg3 - interval size
# arg4 - bigwig file
# arg5 - (optional) test (not working)
genome = Genome('P2MC7N8HCE3K.bedgraph')
Benchmark(genome, 10000, 128, 'P2MC7N8HCE3K.bw')
```