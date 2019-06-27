# pyBedGraph

pyBedGraph is an alternative to pyBigWig.

# Features:

  - Finds mean, approximate mean, max, min, coverage, and standard deviation

# Improvements:
  - Much faster

# Usage:
Create the object:
```python
# Chromosome that is stored in object will be the first that appears in the bedgraph file if not given
genome = Genome('bedgraph_file.bedgraph', 'chr1')

# Specify a bin size (default is 64)
genome = Genome('bedgraph_file.bedgraph', 'chr1', 128)
```

Search from a file:
```python
# Out file is optional, output goes to stdout if not given
# Default statistic to find is mean
genome.stats_from_file('intervals_to_search_for.txt', 'out.txt')
```

Search from a list of intervals:
```python
intervals = [
    [0, 100],
    [101, 200],
    [4, 100],
    [100000, 999999]
]
print(genome.stats(intervals))
# output is [value1, value2, value3, value4]
```

Choose a specific statistic:
  - `mean`
  - `approx_mean` - an approximate mean that is around 4x faster
  - `mod_approx_mean` - a slightly slower but more accurate approximate mean
  - `max`
  - `min`
  - `coverage`
  - `std`
```python
intervals = [
    [0, 100],
    [101, 200],
    [4, 100],
    [100000, 999999]
]
genome.stats(intervals, 'approx_mean')
genome.stats_from_file('intervals_to_search_for.txt', 'out.txt', 'std')
```