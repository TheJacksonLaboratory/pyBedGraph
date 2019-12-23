[![Build Status](https://travis-ci.org/c0ver/pyBedGraph.svg?branch=master)](https://travis-ci.org/c0ver/pyBedGraph)

# pyBedGraph
A Python package for fast operations on 1-dimensional genomic signal tracks.

## Features
- Finds the mean, approx. mean, max, min, coverage, or standard deviation for a given interval in a bedGraph file
- Partly written in Cython for speed improvements
- Can look up exact statistics of 1 million regions in ~0.26 second on a conventional laptop
- An approximate mean for 10,000 regions can be computed in ~0.0012 second w/ error rate of less than 5 percent

## Drawbacks
- Uses memory to load files
    - 16 bytes per line in bedGraph file
    - 4 bytes per basePair in every chromosome loaded
- Loading the bedGraph file can take up to a minute or two
- Only works with sorted bedgraph files

## Installation

Dependency requirements:
- Numpy >= v1.16.4
- pyBigWig >= v0.3.16 (For reading bigWig files)
    - pyBigWig == 0.3.16 (For Benchmarking)

With pip:
```bash
pip3 install pyBedGraph
pip3 install pyBigWig # if using bigwig files
```

With conda:
```bash
conda create -n test
conda activate test
conda install -c bioconda pyBedGraph
conda install -c bioconda pyBigWig # if using bigwig files
```

## Usage

### Download the test files here:
https://thejacksonlaboratory.ent.box.com/s/3jglutwf3d54pnomnp33ivo7a9546vhe

Test files are also included in this Github repository: `test/test_files`.

Enter the directory with the test files.

### Create the object:
```python
from pyBedGraph import BedGraph

# arg1 - chromosome sizes file
# arg2 - bedgraph file
# arg3 - (optional) chromosome_name
# Just load chromosome 'chr1' (uses less memory and takes less time)
bedGraph = BedGraph('myChrom.sizes', 'random_test.bedGraph', 'chr1')

# Load the whole bedGraph file
bedGraph = BedGraph('myChrom.sizes', 'random_test.bedGraph')

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
  - `'approx_mean'` - an approximate mean is faster than exact mean, with < 5% error rate
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

# [-1.    0.9   0.1  -1.    0.82]
print(values)

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

# [-1.          0.9        0.1.         -1.          0.8076923076923077]
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

## Benchmarking pyBedGraph
Actual values are found from the `stats` function in pyBigWig with the `exact` argument being `True`. The error for exact stats will be ~1e-8 due to rounding error of conversion of bigWig and bedGraph files.

Alternatively, one can make actual values be pyBedGraph's exact statistics.

Enter the `graphs` folder in the Github project repository.
```python
from pyBedGraph import BedGraph
from Benchmark import Benchmark

# These files can be downloaded from the link given above
bedGraph = BedGraph('mm10.chrom.sizes', 'ENCFF376VCU.bedGraph', 'chr1')

# Alternatively using a bigwig file
# bedGraph = BedGraph('mm10.chrom.sizes', 'ENCFF376VCU.bigWig', 'chr1')

bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 100)

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
# Test all statistics
result = bench.benchmark(10000, 5000, 'chr1', 100)

for key in result:
    print(key, result[key])

# mean {'run_time': 0.008324861526489258, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'not_included': 0}}
# pyBigWig_mean {'approx_run_time': 1.4333949089050293, 'exact_run_time': 0.7698564529418945, 'error': {'percent_error': 0.06567272540694802, 'ms_error': 0.001222419386871348, 'abs_error': 0.023540340949669364, 'not_included': 79}}
# approx_mean {'run_time': 0.002111673355102539, 'error': {'percent_error': 0.006529644707171326, 'ms_error': 7.858080037556034e-06, 'abs_error': 0.001824641073039555, 'not_included': 4}}
# max {'run_time': 0.005040645599365234, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'not_included': 0}}
# pyBigWig_max {'approx_run_time': 1.2673799991607666, 'exact_run_time': 0.7933700084686279, 'error': {'percent_error': 0.10220448242023446, 'ms_error': 1.2678718593032368, 'abs_error': 0.25865022624731066, 'not_included': 79}}
# min {'run_time': 0.005083560943603516, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'not_included': 0}}
# pyBigWig_min {'approx_run_time': 1.2120039463043213, 'exact_run_time': 0.7468140125274658, 'error': {'percent_error': 0.0001, 'ms_error': 7.109862619931795e-07, 'abs_error': 8.432000130414962e-06, 'not_included': 0}}
# coverage {'run_time': 0.0063626766204833984, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'not_included': 0}}
# pyBigWig_coverage {'approx_run_time': 1.2101118564605713, 'exact_run_time': 0.7483360767364502, 'error': {'percent_error': 0.0, 'ms_error': 0.0, 'abs_error': 0.0, 'not_included': 0}}
# std {'run_time': 0.0422673225402832, 'error': {'percent_error': 9.690484548456011e-05, 'ms_error': 4.764358150024449e-09, 'abs_error': 6.25265457158463e-05, 'not_included': 0}}
# pyBigWig_std {'approx_run_time': 1.219078540802002, 'exact_run_time': 0.7484426498413086, 'error': {'percent_error': 0.04560011737269686, 'ms_error': 0.005008324729263816, 'abs_error': 0.02569405301725115, 'not_included': 79}}
```

## Testing pyBedGraph
Some tests are provided in `test/test.py`. Additional bedgraph and bigwig files for ENCFF376VCU are needed to run extensive_test.py. Build badge comes from a forked repository, [https://github.com/c0ver/pyBedGraph](https://github.com/c0ver/pyBedGraph), that has the same version as this repository.

## Reference 
[pyBedGraph: a Python package for fast operations on 1-dimensional genomic signal tracks](https://www.biorxiv.org/content/10.1101/709683v1), Zhang et al., bioRxiv, 2019

## Bug reports
To report bugs, contact Henry (henrybzhang.99@gmail.com) and Minji (minji.kim@jax.org) or visit the [Issues](https://github.com/TheJacksonLaboratory/pyBedGraph/issues) page. 

## Potential Improvements
Make index array use negative numbers to show index of closest interval instead of just -1. This will help for datasets that have low coverage (RNA-seq)

Use index dictionary instead of index array if there is a memory efficient option.

Store intervals differently as a space between two intervals and length of each interval instead of just the start and end. This might allow usage of 2 byte shorts.
