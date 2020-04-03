import time
import pyBedGraph
from pyBedGraph import BedGraph

print(f'Using {pyBedGraph.__file__}')

DEBUG = False

try:
    bedGraph = BedGraph('test_files/hg38.chrom.sizes', 'test_files/ENCFF376VCU.bigWig', ['chr14'])
    assert False
except RuntimeError:
    print("Passed giving wrong chrom size test!")

start_time = time.time()
bedGraph = BedGraph('test_files/mm10.chrom.sizes', 'test_files/ENCFF376VCU.bedGraph', debug=DEBUG)
print(f"Loading ENCFF376VCU.bedgraph took {time.time() - start_time}")
# Takes 170 seconds on i5-7300HQ

bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 100)

bedGraph.load_chrom_data('chr4')
bedGraph.load_chrom_bins('chr4', 100)

if DEBUG:
    total_num_intervals = 0
    avg_interval_sizes = {
        'chr1': 26.447609,
        'chr10': 25.53135
    }

    for chrom in avg_interval_sizes:
        assert abs(bedGraph.chromosome_map[chrom].avg_interval_size - avg_interval_sizes[chrom]) < 0.0001

test_intervals = [
    ['chr1', 0, 3000000],
    ['chr1', 3000000, 3100000],
    ['chr1', 3100000, 3200000],
    ['chr1', 3000297, 3000387]
]
correct = [0.0, 0.2783213414900005, 0.37037220081612465, 0.733430027961731]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0, 0.27832134, 0.3703722, 0.61754806]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 4.400169849395752, 5.133600234985352, 0.733430027961731]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.0, 0.0, 0.733430027961731]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.5458894746933545, 0.645887192684785, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 27832.1338, 37037.21976, 66.0087]
result = bedGraph.stats('sum', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

test_intervals = [
    ['chr4', 0, 50000],
    ['chr4', 3000000, 3100000],
    ['chr4', 3100000, 3200000],
    ['chr4', 3049939, 3049957]
]
correct = [0.0, 0.16491732227414846, 0.028100075026154518, 0.4027400016784668]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0, 0.16491732, 0.02810007, 0.213423]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 2.8125500679016113, 2.200079917907715, 0.4027400016784668]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.0, 0.0, 0.4027400016784668]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.3551471309944084, 0.15098161797799015, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 16491.7322, 2810.00742, 7.24932]
result = bedGraph.stats('sum', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

print("Passed all bedgraph tests!")

start_time = time.time()
bedGraph = BedGraph('test_files/mm10.chrom.sizes', 'test_files/ENCFF376VCU.bigWig', debug=DEBUG)
print(f"Loading ENCFF376VCU.bigWig took {time.time() - start_time}")
# Takes 69 seconds on i5-7300HQ

bedGraph.load_chrom_data('chr1')
bedGraph.load_chrom_bins('chr1', 100)

bedGraph.load_chrom_data('chr4')
bedGraph.load_chrom_bins('chr4', 100)

if DEBUG:
    avg_interval_sizes = {
        'chr1': 26.447609,
        'chr10': 25.53135
    }

    for chrom in avg_interval_sizes:
        assert abs(bedGraph.chromosome_map[chrom].avg_interval_size - avg_interval_sizes[chrom]) < 0.0001

test_intervals = [
    ['chr1', 0, 3000000],
    ['chr1', 3000000, 3100000],
    ['chr1', 3100000, 3200000],
    ['chr1', 3000297, 3000387]
]
correct = [0.0, 0.2783213414900005, 0.37037220081612465, 0.733430027961731]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0, 0.27832134, 0.3703722, 0.61754806]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 4.400169849395752, 5.133600234985352, 0.733430027961731]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.0, 0.0, 0.733430027961731]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.5458894746933545, 0.645887192684785, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

# LOWERED STANDARD FOR ROUNDING PROBLEMS
correct = [0.0, 27832.1338, 37037.21976, 66.0087]
result = bedGraph.stats('sum', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.1

test_intervals = [
    ['chr4', 0, 50000],
    ['chr4', 3000000, 3100000],
    ['chr4', 3100000, 3200000],
    ['chr4', 3049939, 3049957]
]
correct = [0.0, 0.16491732227414846, 0.028100075026154518, 0.4027400016784668]
result = bedGraph.stats('mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0, 0.16491732, 0.02810007, 0.213423]
result = bedGraph.stats('approx_mean', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 2.8125500679016113, 2.200079917907715, 0.4027400016784668]
result = bedGraph.stats('max', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.0, 0.0, 0.4027400016784668]
result = bedGraph.stats('min', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 0.3551471309944084, 0.15098161797799015, 0.0]
result = bedGraph.stats('std', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

correct = [0.0, 16491.7322, 2810.00742, 7.24932]
result = bedGraph.stats('sum', test_intervals)
for i in range(len(result)):
    assert abs(result[i] - correct[i]) < 0.0001

print("Passed all bigwig tests!")
print("Passed all extensive tests!")
