import sys
import math
sys.path.append("..")
from pyBedGraph.BedGraph import BedGraph
import pyBedGraph

average_interval_size = 5000
num_tests = 10000
chrom_name = 'chr1'
bin_size = int(math.sqrt(average_interval_size))
bin_size = 1000
stats = ['max_index']

print(pyBedGraph.__file__)

bedGraph = BedGraph('mm10.chrom.sizes', 'ENCFF376VCU.bigWig', 'chr1')
# bedGraph = BedGraph('mm10.chrom.sizes', 'ENCFF376VCU.bedGraph', 'chr1')
