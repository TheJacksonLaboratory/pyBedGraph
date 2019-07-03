#!/bin/bash

# human
files=(ENCFF321FZQ ENCFF050CCI ENCFF847JMY ENCFF631HEX)
for file in ${files[*]}
do
	python3 benchmark_graph.py "../pyBedGraph/data/chrom_sizes/hg38.chrom.sizes" "../pyBedGraph/data/$file.bedGraph" "../pyBedGraph/data/$file.bigWig" chr1
done

# mouse
files=(ENCFF376VCU ENCFF643WMY ENCFF384CMP ENCFF770CQD)
for file in ${files[*]}
do
	python3 benchmark_graph.py "../pyBedGraph/data/chrom_sizes/mm10.chrom.sizes" "../pyBedGraph/data/$file.bedGraph" "../pyBedGraph/data/$file.bigWig" chr1
done

# dm3 fruit fly
#files=(ENCFF651QPJ ENCFF487MJL)
#for file in ${files[*]}
#do
#	python3 benchmark_graph.py "../pyBedGraph/data/chrom_sizes/dm3.chrom.sizes" "../pyBedGraph/data/$file.bedGraph" "../pyBedGraph/data/$file.bigWig" chr3R
#done
