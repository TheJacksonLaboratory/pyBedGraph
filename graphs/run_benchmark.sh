#!/bin/bash

DIR='/media/hirwo/extra/pyBedGraph'

# human
files=(ENCFF321FZQ ENCFF050CCI ENCFF847JMY ENCFF631HEX)
for file in ${files[*]}
do
	python3 make_benchmarks.py "$DIR/chrom_sizes/hg38.chrom.sizes" "$DIR/$file.bedGraph" "$DIR/$file.bigWig"
	if [[ $? -ne 0 ]];
	then
	    exit 1
	fi
done

# mouse
files=(ENCFF376VCU ENCFF643WMY ENCFF384CMP ENCFF770CQD)
for file in ${files[*]}
do
	python3 make_benchmarks.py "$DIR/chrom_sizes/mm10.chrom.sizes" "$DIR/$file.bedGraph" "$DIR/$file.bigWig"
	if [[ $? -ne 0 ]];
	then
	    exit 1
	fi
done
