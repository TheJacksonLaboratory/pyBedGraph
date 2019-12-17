#!/bin/bash

DIR='/media/hirow/extra/jax/data'

# hg38
# files=(ENCFF321FZQ ENCFF050CCI ENCFF847JMY ENCFF631HEX)
#files=$(find "$1/hg38" -type f -name "*.bigWig")
#for file in ${files[*]}
#do
#  file=$(basename "${file%.*}")
#  echo "Running benchmark for $file"
#	python3 make_benchmarks.py "$DIR/chrom_sizes/hg38.chrom.sizes" "$1/hg38/$file.bedGraph" "$1/hg38/$file.bigWig"
#	if [[ $? -ne 0 ]];
#	then
#	    exit 1
#	fi
#done

# hg19
files=$(find "$1/hg19" -type f -name "*.bigWig")
for file in ${files[*]}
do
  file=$(basename "${file%.*}")
  echo "Running benchmark for $file"
	python3 make_benchmarks.py "$DIR/chrom_sizes/hg19.chrom.sizes" "$1/hg19/$file.bedGraph" "$1/hg19/$file.bigWig"
	if [[ $? -ne 0 ]];
	then
	    exit 1
	fi
done

# mouse
# files=(ENCFF376VCU ENCFF643WMY ENCFF384CMP ENCFF770CQD)
#files=$(find "$1/mm10" -type f -name "*.bigWig")
#for file in ${files[*]}
#do
#  file=$(basename "${file%.*}")
#  echo "Running benchmark for $file"
#	python3 make_benchmarks.py "$DIR/chrom_sizes/mm10.chrom.sizes" "$1/mm10/$file.bedGraph" "$1/mm10/$file.bigWig"
#	if [[ $? -ne 0 ]];
#	then
#	    exit 1
#	fi
#done
