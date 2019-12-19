#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "You must enter exactly 2 command line arguments:

	1 - Folder location of bedgraph files
	2 - chrom size file abbreviation (ex: mm10, hg38)"
	exit 1
fi

chrom_size="/media/hirow/extra/jax/data/chrom_sizes/$2.chrom.sizes"
if [ ! -f $chrom_size ]; then
	echo "$chrom_size does not exist"
	exit 1
fi

files=$(find "$1/$2" -type f -name "*.bigWig")
echo "Number of files found in \"$1/$2\": $(ls -1 "$1"/"$2" | wc -l)"

for file in ${files[*]}
do
  file=$(basename "${file%.*}")
  echo "Running benchmark for $file"
	python3 make_benchmarks.py "$chrom_size" "$1/$2/$file.bigWig"
	if [[ $? -ne 0 ]];
	then
	    exit 1
	fi
done
