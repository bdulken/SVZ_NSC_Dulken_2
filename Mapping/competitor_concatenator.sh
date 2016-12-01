#!/bin/bash

dir="/Volumes/guacamole/Analyses/Competitors_Dataset/SRX993/"
for f in $dir/*; do
    if [[ -d $f ]]; then
    	base=$(basename $f)
		cat $f/*_1.fastq > "${dir}/${base}_1.fastq"
		cat $f/*_2.fastq > "${dir}/${base}_2.fastq"
		rm -f $f/*_1.fastq
		rm -f $f/*_2.fastq
    fi
done
