#!/bin/bash

dir="/Volumes/guacamole/Analyses/Hippocampus/second"
for f in $dir/*; do
    if [[ -d $f ]]; then
        out=$dir
        for f2 in $f/*; do
        	fastq-dump -I --split-files -O $out $f2/*.sra
        done
    fi
done