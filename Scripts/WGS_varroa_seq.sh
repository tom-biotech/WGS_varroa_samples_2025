#!/bin/bash

############# VARROARESISTENZ WGS 2026 ################ 
#### using bwa-mem2 2.3
#### using samtools 1.23.1
#### using pandepth 2.26
#### using alfred 0.5.3


bam_dir="/home/tomsch/WGS_36/aligned"
mkdir -p "/home/tomsch/WGS_36/unmapped"
unmap_dir="/home/tomsch/WGS_36/unmapped"

# using the unmapped reads from mapping against the Amel_HAv3.1 to mapp against mito genome of Varroa destructor

# extract unmapped reads with samtools

for i in $aligned/*_rmd.bam;
do name=$(basename $i _rmd.bam);

samtools view -b -f 4 -@ 20 $i -o $unmap_dir/${name}_unmapped.bam

done
