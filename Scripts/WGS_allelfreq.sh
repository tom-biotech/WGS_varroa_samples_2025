#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3

aligned_dir="/home/tomsch/WGS_36/aligned"
mkdir -p /home/tomsch/WGS_36/aligned/allele_freq
allel_freq_dir="/home/tomsch/WGS_36/aligned/allele_freq"
genome="/home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

for i in "$aligned_dir"/*_rmd.bam;
do dname=$(dirname ${i}); name=$(basename ${i} _rmd.bam); 
/home/tomsch/grenedalf/bin/grenedalf frequency \
--sam-path "$aligned_dir"/${name}_rmd.bam \
--sam-min-map-qual 60 \
--sam-min-base-qual 30 \
--reference-genome-fasta "$genome" \
--out-dir "$allele_freq_dir" \
--file-prefix ${name}_ \
--compress \
--threads 20 \
--log-file "$allele_freq_dir"/${name}.log \
--write-sample-counts \
--write-sample-read-depth \
--write-sample-alt-freq \
--separator-char tab; 
done

