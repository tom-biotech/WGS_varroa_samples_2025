#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3

aligned_dir="/home/tomsch/WGS_36/aligned"
mkdir -p /home/tomsch/WGS_36/aligned/allele_freq
allel_freq_dir="/home/tomsch/WGS_36/aligned/allele_freq"
genome="/home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

/home/tomsch/grenedalf/bin/grenedalf frequency \
--sync-path sync_files/B5047-SCH-44.sync  sync_files/B5047-SCH-53.sync  sync_files/B5047-SCH-54.sync \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--out-dir b211_allele_freq \
--file-prefix b211_samples \
--compress \
--threads 20 \
--log-file /home/tomsch/WGS_36/aligned/b211_allele_freq/b211_samples_grenedalf_freq.log \
--write-sample-counts \
--write-sample-read-depth \
--write-sample-alt-freq \
--write-total-read-depth \
--write-total-frequency \
--separator-char tab \
--multi-file-locus-set union

