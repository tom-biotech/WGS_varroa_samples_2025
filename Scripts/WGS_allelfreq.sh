#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3

aligned_dir="/home/tomsch/WGS_36/aligned"
mkdir -p /home/tomsch/WGS_36/aligned/allele_freq
allel_freq_dir="/home/tomsch/WGS_36/aligned/allele_freq"
genome="/home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

/home/tomsch/grenedalf/bin/grenedalf frequency \
--sam-path B5047-SCH-25_rmd.bam  B5047-SCH-29_rmd.bam  B5047-SCH-33_rmd.bam  B5047-SCH-37_rmd.bam  B5047-SCH-41_rmd.bam  B5047-SCH-45_rmd.bam  B5047-SCH-49_rmd.bam  B5047-SCH-53_rmd.bam  B5047-SCH-57_rmd.bam B5047-SCH-26_rmd.bam  B5047-SCH-30_rmd.bam  B5047-SCH-34_rmd.bam  B5047-SCH-38_rmd.bam  B5047-SCH-42_rmd.bam  B5047-SCH-46_rmd.bam  B5047-SCH-50_rmd.bam  B5047-SCH-54_rmd.bam  B5047-SCH-58_rmd.bam B5047-SCH-27_rmd.bam  B5047-SCH-31_rmd.bam  B5047-SCH-35_rmd.bam  B5047-SCH-39_rmd.bam  B5047-SCH-43_rmd.bam  B5047-SCH-47_rmd.bam  B5047-SCH-51_rmd.bam  B5047-SCH-55_rmd.bam  B5047-SCH-59_rmd.bam B5047-SCH-28_rmd.bam  B5047-SCH-32_rmd.bam  B5047-SCH-36_rmd.bam  B5047-SCH-40_rmd.bam  B5047-SCH-44_rmd.bam  B5047-SCH-48_rmd.bam  B5047-SCH-52_rmd.bam  B5047-SCH-56_rmd.bam  B5047-SCH-60_rmd.bam  \
--sam-min-map-qual 60 \
--sam-min-base-qual 30 \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--out-dir allele_freq_2 \
--file-prefix all_samples \
--compress \
--threads 15 \
--log-file /home/tomsch/WGS_36/aligned/allele_freq_2/all_samples_grenedalf_freq.log \
--write-sample-counts \
--write-sample-read-depth \
--write-sample-alt-freq \
--write-total-read-depth \
--write-total-frequency \
--separator-char tab \
--multi-file-locus-set union

