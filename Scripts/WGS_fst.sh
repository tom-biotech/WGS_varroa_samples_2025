#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3

# create sync files from bam files
for i in /home/tomsch/WGS_36/aligned/*rmd.bam; \
        do name=$(basename ${i} _rmd.bam); \
        /home/tomsch/grenedalf/bin/grenedalf sync \
        --reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
        --sam-path ${i} \
        --sam-min-map-qual 60 \
        --sam-min-base-qual 30 \
        --compress \
        --log-file ${name}_sync.log \
        --threads 5 \
        --file-prefix ${name}_ \
        --out-dir /home/tomsch/WGS_36/aligned/sync_files;
        done

# first of all we have to check which of the B211 samples is the true one
# so we calculate the fst for samples B11(TGW)24 (B5047-SCH-44), B211T(TGW)24 (B5047-SCH-53) and B211/2(TGW)24 
# (B5047-SCH-54). B211T(TGW)24 is used as the control,
# because we know that this is the true B211(TGW)24, but DNA was extracted from thorax instead as from heads

/home/tomsch/grenedalf/bin/grenedalf fst \
--method unbiased-hudson \
--window-type genome \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-44_sync.sync.gz /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-53_sync.sync.gz /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-54_sync.sync.gz \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 135 \
--filter-sample-max-read-depth 540 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 60 \
--file-prefix B211_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/aligned/allele_freq_2/B211_fst \
--compress \
--log-file B211_fst.log \
--threads 10


