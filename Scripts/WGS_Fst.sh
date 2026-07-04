#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3

# first of all we have to check which of the B211 samples is the true one
# so we calculate the fst for samples B11(TGW)24 (B5047-SCH-44), B211T(TGW)24 (B5047-SCH-53) and B211/2(TGW)24 
# (B5047-SCH-54). B211T(TGW)24 is used as the control,
# because we know that this is the true B211(TGW)24, but DNA was extracted from thorax instead as from heads

/home/tomsch/grenedalf/bin/grenedalf fst \
--method unbiased-hudson \
--window-type genome \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-44.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-53.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-54.sync \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 135 \
--filter-sample-max-read-depth 540 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 60 \
--file-prefix B211_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/aligned/fst_files/B211_fst \
--compress \
--log-file /home/tomsch/WGS_36/aligned/fst_files/B211_fst/B211_fst.log \
--threads 20

# just for the mitochondrium 
# MeanDepth = 24199, Sd = 5010 -> so for Depth filter use mean +- 2*Sd -> 34219 - 14179
###### ATTENTION!!!!!! MeanDepth was calculated with PanDepth and the raw *rmd.bam. (just the -q from fastp was taken into account), 
###### but now the reads are filtert for -q 40 and -Q 20 from samools mpileup. MeanDepth is now different und has to be calculated from 
###### the sync files !!!!!



/home/tomsch/grenedalf/bin/grenedalf fst \
--method unbiased-hudson \
--window-type chromosomes \
--filter-region NC_001566.1:1-16343 \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-44.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-53.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-54.sync \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 14179 \
--filter-sample-max-read-depth 34219 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 30 \
--file-prefix B211_mito_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/aligned/fst_files/B211_mito_fst \
--compress \
--log-file /home/tomsch/WGS_36/aligned/fst_files/B211_mito_fst/B211_mito_fst.log \
--threads 20


# Cathedral-Fst-Plot with grenedalf

/home/tomsch/grenedalf/bin/grenedalf fst-cathedral \
--method unbiased-hudson \
--window-type chromosomes \
--window-region NC_001566.1 \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-44.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-53.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-54.sync \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 14179 \
--filter-sample-max-read-depth 34219 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 30 \
--file-prefix B211_mito_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/aligned/fst_files/B211_mito_fst \
--compress \
--log-file /home/tomsch/WGS_36/aligned/fst_files/B211_mito_fst/B211_mito_fst.log \
--threads 20
