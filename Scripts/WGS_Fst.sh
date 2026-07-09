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

###########################################################################################################################
## Fst calculation for all samples (with B5047-SCH-54 and B5047-SCH-60) over the whole genome
###### ATTENTION!!!!!! MeanDepth was calculated with PanDepth and the raw *rmd.bam. (just the -q from fastp was taken into account), 
###### but now the reads are filtert for -q 40 and -Q 20 from samools mpileup. MeanDepth is now different und has to be calculated from 
###### the sync files !!!!!
###### --> Man Depth over all samples and Chr = 171.6x = 172x

/home/tomsch/grenedalf/bin/grenedalf fst 
--method unbiased-hudson 
--window-type genome 
--write-pi-tables 
--sync-path /home/tomsch/WGS_36/aligned/sync_files 
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna 
--filter-sample-min-count 2 
--filter-sample-min-read-depth 86 
--filter-sample-max-read-depth 344 
--window-average-policy valid-loci 
--filter-total-snp-min-frequency 0.01 
--pool-sizes 60 
--file-prefix all_samples_fst_calculation_ 
--out-dir /home/tomsch/WGS_36/aligned/fst_files/all_samples 
--compress 
--log-file /home/tomsch/WGS_36/aligned/fst_files/all_samples/all_samples_fst.log 
--threads 20

# just for the mitochondrium 
# MeanDepth = 6836, Sd = 1900 -> so for Depth filter use mean +- 2*Sd -> 10636 - 3036
###### ATTENTION!!!!!! MeanDepth was calculated with PanDepth and the raw *rmd.bam. (just the -q from fastp was taken into account), 
###### but now the reads are filtert for -q 40 and -Q 20 from samools mpileup. MeanDepth is now different und has to be calculated from 
###### the sync files !!!!!

/home/tomsch/grenedalf/bin/grenedalf fst \
--method unbiased-hudson \
--window-type chromosomes \
--filter-region NC_001566.1:1-16343 \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/aligned/sync_files \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 3036 \
--filter-sample-max-read-depth 10636 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 30 \
--file-prefix mito_fst_all_samples \
--out-dir /home/tomsch/WGS_36/aligned/fst_files/all_samples_mito \
--compress \
--log-file /home/tomsch/WGS_36/aligned/fst_files/all_samples_mito/mito_fst_all_samples.log \
--threads 20

# Cathedral-Fst-Plot with grenedalf

/home/tomsch/grenedalf/bin/grenedalf fst-cathedral \
--method unbiased-hudson \
--filter-region NC_001566.1:1-16343 \
--sync-path /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-36.sync /home/tomsch/WGS_36/aligned/sync_files/B5047-SCH-43.sync \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 3036 \
--filter-sample-max-read-depth 10636 \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 30 \
--file-prefix 36_43_mito_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/aligned/fst_files/all_samples_mito/cathedral \
--log-file /home/tomsch/WGS_36/aligned/fst_files/all_samples_mito/cathedral/36_43_fst_cathedral.log \
--threads 20

/home/tomsch/grenedalf/bin/grenedalf cathedral-plot \
--json-path B211_mito_fst_calculation_cathedral-plot-B5047-SCH-53.1.B5047-SCH-54.1-NC_001566.1.json \
--file-prefix B211_ \
--threads 20 \
--log-file B211_cathedral_
