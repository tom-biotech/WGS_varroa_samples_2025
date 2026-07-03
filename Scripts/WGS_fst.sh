#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3
## using samtools 1.23.1
## using popoolation2 1.201

# create sync files from bam files

############## DOESN'T WORK ###############################
for i in /home/tomsch/WGS_36/aligned/*rmd.bam; \
        do name=$(basename ${i} _rmd.bam); \
        /home/tomsch/grenedalf/bin/grenedalf sync \
        --reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
        --sam-path ${i} \
        --make-gapless \
        --sam-min-map-qual 60 \
        --sam-min-base-qual 30 \
        --compress \
        --log-file /home/tomsch/WGS_36/aligned/sync_files/${name}_sync.log \
        --threads 5 \
        --file-prefix ${name}_ \
        --out-dir /home/tomsch/WGS_36/aligned/sync_files;
        done
#############################################################

# from samtools mpileup to sync file
for i in B5047-SCH-{44,53,54}_rmd.bam; do name=$(basename ${i} _rmd.bam);
cat /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt | parallel -j 8 \
"samtools mpileup -B -f /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
-q 40 -Q 20 -r {} ${i} > mpileup_files/${name}_{}.mpileup";
done

while read c; do cat B5047-SCH-44_${c}.mpileup; done < /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt > B5047-SCH-44.mpileup


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


