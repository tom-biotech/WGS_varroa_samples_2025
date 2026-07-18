#!/bin/bash

############# VARROARESISTENZ WGS 2026 ################ 
## using grenedalf 0.6.3
## using samtools 1.23.1
## using popoolation2 1.201

# directorys
bam_dir="/home/tomsch/WGS_36/aligned_new"
mpileup_dir="/home/tomsch/WGS_36/aligned_new/mpileup_files"
sync_dir="/home/tomsch/WGS_36/aligned_new/sync_files"

# create sync files from bam files

############## DOESN'T WORK ###############################
'''for i in /home/tomsch/WGS_36/aligned/*rmd.bam; \
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
'''
#############################################################

# from samtools mpileup to sync file
for i in "$bam_dir"/B5047-SCH-{25..60}_rmd.bam; do name=$(basename ${i} _rmd.bam);
cat /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt | parallel -j 20 \
"samtools mpileup -B -f /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
-q 40 -Q 20 -aa -r {} ${i} > "$mpileup_dir"/${name}_{}.mpileup"
while read c; do cat "$mpileup_dir"/${name}_${c}.mpileup; done < /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt > "$mpileup_dir"/${name}.mpileup
rm "$mpileup_dir"/${name}_N*
java -ea -Xmx10g -jar \
/home/tomsch/miniconda3/envs/WGS_36/share/popoolation2-1.201-0/mpileup2sync.jar --input "$mpileup_dir"/${name}.mpileup --output "$sync_dir"/${name}.sync --fastq-type sanger --min-qual 20 --threads 20;
done


