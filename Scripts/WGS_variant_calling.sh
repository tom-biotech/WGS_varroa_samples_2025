#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 

# create all_bams.list
bam_dir="/home/Drives/HDD04_06T_SDF/Tom_WGS_Buckfast_25/fastq/Aligned_GATK"
vcf_dir="/home/Drives/HDD04_06T_SDF/Tom_WGS_Buckfast_25/fastq/vcf"

find $bam_dir -name "*.bam" > $bam_dir/all_bams.list


# Run mpileup per chromosome using parallel 
cat /home/Drives/HDD04_06T_SDF/Tom_WGS_Buckfast_25/ref/chromosomes.txt | parallel -j 40 \
    "bcftools mpileup -Ou -f /home/Drives/HDD04_06T_SDF/Tom_WGS_Buckfast_25/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
    -b $bam_dir/all_bams.list \
    -r {} -q 20 -Q 30 -d XXXXX -a FORMAT/DP,FORMAT/AD \
    -o $vcf_dir/raw_mpileup_{}.bcf"
    
## and run variant calling separately for each chromosome
for i in $(cat /home/Drives/HDD04_06T_SDF/Tom_WGS_Buckfast_25/ref/chromosomes.txt)
do name=$(basename ${i})
    echo "processing $name"
bcftools call -mv -v -Oz --threads 20 -o $vcf_dir/${name}_raw_snps.vcf.gz $vcf_dir/raw_mpileup_${name}.bcf
done

## merge all vcfs from different chromosomes
ls -1 $vcf_dir/*_raw_snps.vcf.gz > $vcf_dir/vcf_list.txt
bcftools concat -Oz $vcf_dir/vcf_list.txt -o $vcf_dir/wgs_36_raw_snp.vcf.gz
tabix -p vcf $vcf_dir/wgs_36_raw_snp.vcf.gz

## check summary information stats from vcf file
bcftools stats  -s - $vcf_dir/wgs_36_raw_snp.vcf.gz > $vcf_dir/wgs_36_raw_snp.sumstats
