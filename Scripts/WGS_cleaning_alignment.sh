
#!/bin/sh

############# VARROARESISTENZ WGS 2026 ################ 
#### using bbmap 39.85
#### using bwa-mem2 2.3
#### using samtools 1.23.1

fastq_raw="/drives/HDD_22TB_RAWDATA/20260605_Genewiz_Varroaresistenz/00_fastq"
mkdir -p "/home/tomsch/WGS_36/Clean"
fastq_clean="/home/tomsch/WGS_36/Clean"
genome="/home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
mkdir -p "/home/tomsch/WGS_36/QC"
qc_dir="/home/tomsch/WGS_36/QC"
prop_aligned="$qc_dir/counts_summary_all_samples.txt"
mkdir -p "/home/tomsch/WGS_36/aligned"
aligned_dir="/home/tomsch/WGS_36/aligned"

################# just for the first 5 samples ############################

for name in B5047-SCH-30 B5047-SCH-31 B5047-SCH-32 B5047-SCH-33 B5047-SCH-34; do
    fastp -i $fastq_raw/${name}_R1_001.fastq.gz \
          -I $fastq_raw/${name}_R2_001.fastq.gz \
          -o $fastq_clean/${name}_clean_1.fastq.gz \
          -O $fastq_clean/${name}_clean_2.fastq.gz \
          -h $fastq_clean/${name}_stats.html \
          -j $fastq_clean/${name}_stats.json \
          -w 10 --dont_overwrite --detect_adapter_for_pe \
          -q 20 -l 30 --correction --cut_right --cut_mean_quality 20
done

####################### fuer alle Proben ##################################
#for i in "$fastq_raw"/*_R1_001.fastq.gz; \
#        do dname=$(dirname ${i}); name=$(basename ${i} _R1_001.fastq.gz); \
#        fastp -i ${dname}/${name}_R1_001.fastq.gz -I ${dname}/${name}_R2_001.fastq.gz -o $fastq_clean/${name}_clean_1.fastq.gz -O $fastq_clean/${name}_clean_2.fastq.gz -h $fastq_clean/#${name}_stats.html -j $fastq_clean/${name}_stats.json -w 10 --dont_overwrite --detect_adapter_for_pe -q 20 -l 30 --correction --cut_right --cut_mean_quality 20 ;
#done

## call the genome (using the curated GCF version), to download here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003254395.2/

## index genome
samtools faidx $genome
bwa-mem2 index $genome

for i in "$fastq_clean"/B5047-SCH-30_clean_1.fastq.gz B5047-SCH-31_clean_1.fastq.gz B5047-SCH-32_clean_1.fastq.gz B5047-SCH-33_clean_1.fastq.gz B5047-SCH-34_clean_1.fastq.gz
        do dname=$(dirname "${i}"); name=$(basename "${i}" _clean_1.fastq.gz)
        echo "name is $name"
   
        in1=${dname}/${name}_clean_1.fastq.gz
        in2=${dname}/${name}_clean_2.fastq.gz
        bam=$aligned_dir/${name}_aligned.bam
        sorted_bam=$aligned_dir/${name}_aligned_sorted.bam
        rmd_bam=$aligned_dir/${name}_rmd.bam
		  
	bwa-mem2 mem -t 30 -R "@RG\tID:${name}\tSM:${name}\tPL:illumina\tLB:lib1\tPU:unit1" $genome $in1 $in2 | samtools view -@ 30 -bSu - > $bam
	samtools sort -@ 30 -o $sorted_bam $bam
	samtools index -@ 30 $sorted_bam    
	gatk --java-options "-Xmx140G" MarkDuplicates I=$sorted_bam O=$rmd_bam REMOVE_DUPLICATES=true M=$qc_dir/${name}.duplicates.txt
	samtools index -@ 30 $rmd_bam
	rm -f "$bam"
	rm -f "$sorted_bam"
	rm -f "${sorted_bam}".bai
done 
## get the total and the mapped only number of reads of a BAM file 
## legend
# samtools view -c = count reads, -F flag (filter out) reads: https://broadinstitute.github.io/picard/explain-flags.html
# -F 260 = filter out unmapped and secondary aligned reads
# -F 256 = filter out only secondary aligned reads
for i in "$aligned_dir"/*_rmd.bam;
do
	base_name=$(basename "$i" _rmd.bam)
	count_total=$(samtools view -c "$i")
	count_mapped=$(samtools view -c -F 260 "$i")
	echo "$base_name, $count_total, $count_mapped" >> "$prop_aligned"
done
