#!/bin/bash

############# VARROARESISTENZ WGS 2026 ################ 
#### using bbmap 39.85
#### using bwa-mem2 2.3
#### using samtools 1.23.1
#### using pandepth 2.26
#### using alfred 0.5.3

# The goal is to minimize the required sequencing data, to get the same allele frequencies like in the full dataset.
# Starting from the clean fastq data

mkdir -p "/home/tomsch/WGS_36/Clean/subsample"
sub_dir="/home/tomsch/WGS_36/Clean/subsample"
fastq_clean="/home/tomsch/WGS_36/Clean"
genome="/home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
mkdir -p "/home/tomsch/WGS_36/QC_sub"
qc_dir="/home/tomsch/WGS_36/QC_sub"
prop_aligned="$qc_dir/counts_summary_all_sub_samples.txt"
mkdir -p "/home/tomsch/WGS_36/sub_aligned"
aligned_dir="/home/tomsch/WGS_36/sub_aligned"

# Subsample to a depth of 30x (genome: 225250884, depth 30x -> 225250884 * 30 ‎ = 6.757.526.520  -> 6757526520/150 (read length) ‎ = 45.050.176,8 -> 45050176,8 /2 (paired)‎ = 22.525.088,4)
'''
subsample_one() {
    local i="$1"
    local fastq_clean="$2"
    local sub_dir="$3"
    local name
    name=$(basename "${i}" _clean_1.fastq.gz)
    echo "name is $name"
    local in1="$fastq_clean/${name}_clean_1.fastq.gz"
    local in2="$fastq_clean/${name}_clean_2.fastq.gz"
    seqtk sample -s100 "$in1" 22525088 | gzip -c > "$sub_dir/${name}_1_sub_30.fastq.gz"
    seqtk sample -s100 "$in2" 22525088 | gzip -c > "$sub_dir/${name}_2_sub_30.fastq.gz"
}
export -f subsample_one
parallel -j 10 subsample_one {} "$fastq_clean" "$sub_dir" ::: "$fastq_clean"/*_clean_1.fastq.gz
'''
for i in "$sub_dir"/*_1_sub_30.fastq.gz
        do dname=$(dirname "${i}"); name=$(basename "${i}" _1_sub_30.fastq.gz)
        echo "name is $name"
   
        in1=${dname}/${name}_1_sub_30.fastq.gz
        in2=${dname}/${name}_2_sub_30.fastq.gz
        bam=$aligned_dir/${name}_aligned_sub_30.bam
        sorted_bam=$aligned_dir/${name}_aligned_sub_30_sorted.bam
        rmd_bam=$aligned_dir/${name}_rmd_sub_30.bam
		  
	bwa-mem2 mem -t 20 -R "@RG\tID:${name}\tSM:${name}\tPL:illumina\tLB:lib1\tPU:unit1" $genome $in1 $in2 | samtools view -@ 20 -bSu - > $bam
	samtools sort -@ 20 -o $sorted_bam $bam
	samtools index -@ 20 $sorted_bam    
	gatk --java-options "-Xmx140G" MarkDuplicates I=$sorted_bam O=$rmd_bam REMOVE_DUPLICATES=true M=$qc_dir/${name}.duplicates.txt
	samtools index -@ 20 $rmd_bam
	rm -f "$bam"
	rm -f "$sorted_bam"
	rm -f "${sorted_bam}".bai
done 
## get the total and the mapped only number of reads of a BAM file 
## legend
# samtools view -c = count reads, -F flag (filter out) reads: https://broadinstitute.github.io/picard/explain-flags.html
# -F 260 = filter out unmapped and secondary aligned reads
# -F 256 = filter out only secondary aligned reads
for i in "$aligned_dir"/*_rmd_sub_30.bam;
do
	base_name=$(basename "$i" _rmd_sub_30.bam)
	count_total=$(samtools view -c "$i")
	count_mapped=$(samtools view -c -F 260 "$i")
	echo "$base_name, $count_total, $count_mapped" >> "$prop_aligned"
done

## statistics for read depth and general bam QC
## Pandepth
for i in "$aligned_dir"/*.bam; 
do 
name=$(basename ${i} _rmd_sub_30.bam);
pandepth -i $i -o "$aligned_dir"/${name}_depth -t 20; 
done

mkdir "$aligned_dir"/bam_qc
bam_qc_dir="/home/tomsch/WGS_36/sub_aligned/bam_qc"

## alfred
parallel -j 20 '
name=$(basename {} _rmd_sub_30.bam)
alfred qc \
  -r "$genome" \
  -j "$bam_qc_dir"/${name}_qc.json.gz \
  -o "$bam_qc_dir"/${name}_qc.tsv.gz \
  {}
' ::: "$aligned_dir"*_rmd_sub_30.bam ::: "$genome" ::: "$bam_qc_dir"
