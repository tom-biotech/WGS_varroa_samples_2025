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
export genome="$genome"
export bam_qc_dir="$bam_qc_dir"

parallel -j 20 '
name=$(basename {} _rmd_sub_30.bam)
alfred qc \
  -r "$genome" \
  -j "$bam_qc_dir"/${name}_qc.json.gz \
  -o "$bam_qc_dir"/${name}_qc.tsv.gz \
  {}
' ::: "$aligned_dir"*_rmd_sub_30.bam

############################
# convert bam to pileup file
############################

mkdir -p "/home/tomsch/WGS_36/sub_aligned/mpileup_files"
mkdir -p "/home/tomsch/WGS_36/sub_aligned/sync_files"
mpileup_dir="/home/tomsch/WGS_36/sub_aligned/mpileup_files"
sync_dir="/home/tomsch/WGS_36/sub_aligned/sync_files"

# from samtools mpileup to sync file
for i in "$aligned_dir"/B5047-SCH-{36..60}_rmd.bam; do name=$(basename ${i} _rmd_sub_30.bam);
cat /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt | parallel -j 20 \
"samtools mpileup -B -f $genome \
-q 40 -Q 20 -aa -r {} ${i} > "$mpileup_dir"/${name}_{}.mpileup"
while read c; do cat "$mpileup_dir"/${name}_${c}.mpileup; done < /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt > "$mpileup_dir"/${name}.mpileup
rm "$mpileup_dir"/${name}_N*
java -ea -Xmx10g -jar \
/home/tomsch/miniconda3/envs/WGS_36/share/popoolation2-1.201-0/mpileup2sync.jar --input "$mpileup_dir"/${name}.mpileup --output "$sync_dir"/${name}.sync --fastq-type sanger --min-qual 20 --threads 20;
done

#####
# Depth in sync files

depth_stats_per_sample.awk
----------
BEGIN{OFS="\t"}
{
    chrom = $1
    n = split($4, counts, ":")
    depth = 0
    for (i = 1; i <= n; i++) depth += counts[i]

    count[chrom]++
    delta = depth - mean[chrom]
    mean[chrom] += delta / count[chrom]
    delta2 = depth - mean[chrom]
    M2[chrom] += delta * delta2
}
END {
    for (chrom in count) {
        n = count[chrom]
        variance = (n > 1) ? M2[chrom] / (n - 1) : 0
        sd = sqrt(variance)
        printf "%s\t%d\t%.4f\t%.4f\n", chrom, n, mean[chrom], sd
    }
}
----------

mkdir -p stats_per_sample

ls /home/tomsch/WGS_36/sync/*.sync | parallel -j 8 \
    'sample=$(basename {} .sync); awk -f depth_stats_per_sample.awk {} > stats_per_sample/${sample}.stats'

echo -e "sample\tchrom\tn_positions\tmean_depth\tsd_depth" > depth_stats_final.tsv

for f in stats_per_sample/*.stats; do
    sample=$(basename "$f" .stats)
    awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}' "$f"
done >> depth_stats_final.tsv

awk 'NR==FNR{order[$1]=NR; next} FNR>1{print $0, order[$2]}' chromosomes.txt depth_stats_final.tsv \
    | sort -k6,6n \
    | cut -d' ' -f1-5 > depth_stats_final_sorted.tsv
#####
# Fst calculation

/home/tomsch/grenedalf/bin/grenedalf fst 
--method unbiased-hudson 
--window-type genome 
--write-pi-tables 
--sync-path /home/tomsch/WGS_36/sub_aligned/sync_files 
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna 
--filter-sample-min-count 2 
--filter-sample-min-read-depth 86 
--filter-sample-max-read-depth 344 
--window-average-policy valid-loci 
--filter-total-snp-min-frequency 0.01 
--pool-sizes 60 
--file-prefix all_samples_fst_calculation_ 
--out-dir /home/tomsch/WGS_36/sub_aligned/fst_files/all_samples 
--compress 
--log-file /home/tomsch/WGS_36/sub_aligned/fst_files/all_samples/all_samples_fst.log 
--threads 20
