#!/bin/bash

############# VARROARESISTENZ WGS 2026 ################ 
#### using bwa-mem2 2.3
#### using samtools 1.23.1
#### using pandepth 2.26
#### using alfred 0.5.3


bam_dir="/home/tomsch/WGS_36/aligned"
mkdir -p "/home/tomsch/WGS_36/Clean/subsample/"
fastq_dir="/home/tomsch/WGS_36/Clean/subsample/"
genome="/home/tomsch/WGS_36/Vdes_3.0/Varroa_destructor_mito.fasta"
mkdir -p "/home/tomsch/WGS_36/unmapped_bam"
unmap_dir="/home/tomsch/WGS_36/unmapped_bam"

# using the unmapped reads from mapping against the Amel_HAv3.1 to mapp against mito genome of Varroa destructor

# extract unmapped reads with samtools

for i in $bam_dir/*_rmd.bam;
do name=$(basename $i _rmd.bam);

samtools view -f 4 -@ 20 $i | samtools fastq -c 3 -1 $fastq_dir/${name}_1.fastq.gz -2 $fastq_dir/${name}_2.fastq.gz -@ 10 -s $fastq_dir/${name}_singletons.fastq.gz 

done

## index genome
samtools faidx $genome
bwa-mem2 index $genome

for i in "$fastq_dir"/*_1_sub_30.fastq.gz
        do dname=$(dirname "${i}"); name=$(basename "${i}" _1_sub_30.fastq.gz)
        echo "name is $name"
   
        in1=${dname}/${name}_1_sub_30.fastq.gz
        in2=${dname}/${name}_2_sub_30.fastq.gz
        bam=$unmap_dir/${name}_aligned.bam
        sorted_bam=$unmap_dir/${name}_aligned_sorted.bam
        rmd_bam=$unmap_dir/${name}_varroa_mapped.bam
		  
	bwa-mem2 mem -t 30 -R "@RG\tID:${name}\tSM:${name}\tPL:illumina\tLB:lib1\tPU:unit1" $genome $in1 $in2 | samtools view -@ 30 -bSu - > $bam
	samtools sort -@ 30 -o $sorted_bam $bam
	samtools index -@ 30 $sorted_bam    
	gatk --java-options "-Xmx140G" MarkDuplicates I=$sorted_bam O=$rmd_bam REMOVE_DUPLICATES=true M=$unmap_dir/${name}.duplicates.txt
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
for i in "$unmap_dir"/*_varroa_mapped.bam;
do
	base_name=$(basename "$i" _varroa_mapped.bam)
	count_total=$(samtools view -c "$i")
	count_mapped=$(samtools view -c -F 260 "$i")
	echo "$base_name, $count_total, $count_mapped" >> "$unmap_dir/all_count_stats.txt"
done

## statistics for read depth and general bam QC
## Pandepth
for i in "$unmap_dir"/*.bam; 
do 
name=$(basename ${i} _rmd.bam);
pandepth -i $i -o "$unmap_dir"/${name}_depth -t 20; 
done

mkdir "$unmap_dir"/bam_qc
bam_qc_dir="/home/tomsch/WGS_36/unmapped_bam/bam_qc"

## alfred
parallel -j 20 '
name=$(basename {} _varroa_mapped.bam)
alfred qc \
  -r "$genome" \
  -j "$bam_qc_dir"/${name}_qc.json.gz \
  -o "$bam_qc_dir"/${name}_qc.tsv.gz \
  {}
' ::: "$unmap_dir"*_varroa_mapped.bam
