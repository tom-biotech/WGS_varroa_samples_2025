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

'''
echo -e "Sample\tCoverage\tMeanDepth" > summary.tsv

for f in B5047-SCH-*.chr.stat.gz
do
    sample=$(basename "$f" .chr.stat.gz)
    zgrep '^##' "$f" | tail -n1 | \
    awk -v sample="$sample" '
    {
        for(i=1;i<=NF;i++){
            if($i ~ /^Coverage/) cov=$(i+1)
            if($i ~ /^MeanDepth/) md=$(i+1)
        }
        print sample "\t" cov "\t" md
    }' >> pan_depth_summary.tsv
done
'''

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
# convert bam to pileup and sync file
############################

mkdir -p "/home/tomsch/WGS_36/sub_aligned/mpileup_files"
mkdir -p "/home/tomsch/WGS_36/sub_aligned/sync_files"
mpileup_dir="/home/tomsch/WGS_36/sub_aligned/mpileup_files"
sync_dir="/home/tomsch/WGS_36/sub_aligned/sync_files"

# from samtools mpileup to sync file
for i in "$aligned_dir"/B5047-SCH-{25..60}_rmd_sub_30.bam; do name=$(basename ${i} _rmd_sub_30.bam);
cat /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt | parallel -j 20 \
"samtools mpileup -B -f $genome \
-q 40 -Q 20 -aa -r {} ${i} > "$mpileup_dir"/${name}_{}.mpileup"
while read c; do cat "$mpileup_dir"/${name}_${c}.mpileup; done < /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt > "$mpileup_dir"/${name}.mpileup
rm "$mpileup_dir"/${name}_N*
java -ea -Xmx10g -jar \
/home/tomsch/miniconda3/envs/WGS_36/share/popoolation2-1.201-0/mpileup2sync.jar --input "$mpileup_dir"/${name}.mpileup --output "$sync_dir"/${name}.sync --fastq-type sanger --min-qual 20 --threads 20;
done
'''
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
# Depth stats per sample and chromosom
mkdir -p stats_per_sample

ls /home/tomsch/WGS_36/sub_aligned/sync_files/*.sync | parallel -j 10 \
    'sample=$(basename {} .sync); awk -f depth_stats_per_sample.awk {} > stats_per_sample/${sample}.stats'

echo -e "sample\tchrom\tn_positions\tmean_depth\tsd_depth" > depth_stats_final.tsv

for f in stats_per_sample/*.stats; do
    sample=$(basename "$f" .stats)
    awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}' "$f"
done >> stats_per_sample/depth_stats_final.tsv

awk 'NR==FNR{order[$1]=NR; next} FNR>1{print $0, order[$2]}' /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt depth_stats_final.tsv \
    | sort -k6,6n \
    | cut -d' ' -f1-5 > depth_stats_final_sorted.tsv

# Mean depth per sample
mean_depth_per_sample.awk:
-------------
{
    n = split($4, counts, ":")
    depth = 0
    for (i = 1; i <= n; i++) depth += counts[i]

    count++
    delta = depth - mean
    mean += delta / count
    delta2 = depth - mean
    M2 += delta * delta2

    if (depth > 0) covered++
}
END {
    variance = (count > 1) ? M2 / (count - 1) : 0
    sd = sqrt(variance)
    breadth = (covered / count) * 100
    printf "%d\t%.4f\t%.4f\t%d\t%.4f\n", count, mean, sd, covered, breadth
}
-------------

ls /home/tomsch/WGS_36/sub_aligned/sync_files/*.sync | parallel -j 10 \
    'sample=$(basename {} .sync); echo -e "$sample\t$(awk -f mean_depth_per_sample.awk {})"' \
    > stats_per_sample/mean_depth_stats_per_sample.tsv
	
awk '{ total += $3 } END { print total/NR }' mean_depth_stats_per_sample.tsv
'''
#####
# Fst calculation for all subsampled sync files (one per sample)
# Mean Depth: 24.1

/home/tomsch/grenedalf/bin/grenedalf fst \
--method unbiased-hudson \
--window-type genome \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/sub_aligned/sync_files \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 12 \
--filter-sample-max-read-depth 48 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--pool-sizes 60 \
--file-prefix all_sub_30_samples_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/sub_aligned/fst_files/all_samples \
--compress \
--log-file /home/tomsch/WGS_36/sub_aligned/fst_files/all_samples/all_samples_fst.log \
--threads 20

# calculating significance of allele frequency differences (with popoolation2, Fisher's Exact Test)

######################################################################################################
# for the different subsamples of B5047-SCH-36 (seqtk -s100 - -s109)
######################################################################################################
# calculate depth elbow to mask sync files to min 5 and max Elbow depth
# delete Mitochondria from depth file (because it have a way higher depth than the normal chromosomes)

for i in *_depth.txt; do name=$(basename "$i" _depth.txt); awk '$1 != "NC_001566.1"' $i | gzip > ${name}_depth_without_mito.txt.gz; echo "Probe ${name} bearbeitet..."; done

# from there on with an R-Script calculate elbow depth

library(data.table)
library(ggplot2)
library(scales)

setDTthreads(20)

# --- Konfiguration ---
depth_dir   <- "/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync"
out_dir     <- "/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync/density"
dir.create(out_dir, showWarnings = FALSE)

depth_files <- list.files(depth_dir, pattern = "_depth_without_mito\\.txt\\.gz$", full.names = TRUE)
sample_names <- sub("_depth_without_mito\\.txt\\.gz$", "", basename(depth_files))

cat("Gefundene Samples:", length(depth_files), "\n")

# --- Funktion: Ellbogen via Kneedle auf einer density()-Schätzung ---
find_elbow <- function(d) {
  peak_idx <- which.max(d$y)
  x_desc <- d$x[peak_idx:length(d$x)]
  y_desc <- d$y[peak_idx:length(d$x)]

  x_norm <- (x_desc - min(x_desc)) / (max(x_desc) - min(x_desc))
  y_norm <- (y_desc - min(y_desc)) / (max(y_desc) - min(y_desc))

  line_y <- y_norm[1] + (y_norm[length(y_norm)] - y_norm[1]) *
            (x_norm - x_norm[1]) / (x_norm[length(x_norm)] - x_norm[1])
  dist_to_line <- line_y - y_norm

  x_desc[which.max(dist_to_line)]
}

# --- Container für Ergebnisse ---
density_list <- vector("list", length(depth_files))
stats_list   <- vector("list", length(depth_files))

# --- Hauptschleife: ein Sample nach dem anderen ---
for (i in seq_along(depth_files)) {

  f <- depth_files[i]
  s <- sample_names[i]
  cat(sprintf("[%d/%d] Verarbeite %s ...\n", i, length(depth_files), s))

  # Nur die Depth-Spalte extrahieren (awk $3) -> deutlich weniger RAM als 3 Spalten
  depth_vec <- fread(
    cmd = paste0("zcat ", f, " | grep -v '^#' | awk '{print $3}'"),
    header = FALSE
  )[[1]]
  depth_vec <- as.numeric(depth_vec)
  depth_vec <- depth_vec[!is.na(depth_vec)]
  depth_vec <- depth_vec[depth_vec > 0]

  # Kennzahlen
  mean_d   <- mean(depth_vec)
  median_d <- median(depth_vec)
  q95      <- quantile(depth_vec, 0.95, names = FALSE)
  q99      <- quantile(depth_vec, 0.99, names = FALSE)
  upper_cut <- quantile(depth_vec, 0.999, names = FALSE)

  # Density nur bis 99.9%-Quantil schätzen (Tail sonst verzerrend)
  d <- density(depth_vec, from = 0, to = upper_cut, n = 4096)

  modus_d <- d$x[which.max(d$y)]
  elbow_d <- find_elbow(d)

  stats_list[[i]] <- data.table(
    sample      = s,
    elbow_depth = elbow_d,
    q95         = q95,
    q99         = q99,
    mean        = mean_d,
    median      = median_d,
    modus       = modus_d
  )

  density_list[[i]] <- data.table(sample = s, x = d$x, y = d$y)

  # Rohdaten sofort verwerfen
  rm(depth_vec, d)
  gc(verbose = FALSE)
}

# --- Zusammenführen ---
stats_dt   <- rbindlist(stats_list)
density_dt <- rbindlist(density_list)

fwrite(stats_dt, file.path(out_dir, "depth_summary_stats.csv"))

print(stats_dt)




# merge sync files in on sync file


perl <popoolation2-path>/fisher-test.pl --input p1_p2.sync --output p1_p2.fet --min-count 6 --min-coverage 10 --max-coverage 200 --suppress-noninformative
