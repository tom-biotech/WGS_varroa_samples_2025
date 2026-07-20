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
--filter-sample-min-count 3 \
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

    hist[chrom, depth]++
    if (depth > maxdepth[chrom]) maxdepth[chrom] = depth
}
END {
    for (chrom in count) {
        n = count[chrom]
        variance = (n > 1) ? M2[chrom] / (n - 1) : 0
        sd = sqrt(variance)

        target1 = int((n + 1) / 2)
        target2 = int(n / 2 + 1)
        cum = 0
        med1 = -1
        med2 = -1
        for (d = 0; d <= maxdepth[chrom]; d++) {
            if (!((chrom, d) in hist)) continue
            cum += hist[chrom, d]
            if (med1 == -1 && cum >= target1) med1 = d
            if (med2 == -1 && cum >= target2) med2 = d
            if (med1 != -1 && med2 != -1) break
        }
        median = (n % 2 == 1) ? med1 : (med1 + med2) / 2

        printf "%s\t%d\t%.4f\t%.4f\t%.4f\n", chrom, n, mean[chrom], sd, median
    }
}
----------
# Depth stats per sample and chromosom
mkdir -p stats_per_sample
ls /home/tomsch/WGS_36/aligned_new/sync_files/*.sync | parallel -j 10 \
    'sample=$(basename {} .sync); awk -f depth_stats_per_sample.awk {} > stats_per_sample/${sample}.stats'

echo -e "sample\tchrom\tn_positions\tmean_depth\tsd_depth\tmedian_depth" > depth_stats_final.tsv
for f in stats_per_sample/*.stats; do
    sample=$(basename "$f" .stats)
    awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}' "$f"
done >> stats_per_sample/depth_stats_final.tsv

awk 'NR==FNR{order[$1]=NR; next} FNR>1{print $0, order[$2]}' /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt depth_stats_final.tsv \
    | sort -k7,7n \
    | cut -d' ' -f1-6 > depth_stats_final_sorted.tsv

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

    hist[depth]++
    if (depth > maxdepth) maxdepth = depth
}
END {
    variance = (count > 1) ? M2 / (count - 1) : 0
    sd = sqrt(variance)
    breadth = (covered / count) * 100

    # Median über kumulative Häufigkeiten im Histogramm bestimmen
    target1 = int((count + 1) / 2)
    target2 = int(count / 2 + 1)
    cum = 0
    med1 = -1
    med2 = -1
    for (d = 0; d <= maxdepth; d++) {
        if (!(d in hist)) continue
        cum += hist[d]
        if (med1 == -1 && cum >= target1) med1 = d
        if (med2 == -1 && cum >= target2) med2 = d
        if (med1 != -1 && med2 != -1) break
    }
    median = (count % 2 == 1) ? med1 : (med1 + med2) / 2

    printf "%d\t%.4f\t%.4f\t%d\t%.4f\t%.4f\n", count, mean, sd, covered, breadth, median
}
-------------

ls /home/tomsch/WGS_36/aligned_new/sync_files/*.sync | parallel -j 20 \
    'sample=$(basename {} .sync); echo -e "$sample\t$(awk -f mean_depth_per_sample.awk {})"' \
    > stats_per_sample/mean_depth_stats_per_sample.tsv
	
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' mean_depth_stats_per_sample.tsv 

###### --> Mean Depth over all samples and Chr = 175.604x
###### --> Median Depth = 184.528x

/home/tomsch/grenedalf/bin/grenedalf fst \
--method unbiased-hudson \
--window-type genome \
--write-pi-tables \
--sync-path /home/tomsch/WGS_36/aligned_new/sync_files \
--reference-genome-fasta /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
--filter-sample-min-count 3 \
--filter-sample-min-read-depth 20 \
--filter-sample-max-read-depth 370 \
--window-average-policy valid-loci \
--filter-total-snp-min-frequency 0.01 \
--filter-total-only-biallelic-snps \
--pool-sizes 60 \
--file-prefix all_samples_fst_calculation_ \
--out-dir /home/tomsch/WGS_36/aligned_new/fst_files/all_samples \
--compress \
--log-file /home/tomsch/WGS_36/aligned_new/fst_files/all_samples/all_samples_fst.log \
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
