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
        --sam-min-map-qual 40 \
        --sam-min-base-qual 20 \
        --compress \
        --log-file /home/tomsch/WGS_36/aligned/sync_files/${name}_sync.log \
        --threads 5 \
        --file-prefix ${name}_ \
        --out-dir /home/tomsch/WGS_36/aligned/sync_files;
        done
'''
#############################################################
# from samtools mpileup to sync file
cat /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt | parallel -j 20 \
"samtools mpileup -B -f /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
-q 40 -Q 20 -aa -r {} -b bam.list > mpileup_files/WGS_36.{}.mpileup"

# danach in der richtigen Reihenfolge zusammenführen (Reihenfolge aus chromosomes.txt!)
while read chr; do
    cat mpileup_files/WGS_36.${chr}.mpileup
done < /home/tomsch/WGS_36/Amel_HAv3.1/ncbi_dataset/data/GCF_003254395.2/chromosomes.txt > mpileup_files/WGS_36.mpileup

/home/tomsch/miniconda3/envs/WGS_36/share/popoolation2-1.201-0/mpileup2sync.jar --input /home/tomsch/WGS_36/aligned_new/mpileup_files/WGS_36.mpileup --output /home/tomsch/WGS_36/aligned_new/sync_files/WGS_36.sync --fastq-type sanger --min-qual 20 --threads 20

''' for every bam file make their own mpileup
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
'''
#############################################################
# find and filter indels with popoolation 2

/home/tomsch/miniconda3/envs/WGS_36/bin/identify-indel-regions.pl --input WGS_36.mpileup --output WGS_36_indel_regions.gtf --min-count 2 --indel-window 5

/home/tomsch/miniconda3/envs/WGS_36/bin/filter-sync-by-gtf.pl --input $i --output ${name}_rm_indel.sync --gtf ${name}_indel_regions.gtf


#############################################################


# create depth file from sync (same as samtools depth)
# sync_to_depth.sh

#!/bin/bash
SYNC_DIR="/home/tomsch/WGS_36/aligned_new/sync_files"
OUTDIR="/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync"
mkdir -p "$OUTDIR"

export OUTDIR

process_one() {
    local sync_file="$1"
    local sample
    sample=$(basename "$sync_file" .sync)
    local outfile="${OUTDIR}/${sample}_depth.txt"

    {
      echo -e "#CHROM\tPOS\t${sample}"
      gawk 'BEGIN{FS="\t"} {
          split($4, cnt, ":")
          depth = cnt[1]+cnt[2]+cnt[3]+cnt[4]+cnt[5]+cnt[6]
          print $1"\t"$2"\t"depth
      }' "$sync_file"
    } > "$outfile"
}
export -f process_one

parallel -j 5 process_one ::: "$SYNC_DIR"/*.sync

# delete Mitochondria from depth file (because it have a way higher depth than the normal chromosomes)

for i in *_depth.txt; do name=$(basename "$i" _depth.txt); awk '$1 != "NC_001566.1"' $i | gzip > ${name}_depth_without_mito.txt.gz; echo "Probe ${name} bearbeitet..."; done

# OR take just the mito genome depth -> for analysis of mitochondrium

for i in *_depth.txt.gz; do name=$(basename "$i" _depth.txt.gz); zcat $i | awk '$1 == "NC_001566.1"' | gzip > ${name}_depth_mito.txt.gz; echo "Probe ${name} bearbeitet..."; done

###########################################################################################################################
# from there on with an R-Script calculate elbow depth of genome WITHOUT MITO
###########################################################################################################################

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


###########################################################################################################################
# from there on with an R-Script calculate elbow depth of MITO GENOME
###########################################################################################################################

library(data.table)
library(ggplot2)
library(scales)

setDTthreads(20)

# --- Konfiguration ---
depth_dir   <- "/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync"
out_dir     <- "/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync/density"
dir.create(out_dir, showWarnings = FALSE)

depth_files <- list.files(depth_dir, pattern = "_depth_mito\\.txt\\.gz$", full.names = TRUE)
sample_names <- sub("_depth_mito\\.txt\\.gz$", "", basename(depth_files))

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

fwrite(stats_dt, file.path(out_dir, "depth_summary_stats_mito.csv"))

print(stats_dt)




###########################################################################################################################
# create density plot of all samples and the genome WITHOUT MITO
###########################################################################################################################

library(data.table)
library(ggplot2)
library(scales)
library(plotly)
library(htmlwidgets)

base_dir   <- "/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync"
out_dir    <- file.path(base_dir, "density")
dir.create(out_dir, showWarnings = FALSE)

sample_numbers <- 25:60
density_list   <- vector("list", length(sample_numbers))

for (i in seq_along(sample_numbers)) {

  n <- sample_numbers[i]
  sample_id  <- paste0("B5047-SCH-", n)
  depth_file <- file.path(base_dir, paste0(sample_id, "_depth_without_mito.txt.gz"))

  if (!file.exists(depth_file)) {
    message("Datei fehlt, überspringe: ", depth_file)
    next
  }

  message("Verarbeite ", sample_id, " ...")

  depth_dt <- fread(cmd = paste0("zcat ", depth_file, " | grep -v '^#' | awk '{print $3}'"),
                   header = FALSE)$V1          # <- extract the column, not the whole data.table

  depth_dt <- as.numeric(depth_dt)
  depth_dt <- depth_dt[!is.na(depth_dt)]
  depth_dt <- depth_dt[depth_dt > 0]

  upper_cut <- quantile(depth_dt, 0.999)
  
  # Dichte berechnen (kompakte Repräsentation, n=4096 Stützstellen)
  dens <- density(depth_dt, n = 4096, from = 0, to = upper_cut, na.rm = TRUE)

  density_list[[i]] <- data.table(
    sample = sample_id,
    x = dens$x,
    y = dens$y
  )

  # Rohdaten sofort freigeben, bevor die nächste Probe geladen wird
  rm(depth_dt, dens, upper_cut)
  gc()
}

# Alle Dichtekurven zusammenführen (leichtgewichtig: 36 x 4096 Zeilen statt Milliarden)
all_dens <- rbindlist(density_list, use.names = TRUE)

p <- ggplot(
  all_dens,
  aes(
    x = x,
    y = y,
    color = sample,
    group = sample,
    text = paste0(
      "Probe: ", sample,
      "<br>Depth: ", round(x, 1),
      "<br>Density: ", signif(y, 4)
    )
  )
) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    #trans = pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 20, 100, 200, 300, 400, 500, 600),
    labels = c(0, 20, 100, 200, 300, 400, 500, 600)
  ) +
  coord_cartesian(xlim = c(0, 600)) +
  labs(
    title = "Coverage Depth Distribution - all samples",
    x = "Depth",
    y = "Density",
    color = "Probe"
  ) +
  theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm"))

p_html <- ggplotly(p, tooltip = "text")

saveWidget(
  p_html,
  file = file.path(out_dir, "all_samples_depth_density_log.html"),
  selfcontained = TRUE
)

ggsave(
  filename = file.path(out_dir, "all_samples_depth_density_log.png"),
  plot = p, width = 10, height = 6, dpi = 300
)

message("Fertig – kombinierter Plot liegt in: ", out_dir)


###########################################################################################################################
# create density plot of all samples and the MITO GENOME
###########################################################################################################################

library(data.table)
library(ggplot2)
library(scales)
library(plotly)
library(htmlwidgets)

base_dir   <- "/home/tomsch/WGS_36/aligned_new/sync_files/depth_from_sync"
out_dir    <- file.path(base_dir, "density")
dir.create(out_dir, showWarnings = FALSE)

sample_numbers <- 25:60
density_list   <- vector("list", length(sample_numbers))

for (i in seq_along(sample_numbers)) {

  n <- sample_numbers[i]
  sample_id  <- paste0("B5047-SCH-", n)
  depth_file <- file.path(base_dir, paste0(sample_id, "_depth_mito.txt.gz"))

  if (!file.exists(depth_file)) {
    message("Datei fehlt, überspringe: ", depth_file)
    next
  }

  message("Verarbeite ", sample_id, " ...")

  depth_dt <- fread(cmd = paste0("zcat ", depth_file, " | grep -v '^#' | awk '{print $3}'"),
                   header = FALSE)$V1          # <- extract the column, not the whole data.table

  depth_dt <- as.numeric(depth_dt)
  depth_dt <- depth_dt[!is.na(depth_dt)]
  depth_dt <- depth_dt[depth_dt > 0]

  upper_cut <- quantile(depth_dt, 0.999)
  
  # Dichte berechnen (kompakte Repräsentation, n=4096 Stützstellen)
  dens <- density(depth_dt, n = 4096, from = 0, to = upper_cut, na.rm = TRUE)

  density_list[[i]] <- data.table(
    sample = sample_id,
    x = dens$x,
    y = dens$y
  )

  # Rohdaten sofort freigeben, bevor die nächste Probe geladen wird
  rm(depth_dt, dens, upper_cut)
  gc()
}

# Alle Dichtekurven zusammenführen (leichtgewichtig: 36 x 4096 Zeilen statt Milliarden)
all_dens <- rbindlist(density_list, use.names = TRUE)

p <- ggplot(
  all_dens,
  aes(
    x = x,
    y = y,
    color = sample,
    group = sample,
    text = paste0(
      "Probe: ", sample,
      "<br>Depth: ", round(x, 1),
      "<br>Density: ", signif(y, 4)
    )
  )
) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    #trans = pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 20, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000),
    labels = c(0, 20, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000)
  ) +
  coord_cartesian(xlim = c(0, 9000)) +
  labs(
    title = "Coverage Depth Distribution - all samples mitochondria",
    x = "Depth",
    y = "Density",
    color = "Probe"
  ) +
  theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm"))

p_html <- ggplotly(p, tooltip = "text")

saveWidget(
  p_html,
  file = file.path(out_dir, "all_samples_MITO_depth_density_log.html"),
  selfcontained = TRUE
)

ggsave(
  filename = file.path(out_dir, "all_samples_MITO_depth_density_log.png"),
  plot = p, width = 10, height = 6, dpi = 300
)

message("Fertig – kombinierter Plot liegt in: ", out_dir)

###########################################################################################################################
# mask_sync_by_elbow.sh
# mask sync files on the basis of min depth of 20 and max depth of elbow value for every sample
###########################################################################################################################

#!/bin/bash
#
# mask_sync_by_elbow.sh
#
# Maskiert Positionen in den Pro-Sample-sync-Dateien, deren Depth
# NICHT im Bereich [min_depth, elbow_depth(sample)] liegt.
# Maskierte Positionen bekommen "0:0:0:0:0:0" im Zählfeld.
#
# Usage:
#   ./mask_sync_by_elbow.sh <depth_summary.csv> <sync_dir> <out_dir> [min_depth]
#
# depth_summary.csv  = deine R-Ausgabe mit Spalten "sample" und "elbow_depth"
#                       (Spaltenreihenfolge egal, muss aber diese Header haben)
# sync_dir            = Ordner mit den originalen Pro-Sample sync-Dateien
# out_dir             = Zielordner (Standard-Empfehlung: masked_sync/)
# min_depth           = optional, Standard = 20
#
set -euo pipefail

summary_csv="$1"
sync_dir="$2"
out_dir="$3"
min_depth="${4:-20}"

mkdir -p "$out_dir"

# --- 1. sample -> elbow_depth Lookup-Tabelle bauen (tab-separiert, ohne Header) ---
mapping_file=$(mktemp)

gawk -F',' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/"/, "", $i)
        if ($i == "sample")      s_col = i
        if ($i == "elbow_depth") e_col = i
    }
    if (!s_col || !e_col) {
        print "ERROR: Spalten \"sample\" und/oder \"elbow_depth\" nicht gefunden" > "/dev/stderr"
        exit 1
    }
    next
}
{
    gsub(/"/, "")
    print $s_col "\t" $e_col
}' "$summary_csv" > "$mapping_file"

n_samples=$(wc -l < "$mapping_file")
echo "Gefunden: $n_samples Samples in $summary_csv"
echo "min_depth = $min_depth, out_dir = $out_dir"
echo

# --- 2. Maskier-Funktion, einmal pro Sample aufgerufen ---
mask_sample() {
    local sample="$1"
    local elbow="$2"
    local min="$3"
    local sdir="$4"
    local odir="$5"

    local infile
    infile=$(find "$sdir" -maxdepth 1 -type f -iname "${sample}*sync*" | head -n 1)

    if [[ -z "$infile" ]]; then
        echo "WARNUNG: keine sync-Datei fuer $sample in $sdir gefunden" >&2
        return 1
    fi

    local base
    base=$(basename "$infile")
    local outfile="${odir}/${base%.sync}_masked.sync"

    gawk -v min="$min" -v max="$elbow" '
    BEGIN { FS = "\t"; OFS = "\t" }
    {
        n = split($4, a, ":")
        depth = 0
        for (i = 1; i <= n; i++) depth += a[i]
        if (depth < min || depth > max) {
            $4 = "0:0:0:0:0:0"
        }
        print
    }' "$infile" > "$outfile"

    echo "OK: $sample -> $(basename "$outfile")  (min=$min, max=$elbow)"
}
export -f mask_sample

# --- 3. Parallel ueber alle Samples ---
parallel --colsep '\t' --jobs 8 \
    mask_sample {1} {2} "$min_depth" "$sync_dir" "$out_dir" \
    :::: "$mapping_file"

rm -f "$mapping_file"

echo
echo "Fertig. Maskierte Dateien liegen in: $out_dir"

