#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=256G
#SBATCH --partition=HTC
#SBATCH --mail-user=jah550@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# -----------------------------
# Environment setup
# -----------------------------
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate diffbind
module load bedtools

OUTPUT_DIR="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/HOMER"

# =====================================================
# STEP 1: Define file paths to bed and fam files 
# =====================================================

Rscript - << 'EOF'

# Load all required libraries
library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggrepel)
library(ggplot2)

#Define output directory:
output_dir <- "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/HOMER"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define file paths on your computer
peak_files <- c(
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep1/WTrep1_peaks.narrowPeak", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep2/WTrep2_peaks.narrowPeak", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep3/WTrep3_peaks.narrowPeak"
)

bam_files <- c(
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep1/WTrep1_final.bam", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep2/WTrep2_final.bam", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep3/WTrep3_final.bam"
)

# Define sample metadata
sample_names <- c("WT1","WT2","WT3")
treatment <- c("WT", "WT", "WT")

# Create the sample sheet data frame
sample_sheet <- data.frame(
  SampleID = sample_names,       # Unique identifier for each sample
  Condition = treatment,         # Experimental condition/group
  bamReads = bam_files,          # Path to aligned reads (BAM)
  Peaks = peak_files,            # Path to called peaks
  PeakCaller = "narrow"          # Type of peak files (narrow or broad)
)

# Optional: Save the sample sheet for future reference
write.csv(sample_sheet, file.path(output_dir, "atac_seq_sample_sheet.csv"), row.names = FALSE)


# =====================================================
# STEP 2: Create and Process DiffBind Object
# =====================================================

# Create the DBA object from the sample sheet
atac_wt_db <- dba(sampleSheet = sample_sheet)
 
# Examine the DBA object
dba.show(atac_wt_db)

# =====================================================
# STEP 3: Counting Reads in Peaks
# =====================================================

# Count reads overlapping peaks
atac_wt_db <- dba.count(atac_wt_db,
                     summits = 250,       # Center peaks on summits, extend +/- 250bp
                     minOverlap = 2,      # Require at least 2 samples to share a peak
                     bUseSummarizeOverlaps = TRUE)  # Use GenomicRanges::summarizeOverlaps()
 
# View updated DBA object with count information
dba.show(atac_wt_db)
 
# Extract the raw count matrix (optional)
counts <- dba.peakset(atac_wt_db, bRetrieve = TRUE)
head(counts)

# Extract consensus peaks (shared in ≥2 WT samples)
consensus <- dba.peakset(atac_wt_db, bRetrieve = TRUE)


# Convert to BED (0-based, half-open)
consensus_df <- data.frame(
  seqnames = as.character(seqnames(consensus)),
  start_bed = pmax(start(consensus) - 1, 0),  # 0-based start
  end = end(consensus)
)


# =====================================================
# STEP 4: Write BED file
# =====================================================

bed_out <- file.path(output_dir, "WT_consensus.bed")
fwrite(consensus_df,
       bed_out,
       sep = "\t",
       col.names = FALSE)

cat("Consensus BED saved to:", bed_out, "\n")

EOF

conda deactivate
conda activate env_3

# =====================================================
# STEP 2 — Paths to differential BED files
# =====================================================
WT_BG="$OUTPUT_DIR/WT_consensus.bed"
KO1_UP="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/KO1_up.bed"
KO1_DOWN="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/KO1_down.bed"
KO2_UP="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/KO2_up.bed"
KO2_DOWN="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/KO2_down.bed"
CON_UP="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/consensus_up.bed"
CON_DOWN="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/consensus_down.bed"

# =====================================================
# STEP 3 — Run HOMER motif finding
# =====================================================
mkdir -p $OUTPUT_DIR

run_homer () {
    INPUT=$1
    NAME=$2

    findMotifsGenome.pl "$INPUT" hg38 "$OUTPUT_DIR/$NAME" \
        -bg "$WT_BG" \
        -size 200 \
        -p 60
}

run_homer "$KO1_UP"   "KO1_UP_motifs"
run_homer "$KO1_DOWN" "KO1_DOWN_motifs"
run_homer "$KO2_UP"   "KO2_UP_motifs"
run_homer "$KO2_DOWN" "KO2_DOWN_motifs"
run_homer "$CON_UP"   "Consensus_UP_motifs"
run_homer "$CON_DOWN" "Consensus_DOWN_motifs"

echo "[INFO] 🎉 All HOMER motif analyses completed successfully."
