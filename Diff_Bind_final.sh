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
output_dir <- "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define file paths on your computer
peak_files <- c(
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep1/WTrep1_peaks.narrowPeak", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep2/WTrep2_peaks.narrowPeak", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep3/WTrep3_peaks.narrowPeak", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep1/KO1rep1_peaks.narrowPeak",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep2/KO1rep2_peaks.narrowPeak",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep3/KO1rep3_peaks.narrowPeak",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep1/KO2rep1_peaks.narrowPeak",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep2/KO2rep2_peaks.narrowPeak",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep3/KO2rep3_peaks.narrowPeak"
)

bam_files <- c(
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep1/WTrep1_final.bam", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep2/WTrep2_final.bam", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep3/WTrep3_final.bam", 
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO1rep1/KO1rep1_final.bam",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO1rep2/KO1rep2_final.bam",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO1rep3/KO1rep3_final.bam",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO2rep1/KO2rep1_final.bam",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO2rep2/KO2rep2_final.bam",
  "/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO2rep3/KO2rep3_final.bam"
)

# Define sample metadata
sample_names <- c("WT1","WT2","WT3","KO1_1","KO1_2","KO1_3","KO2_1","KO2_2","KO2_3")
treatment <- c("WT", "WT", "WT", "KO1", "KO1", "KO1", "KO2", "KO2", "KO2")

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
atac_db <- dba(sampleSheet = sample_sheet)
 
# Examine the DBA object
dba.show(atac_db)

# =====================================================
# STEP 3: Counting Reads in Peaks
# =====================================================

# Count reads overlapping peaks
atac_db <- dba.count(atac_db,
                     summits = 250,       # Center peaks on summits, extend +/- 250bp
                     minOverlap = 2,      # Require at least 2 samples to share a peak
                     bUseSummarizeOverlaps = TRUE)  # Use GenomicRanges::summarizeOverlaps()
 
# View updated DBA object with count information
dba.show(atac_db)
 
# Extract the raw count matrix (optional)
counts <- dba.peakset(atac_db, bRetrieve = TRUE)
head(counts)

# =====================================================
# STEP 4: Exploratory Data Analysis
# =====================================================
# Generate correlation heatmap
pdf(file.path(output_dir, "correlation_heatmap.pdf"), width = 8, height = 8)
dba.plotHeatmap(atac_db)
dev.off()


# Create PCA plot to visualize sample relationships
pdf(file.path(output_dir, "PCA_plot.pdf"), width = 8, height = 8)
dba.plotPCA(atac_db, DBA_CONDITION, label = DBA_ID)
dev.off()

# =====================================================
# STEP 5: Defining Contrasts 
# =====================================================

# WT vs KO1
atac_db <- dba.contrast(atac_db, 
                        categories = DBA_CONDITION,
                        minMembers = 3,
                        contrast = c("Condition", "WT", "KO1"))
# WT vs KO2
atac_db <- dba.contrast(atac_db, 
                        categories = DBA_CONDITION,
                        minMembers = 3,
                        contrast = c("Condition", "WT", "KO2"))

dba.show(atac_db, bContrasts = TRUE)

# =====================================================
# STEP 6: Running Differential Analysis
# =====================================================

#Run differential analysis 
atac_db <- dba.analyze(atac_db)
dba.show(atac_db, bContrasts = TRUE)

# =====================================================
# STEP 7a: Retrieving Differential Peaks for KO1 
# =====================================================

#Get significant results (FDR < 0.05, log fold change > 1)

diff_peaks_WT_vs_KO1 <- dba.report(atac_db, contrast = 1, th = 0.05, fold = 1)
length(diff_peaks_WT_vs_KO1)
diff_df1 <- as.data.frame(diff_peaks_WT_vs_KO1)
fwrite(diff_df1, file.path(output_dir, "differential_peaks_WT_vs_KO1.csv"))

#Annotate all significant peaks
peak_anno <- annotatePeak(diff_peaks_WT_vs_KO1,
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb = "org.Hs.eg.db",
                          tssRegion = c(-3000, 3000),
                          verbose = FALSE)

#Convert annotated object to data.frame
peak_df <- as.data.frame(peak_anno)

#Rename columns to match your volcano plot code
peak_df$logFC <- peak_df$Fold
peak_df$adj.P.Val <- peak_df$FDR
peak_df$Gene <- peak_df$SYMBOL

#Thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Assign regulation
peak_df$diffexpressed <- "Not Significant"
peak_df$diffexpressed[peak_df$logFC < -lfc_cutoff & peak_df$adj.P.Val < padj_cutoff] <- "KO Up-regulated"
peak_df$diffexpressed[peak_df$logFC > lfc_cutoff & peak_df$adj.P.Val < padj_cutoff] <- "KO Down-regulated"

#Save annotated table
fwrite(peak_df, file.path(output_dir, "differential_peaks_WT_vs_KO1_ANNOTATED.csv"))

#Save upregulated/downregulated separately

fwrite(
        peak_df[peak_df$diffexpressed == "KO Up-regulated",],
        file.path(output_dir, "differential_KO1_up.csv")
       )
fwrite(
        peak_df[peak_df$diffexpressed == "KO Down-regulated",],
        file.path(output_dir, "differential_KO1_down.csv")
       )

# =====================================================
# Volcano Plot — KO1 (Top 40 unique genes up & down)
# =====================================================

# Add ranks for combined ranking
peak_df <- peak_df %>%
  mutate(
    rank_padj = rank(adj.P.Val, ties.method = "average"),
    rank_logFC = rank(-abs(logFC), ties.method = "average"),
    combined_rank = (rank_padj + rank_logFC) / 2,
    peak_id = row_number()
  )

# Function: pick top N unique genes
get_top_unique_genes <- function(df, diff_type, n = 40) {
  df %>%
    filter(diffexpressed == diff_type, !is.na(Gene), Gene != "") %>%
    group_by(Gene) %>%
    slice_min(combined_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(combined_rank) %>%
    slice_head(n = n)
}

# Select top genes
top_up_KO1   <- get_top_unique_genes(peak_df, "KO Up-regulated", 40)
top_down_KO1 <- get_top_unique_genes(peak_df, "KO Down-regulated", 40)

# Label only the exact selected peaks
selected_ids_KO1 <- c(top_up_KO1$peak_id, top_down_KO1$peak_id)
peak_df$label <- ifelse(peak_df$peak_id %in% selected_ids_KO1, peak_df$Gene, NA)

# Flip logFC for plotting
peak_df$logFC_plot <- -peak_df$logFC

# Volcano plot
volcano_KO1 <- ggplot(peak_df, aes(x = logFC_plot, y = -log10(adj.P.Val), color = diffexpressed)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c(
    "Not Significant" = "grey",
    "KO Up-regulated" = "blue",
    "KO Down-regulated" = "red"
  )) +
  geom_text_repel(aes(label = label),
                  size = 3, max.overlaps = 20,
                  box.padding = 0.5, segment.color = "grey50", na.rm = TRUE) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  labs(
    title = "Differential Binding: KO1 vs WT",
    subtitle = paste0(
      "Blue: Up (", sum(peak_df$diffexpressed == "KO Up-regulated"), "), ",
      "Red: Down (", sum(peak_df$diffexpressed == "KO Down-regulated"), ")"
    ),
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggsave(file.path(output_dir, "Volcano_Plot_WT_vs_KO1.pdf"),
       volcano_KO1, width = 10, height = 8)



# =====================================================
# STEP 7b: Retrieving Differential Peaks for KO2
# =====================================================

#Get significant results (FDR < 0.05, log fold change > 1)
diff_peaks_WT_vs_KO2 <- dba.report(atac_db, contrast = 2, th = 0.05, fold = 1)
length(diff_peaks_WT_vs_KO2)
diff_df2 <- as.data.frame(diff_peaks_WT_vs_KO2)
fwrite(diff_df2, file.path(output_dir, "differential_peaks_WT_vs_KO2.csv"))

#Annotate all significant peaks
peak_anno2 <- annotatePeak(diff_peaks_WT_vs_KO2,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           annoDb = "org.Hs.eg.db",
                           tssRegion = c(-3000, 3000),
                           verbose = FALSE)

#Convert annotated object to data.frame
peak_df2 <- as.data.frame(peak_anno2)

#Rename columns to match your volcano plot code
peak_df2$logFC <- peak_df2$Fold
peak_df2$adj.P.Val <- peak_df2$FDR
peak_df2$Gene <- peak_df2$SYMBOL

#Thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Assign regulation
peak_df2$diffexpressed <- "Not Significant"
peak_df2$diffexpressed[peak_df2$logFC < -lfc_cutoff & peak_df2$adj.P.Val < padj_cutoff] <- "KO2 Up-regulated"
peak_df2$diffexpressed[peak_df2$logFC > lfc_cutoff & peak_df2$adj.P.Val < padj_cutoff] <- "KO2 Down-regulated"

#Save annotated table
fwrite(peak_df2, file.path(output_dir, "differential_peaks_WT_vs_KO2_ANNOTATED.csv"))

#Save upregulated/downregulated separately
fwrite(
        peak_df2[peak_df2$diffexpressed == "KO2 Up-regulated",],
        file.path(output_dir, "differential_KO2_up.csv")
       )
fwrite(
        peak_df2[peak_df2$diffexpressed == "KO2 Down-regulated",],
        file.path(output_dir, "differential_KO2_down.csv")
       )


# =====================================================
# Volcano Plot — KO2 (Top 40 unique genes up & down)
# =====================================================

# Add ranks
peak_df2 <- peak_df2 %>%
  mutate(
    rank_padj = rank(adj.P.Val, ties.method = "average"),
    rank_logFC = rank(-abs(logFC), ties.method = "average"),
    combined_rank = (rank_padj + rank_logFC) / 2,
    peak_id = row_number()
  )

# Select top unique genes
top_up_KO2   <- get_top_unique_genes(peak_df2, "KO2 Up-regulated", 50)
top_down_KO2 <- get_top_unique_genes(peak_df2, "KO2 Down-regulated", 50)

# Label only the selected peaks
selected_ids_KO2 <- c(top_up_KO2$peak_id, top_down_KO2$peak_id)
peak_df2$label <- ifelse(peak_df2$peak_id %in% selected_ids_KO2, peak_df2$Gene, NA)

# Flip fold change
peak_df2$logFC_plot <- -peak_df2$logFC

# Volcano plot
volcano_KO2 <- ggplot(peak_df2, aes(x = logFC_plot, y = -log10(adj.P.Val), color = diffexpressed)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c(
    "Not Significant" = "grey",
    "KO2 Up-regulated" = "blue",
    "KO2 Down-regulated" = "red"
  )) +
  geom_text_repel(aes(label = label),
                  size = 3, max.overlaps = 50,
                  box.padding = 0.5, segment.color = "grey50", na.rm = TRUE) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  labs(
    title = "Differential Binding: KO2 vs WT",
    subtitle = paste0(
      "Blue: Up (", sum(peak_df2$diffexpressed == "KO2 Up-regulated"), "), ",
      "Red: Down (", sum(peak_df2$diffexpressed == "KO2 Down-regulated"), ")"
    ),
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggsave(file.path(output_dir, "Volcano_Plot_WT_vs_KO2.pdf"),
       volcano_KO2, width = 10, height = 8)



# =====================================================
# STEP 9: Generate BED files for KO1 and KO2 (corrected)
# =====================================================

library(data.table)

create_bed_files <- function(df, prefix, up_label, down_label, output_dir) {
  # Adjust start coordinate for 0-based BED format
  df$start_bed <- df$start - 1
  
  # Up-regulated peaks
  bed_up <- df[df$diffexpressed == up_label, c("seqnames", "start_bed", "end")]
  fwrite(bed_up, file.path(output_dir, paste0(prefix, "_up.bed")), 
         sep = "\t", col.names = FALSE)
  
  # Down-regulated peaks
  bed_down <- df[df$diffexpressed == down_label, c("seqnames", "start_bed", "end")]
  fwrite(bed_down, file.path(output_dir, paste0(prefix, "_down.bed")), 
         sep = "\t", col.names = FALSE)
  
  cat(sprintf("BED files for %s created: %d up, %d down peaks.\n", 
              prefix, nrow(bed_up), nrow(bed_down)))
}

# Generate BED files for KO1
create_bed_files(
  df = peak_df, 
  prefix = "KO1",
  up_label = "KO Up-regulated",
  down_label = "KO Down-regulated",
  output_dir = output_dir
)

# Generate BED files for KO2
create_bed_files(
  df = peak_df2, 
  prefix = "KO2",
  up_label = "KO2 Up-regulated",
  down_label = "KO2 Down-regulated",
  output_dir = output_dir
)
EOF

# =====================================================
# STEP 10: Generate Consensus BED files for KO1 and KO2 
# =====================================================
OUTPUT_DIR="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind"
KO1_UP="$OUTPUT_DIR/KO1_up.bed"
KO1_DOWN="$OUTPUT_DIR/KO1_down.bed"
KO2_UP="$OUTPUT_DIR/KO2_up.bed"
KO2_DOWN="$OUTPUT_DIR/KO2_down.bed"

CONSENSUS_KO_UP="$OUTPUT_DIR/consensus_up.bed"
CONSENSUS_KO_DOWN="$OUTPUT_DIR/consensus_down.bed"

# Generate consensus BED files using bedtools intersect
bedtools intersect -a "$KO1_UP" -b "$KO2_UP" -u > "$CONSENSUS_KO_UP"
bedtools intersect -a "$KO1_DOWN" -b "$KO2_DOWN" -u > "$CONSENSUS_KO_DOWN"

echo "[INFO] ✅ All processing complete."
