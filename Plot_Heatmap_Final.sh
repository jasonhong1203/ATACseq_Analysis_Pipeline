#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=128G
#SBATCH --partition=HTC
#SBATCH --mail-user=jah550@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# -----------------------------
# User variables
# -----------------------------
OUTPUT_DIR="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Heatmap"
cd "$OUTPUT_DIR" || exit 1

# Activate environment
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate env_3

THREADS=60

# Input files
TSS_BED="/ix1/bnacev/JJH/Promoter_bed/hg38_promoters_final.bed"
UP_BED="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/consensus_up.bed"
DOWN_BED="/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/diffbind/consensus_down.bed"

# -----------------------------
# Step 1: Merge replicates into averaged BigWigs
# -----------------------------
echo "[INFO] Averaging replicate BigWigs..."

bigwigAverage \
  -b /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep1/WTrep1.bw \
     /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep2/WTrep2.bw \
     /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep3/WTrep3.bw \
  --binSize 50 --skipNonCoveredRegions --numberOfProcessors $THREADS \
  --outFileFormat bigwig -o WT_merged.bw

bigwigAverage \
  -b /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep1/KO1rep1.bw \
     /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep2/KO1rep2.bw \
     /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep3/KO1rep3.bw \
  --binSize 50 --skipNonCoveredRegions --numberOfProcessors $THREADS \
  --outFileFormat bigwig -o KO1_merged.bw

bigwigAverage \
  -b /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep1/KO2rep1.bw \
     /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep2/KO2rep2.bw \
     /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep3/KO2rep3.bw \
  --binSize 50 --skipNonCoveredRegions --numberOfProcessors $THREADS \
  --outFileFormat bigwig -o KO2_merged.bw

echo "[INFO] Done averaging replicates."

# -----------------------------
# Step 2: Compute matrices(±3 kb around center)
# -----------------------------
echo "[INFO] Computing computeMatrix"

computeMatrix reference-point \
  -S WT_merged.bw KO1_merged.bw KO2_merged.bw \
  -R $UP_BED \
  --referencePoint center \
  --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
  --skipZeros -p "$THREADS" \
  -o KO_up_matrix.gz \
  --outFileSortedRegions KO_up_regions.bed

computeMatrix reference-point \
  -S WT_merged.bw KO1_merged.bw KO2_merged.bw \
  -R $DOWN_BED \
  --referencePoint center \
  --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
  --skipZeros -p "$THREADS" \
  -o KO_down_matrix.gz \
  --outFileSortedRegions KO_down_regions.bed

# -----------------------------
# Step 3: Relabel groups and merge
# -----------------------------
echo "[INFO] Relabeling and merging matrices..."

computeMatrixOperations relabel -m KO_up_matrix.gz -o KO_up_labeled.mat.gz --groupLabels Upregulated
computeMatrixOperations relabel -m KO_down_matrix.gz -o KO_down_labeled.mat.gz --groupLabels Downregulated
computeMatrixOperations rbind -m KO_up_labeled.mat.gz KO_down_labeled.mat.gz -o KO_combined_matrix.mat.gz

# -----------------------------
# Step 4a: Plot heatmap
# -----------------------------
echo "[INFO] Generating heatmap..."

plotHeatmap -m KO_combined_matrix.mat.gz \
  -out KO_updown_heatmap_final.pdf \
  --colorList "white,red" \
  --samplesLabel WT KO1 KO2 \
  --sortRegions descend --sortUsing mean \
  --zMin 0 --zMax 30 \
  --missingDataColor "#FFFFFF" \
  --heatmapHeight 12 --heatmapWidth 6 \
  --refPointLabel "Center" \
  --plotTitle "ATAC-seq centered ±3kb Heatmap"

# -----------------------------
# Step 4b: Plot average profile
# -----------------------------
echo "[INFO] Generating average profile..."

plotProfile -m KO_combined_matrix.mat.gz \
  -out KO_updown_profile_final.pdf \
  --perGroup \
  --colors "black" "red" "darkred" \
  --samplesLabel WT KO1 KO2 \
  --refPointLabel "Center" \
  --plotTitle "ATAC-seq centered ±3kb Profile" \
  --averageType mean

echo "[INFO] DONE. Heatmaps and profiles are in $OUTPUT_DIR"
