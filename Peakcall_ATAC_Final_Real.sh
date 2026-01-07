#!/bin/bash
#SBATCH --time=12:00:00           # Request 12 hours of runtime
#SBATCH --ntasks=1                # Request 1 task
#SBATCH --cpus-per-task=60        # Request 60 CPU cores
#SBATCH --mem=256G                # Request 256 GB of memory
#SBATCH --partition=HTC           # Partition to use
#SBATCH --mail-user=jah550@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load Picard module
module purge
module load picard/3.0.0

# Source your Miniconda installation
source $HOME/miniconda/etc/profile.d/conda.sh

# Activate environment
conda activate env_3

# Define paths
DATA_PATH=/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Temporary
INDEX_PATH=/ix1/bnacev/JJH/ref/hg38/hg38
BLACKLIST=/ix1/bnacev/JJH/BlacklistFile/hg38-blacklist.v2.bed
THREADS=60
OUTPUT_DIR=/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed
mkdir -p "$OUTPUT_DIR"

# === Step 2: Convert BAM → BigWig ===
for folder in "$DATA_PATH"/*; do
    SAMPLE=$(basename "$folder")
    mkdir -p "$OUTPUT_DIR/$SAMPLE"

    BAM_FILE="$folder/${SAMPLE}_final.bam"
    BIGWIG_FILE="$OUTPUT_DIR/$SAMPLE/${SAMPLE}.bw"

    echo "Generating BigWig for $SAMPLE..."
    bamCoverage -b "$BAM_FILE" \
        -o "$BIGWIG_FILE" \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --centerReads \
        -p "$THREADS"
done

# === Step 3: Peak calling with MACS2 ===
for folder in "$DATA_PATH"/*; do
    SAMPLE=$(basename "$folder")
    BAM_FILE="$folder/${SAMPLE}_final.bam"
    MACS2_PREFIX="$OUTPUT_DIR/$SAMPLE/${SAMPLE}"

    echo "Calling peaks for $SAMPLE..."
    macs2 callpeak \
        -t "$BAM_FILE" \
        -f BAMPE \
        -g hs \
        -n "$MACS2_PREFIX" \
        --outdir "$OUTPUT_DIR/$SAMPLE" \
        --nomodel \
        --nolambda \
        --keep-dup all \
        --call-summits
done

# === Step 4: Calculate FRiP scores ===
FRIP_FILE="$OUTPUT_DIR/FRiP_scores.txt"
echo "Sample\tFRiP" > "$FRIP_FILE"

for folder in "$DATA_PATH"/*; do
    SAMPLE=$(basename "$folder")
    BAM_FILE="$folder/${SAMPLE}_final.bam"
    PEAK_FILE="$OUTPUT_DIR/$SAMPLE/${SAMPLE}_peaks.narrowPeak"

    if [ -f "$PEAK_FILE" ]; then
        TOTAL_READS=$(samtools view -c "$BAM_FILE")
        READS_IN_PEAKS=$(bedtools intersect -u -a <(samtools view -b "$BAM_FILE") -b "$PEAK_FILE" | samtools view -c)
        FRIP=$(echo "scale=4; $READS_IN_PEAKS / $TOTAL_READS" | bc)
        echo -e "$SAMPLE\t$FRIP" >> "$FRIP_FILE"
    else
        echo -e "$SAMPLE\tNA" >> "$FRIP_FILE"
    fi
done

# === Step 5: Clean up intermediate BAMs ===
for folder in "$DATA_PATH"/*; do
    SAMPLE=$(basename "$folder")
    echo "Removing intermediate BAMs for $SAMPLE..."
    rm -f "$folder/"*_sorted.bam
    rm -f "$folder/"*_dedup.bam
    rm -f "$folder/"*_dedup_filtered.bam
done

echo "All processing complete."
