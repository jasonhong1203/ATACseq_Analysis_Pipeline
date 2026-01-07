#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=256G
#SBATCH --partition=HTC
#SBATCH --mail-user=jah550@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Conda environment
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate env_3
THREADS=60

#-----------------------------------------------
# STEP 1: Copy and merge bam and bed files 
#-----------------------------------------------

#Make directories for merged bam and peak files 

mkdir -p /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/WT \
         /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO1 \
         /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO2 \
         /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT \
         /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1 \
         /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2

#Copy BAM files for the WT to the corresponding folder. 

cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep1/WTrep1_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/WT
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep2/WTrep2_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/WT
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/WTrep3/WTrep3_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/WT

#Copy BAM files for KO1 to the corresponding folder 

cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO1rep1/KO1rep1_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO1
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO1rep2/KO1rep2_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO1
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO1rep3/KO1rep3_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO1

#Copy BAM files for KO2 to the corresponding folder 

cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO2rep1/KO2rep1_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO2
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO2rep2/KO2rep2_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO2
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData/KO2rep3/KO2rep3_final.bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO2

#Copy peak files for WT to the corresponding folder 

cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep1/WTrep1_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep2/WTrep2_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/WTrep3/WTrep3_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT

#Copy peak files for KO1 to the corresponding folder 

cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep1/KO1rep1_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep2/KO1rep2_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO1rep3/KO1rep3_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1

#Copy peak files for KO2 to the corresponding folder 

cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep1/KO2rep1_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep2/KO2rep2_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2
cp /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/KO2rep3/KO2rep3_peaks.narrowPeak /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2

#Merge BAM files for the WT 
cd /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/WT
samtools merge -@ $THREADS WT_merged.bam WTrep*_final.bam
rm WTrep*_final.bam  # Remove individual files to save space

# Merge BAM files for KO1
cd /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO1
samtools merge -@ $THREADS KO1_merged.bam KO1rep*_final.bam
rm KO1rep*_final.bam

# Merge BAM files for KO2
cd /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO2
samtools merge -@ $THREADS KO2_merged.bam KO2rep*_final.bam
rm KO2rep*_final.bam

# Merge peak files for WT
cd /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT
cat WTrep*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > WT_merged.narrowPeak
rm WTrep*_peaks.narrowPeak  # Remove individual files to save space

# Merge peak files for KO1
cd /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1
cat KO1rep*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > KO1_merged.narrowPeak
rm KO1rep*_peaks.narrowPeak  # Remove individual files to save space

# Merge peak files for KO2
cd /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2
cat KO2rep*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > KO2_merged.narrowPeak
rm KO2rep*_peaks.narrowPeak  # Remove individual files to save space

#-----------------------------------------------
# STEP 2: Correct sequence bias in ATAC-seq data 
#-----------------------------------------------

conda deactivate
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate TOBIAS_ENV

mkdir -p /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/WT
mkdir -p /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO1
mkdir -p /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO2

# Correct bias for WT
TOBIAS ATACorrect \
    --bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/WT/WT_merged.bam \
    --genome /ix1/bnacev/JJH/ref/hg38/hg38.fa \
    --peaks /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT/WT_merged.narrowPeak \
    --blacklist /ix1/bnacev/JJH/BlacklistFile/hg38-blacklist.v2.bed \
    --outdir /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/WT \
    --cores $THREADS

# Correct bias for KO1
TOBIAS ATACorrect \
    --bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO1/KO1_merged.bam \
    --genome /ix1/bnacev/JJH/ref/hg38/hg38.fa \
    --peaks /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1/KO1_merged.narrowPeak \
    --blacklist /ix1/bnacev/JJH/BlacklistFile/hg38-blacklist.v2.bed \
    --outdir /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO1 \
    --cores $THREADS

# Correct bias for KO2
TOBIAS ATACorrect \
    --bam /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/bam/KO2/KO2_merged.bam \
    --genome /ix1/bnacev/JJH/ref/hg38/hg38.fa \
    --peaks /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2/KO2_merged.narrowPeak \
    --blacklist /ix1/bnacev/JJH/BlacklistFile/hg38-blacklist.v2.bed \
    --outdir /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO2 \
    --cores $THREADS

#-----------------------------------------------
# STEP 3: Calculate footprinting scores
#-----------------------------------------------

# Calculate scores for WT
TOBIAS FootprintScores \
    --signal /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/WT/WT_merged_corrected.bw \
    --regions /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT/WT_merged.narrowPeak \
    --output /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/WT/WT_footprint.bw \
    --cores $THREADS
 
# Calculate scores for KO1
TOBIAS FootprintScores \
    --signal /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO1/KO1_merged_corrected.bw \
    --regions /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1/KO1_merged.narrowPeak \
    --output /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO1/KO1_footprint.bw \
    --cores $THREADS

# Calculate scores for KO2
TOBIAS FootprintScores \
    --signal /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO2/KO2_merged_corrected.bw \
    --regions /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2/KO2_merged.narrowPeak \
    --output /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO2/KO2_footprint.bw \
    --cores $THREADS

#------------------------------------------------------------------------------
# STEP 4a: Prepare unified peak set for comparative analysis between WT and KO1
#------------------------------------------------------------------------------
 
# Create a directory for comparative analysis
mkdir -p /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1
 
# Combine peaks from both conditions into a unified set
cat /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT/WT_merged.narrowPeak \
    /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO1/KO1_merged.narrowPeak | \
    sort -k1,1 -k2,2n | \
    bedtools merge > /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined.narrowPeak
 
# Download genome annotation file (human)
cd /ix1/bnacev/JJH/ref/genome_annotation
if [ ! -f gencode.v48.annotation.gtf ]; then
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
    gunzip gencode.v48.annotation.gtf.gz
fi

 
# Annotate combined peaks with nearby genes
uropa \
    --bed /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined.narrowPeak \
    --gtf /ix1/bnacev/JJH/ref/genome_annotation/gencode.v48.annotation.gtf \
    --show_attributes gene_id gene_name \
    --feature_anchor start \
    --distance 20000 10000 \
    --feature gene \
    -o /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1
 
# Extract header for the annotated peaks
cut -f 1-6,16-17 /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits.txt | \
    head -n 1 > /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits_header.txt
 
# Match the header and the annotated peaks
cut -f 1-6,16-17 /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits.txt > \
/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits_2.txt
 
# Filter for standard chromosomes only (removes contigs and alternative assemblies)
grep -E '^chr[0-9XY]+\s' /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits_2.txt > \
    /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits_3.bed

#------------------------------------------------------------------------------
# STEP 4b: Prepare unified peak set for comparative analysis between WT and KO2
#------------------------------------------------------------------------------

# Create a directory for comparative analysis
mkdir -p /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2

# Combine peaks from both conditions into a unified set
cat /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/WT/WT_merged.narrowPeak \
    /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/merged/peaks/KO2/KO2_merged.narrowPeak | \
    sort -k1,1 -k2,2n | \
    bedtools merge > /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined.narrowPeak

# Annotate combined peaks with nearby genes
uropa \
    --bed /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined.narrowPeak \
    --gtf /ix1/bnacev/JJH/ref/genome_annotation/gencode.v48.annotation.gtf \
    --show_attributes gene_id gene_name \
    --feature_anchor start \
    --distance 20000 10000 \
    --feature gene \
    -o /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2

# Extract header for the annotated peaks
cut -f 1-6,16-17 /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits.txt | \
    head -n 1 > /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits_header.txt

# Match the header and the annotated peaks
cut -f 1-6,16-17 /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits.txt > \
/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits_2.txt

# Filter for standard chromosomes only (removes contigs and alternative assemblies)
grep -E '^chr[0-9XY]+\s' /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits_2.txt > \
/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits_3.bed

#------------------------------------------------------------
# STEP 5a: Detect differential TF binding between WT and KO1
#------------------------------------------------------------
 
# Identify differential TF binding sites between WT and KO1
TOBIAS BINDetect \
    --motifs /ix1/bnacev/JJH/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
    --signals /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/WT/WT_footprint.bw \
             /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO1/KO1_footprint.bw \
    --genome /ix1/bnacev/JJH/ref/hg38/hg38.fa \
    --peaks /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits_3.bed \
    --peak_header /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1/WT_KO1_combined_finalhits_header.txt \
    --outdir /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO1 \
    --cond_names WT KO1 \
    --cores $THREADS

#------------------------------------------------------------
# STEP 5b: Detect differential TF binding between WT and KO2
#------------------------------------------------------------

# Identify differential TF binding sites between WT and KO2
TOBIAS BINDetect \
    --motifs /ix1/bnacev/JJH/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
    --signals /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/WT/WT_footprint.bw \
             /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/atacorrect/KO2/KO2_footprint.bw \
    --genome /ix1/bnacev/JJH/ref/hg38/hg38.fa \
    --peaks /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits_3.bed \
    --peak_header /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2/WT_KO2_combined_finalhits_header.txt \
    --outdir /ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/Processed/Tobias/WT_KO2 \
    --cond_names WT KO2 \
    --cores $THREADS

echo "[INFO] TOBIAS motif analysis complete"