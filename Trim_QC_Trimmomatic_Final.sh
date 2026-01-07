#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=256G
#SBATCH --partition=HTC
#SBATCH --mail-user=jah550@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load modules
module purge
module load trimmomatic/0.38

# Activate Conda environment
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate env_3

# Define paths
DATA_PATH=/ix1/bnacev/JJH/Srishti/Srishti_LMS_ATAC/01.RawData
THREADS=60
TRIM_FA=$HOME/trimmomatic_adapters/NexteraPE-PE.fa

for folder in "$DATA_PATH"/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        echo "Processing sample: $folder_name"

        R1_FILE="$folder/${folder_name}_CKDL250025412-1A_22YY7NLT4_L6_1.fq.gz"
        R2_FILE="$folder/${folder_name}_CKDL250025412-1A_22YY7NLT4_L6_2.fq.gz"

        # Pre-trim FastQC
        fastqc --threads $THREADS "$R1_FILE" "$R2_FILE" -o "$folder"

        # Trimmomatic PE
        trimmomatic PE -threads $THREADS -phred33 \
            "$R1_FILE" "$R2_FILE" \
            "$folder/${folder_name}_R1_paired.fq.gz" "$folder/${folder_name}_R1_unpaired.fq.gz" \
            "$folder/${folder_name}_R2_paired.fq.gz" "$folder/${folder_name}_R2_unpaired.fq.gz" \
            ILLUMINACLIP:"$TRIM_FA":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        # Post-trim FastQC
        fastqc --threads $THREADS \
               "$folder/${folder_name}_R1_paired.fq.gz" "$folder/${folder_name}_R2_paired.fq.gz" \
               -o "$folder"
    fi
done
