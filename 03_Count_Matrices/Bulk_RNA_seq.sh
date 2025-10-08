#!/bin/bash
# RNA-seq Analysis Pipeline Script
# Author: Javeria Farooq
# Description: Automated pipeline for RNA-seq data processing including QC, trimming, alignment, and count generation.

set -e  # Exit immediately if a command exits with a non-zero status
exec > >(tee -i pipeline.log)
exec 2>&1

# Step 1: Install Dependencies
sudo apt update && sudo apt upgrade -y
sudo apt install fastqc fastp hisat2 samtools sra-toolkit subread wget curl default-jdk -y

# Step 2: Create Directory Structure
mkdir -p data/raw data/trimmed alignment/sam alignment/bam qc_reports counts reference

# Step 3: Check for HISAT2 Index
if [ ! -f data/reference/genome.1.ht2 ]; then
    echo "Building HISAT2 index..."
    hisat2-build data/reference/genome.fna data/reference/genome
else
    echo "HISAT2 index already exists. Skipping build."
fi

# Step 4: Main Loop for Processing Samples
for ID in ERR5060645 ERR5060646
do
    echo "Processing sample: $ID"

    # Download Data
    prefetch $ID
    fastq-dump --split-files ${ID}.sra -O data/raw/

    # Quality Control
    fastqc data/raw/${ID}_1.fastq data/raw/${ID}_2.fastq -o qc_reports/

    # Trimming
    fastp -i data/raw/${ID}_1.fastq -I data/raw/${ID}_2.fastq \
          -o data/trimmed/${ID}_1_trim.fastq.gz -O data/trimmed/${ID}_2_trim.fastq.gz

    # Alignment
    hisat2 -p 4 -x data/reference/genome \
      -1 data/trimmed/${ID}_1_trim.fastq.gz -2 data/trimmed/${ID}_2_trim.fastq.gz \
      -S alignment/sam/${ID}.sam

    # Convert and Sort
    samtools view -bS alignment/sam/${ID}.sam > alignment/bam/${ID}.bam
    samtools sort alignment/bam/${ID}.bam -o alignment/bam/${ID}_sorted.bam
    samtools index alignment/bam/${ID}_sorted.bam

    echo "Finished processing $ID"
done

# Step 5: Generate Read Count Matrix
featureCounts -T 4 -a data/reference/annotation.gtf -o counts/gene_counts.txt alignment/bam/*_sorted.bam

echo "RNA-seq pipeline completed successfully. Results saved in 'counts/gene_counts.txt'."
