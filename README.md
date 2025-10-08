# Bulk RNA-seq Analysis Repository

This repository contains quality control reports, count matrices, scripts, and visualizations from bulk RNA-seq analysis.

## Repository Structure
- **01_QC_Reports/**: FastQC HTML reports for raw sequencing reads.
- **02_QC_Plots/**: PNG visualizations for QC metrics.
- **03_Count_Matrices/**: Transcript and gene count matrices with scripts (`prepDE.py` and `Bulk_RNA_seq.sh`) for generating them.

## Tools Used
- FastQC for quality control
- Python (`prepDE.py`) for count matrix extraction
- Bash (`Bulk_RNA_seq.sh`) for workflow automation
- Visualization for DEGS performed in IDEp
