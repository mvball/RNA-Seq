## Project Objectives:
This repository contains the analysis and findings from an RNA sequencing (RNA-Seq) study on breast cancer subtypes, focusing on Triple-Negative Breast Cancer (TNBC), HER2-positive, and non-TNBC. The study aims to characterize gene expression profiles and identify differentially expressed genes (DEGs) based on publicly available data from Eswaran et al., 2012 (The raw RNA-Seq data can be downloaded from NCBI GEO using the accession number GSE52194).

## Main Objectives:
1. Obtain RNA-Seq data from the Gene Expression Omnibus (GEO).
2. Perform quality control (QC) using FastQC and MultiQC.
3. Align sequencing reads to the human reference genome (GRCh38) using HISAT2.
4. Process alignments (SAM → BAM conversion, sorting, indexing) using Samtools.
5. Quantify gene expression using featureCounts.
6. Conduct differential expression analysis using DESeq2 in R.
7. Perform functional enrichment analysis (Gene Ontology) and visualize results.

## Installation & Dependencies
To run this pipeline, ensure you have the following software installed:

Command-line tools:

1. FastQC (v0.11.9) – Quality control of raw reads.
2. MultiQC (v1.19) – Summarizes FastQC reports.
3. Fastp (v0.20.0) – Trims adapters and filters low-quality reads.
4. HISAT2 (v2.1.0) – Aligns RNA-Seq reads to the genome.
5. Samtools (v1.21) – Processes alignment files.
6. FeatureCounts (v2.0.7) – Counts gene-level expression.

R Packages:
1. DESeq2 (v3.20) – Differential gene expression analysis.
2. ggplot2 – Visualization of results.
3. Pheatmap – Heatmaps of DEGs.
4. clusterProfiler – Gene Ontology (GO) enrichment.

All scripts contain commentary on each step and parameters used. 

## Reference
Eswaran J, et al. Transcriptomic landscape of breast cancers through mRNA sequencing. Sci Rep 2, 264 (2012). https://doi.org/10.1038/srep00264