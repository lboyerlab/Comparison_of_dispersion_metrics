# Comparison_dispersion_metrics
Supporting information for manuscript submitted to PLoS Computational Biology titled "Comparison of dispersion metrics for estimating transcriptional noise in single-cell RNA-seq data and applications to cardiomyocyte biology"

"Simulate Counts.R" provides code for simulating single-cell counds from the Poisson, negative binomial, beta-Poisson, and uniform distributions, with varying levels of dispersion and dataset size. This script calculates dispersion in the resulting simulated counts matrices using the Gini index, variance-to-mean ratio (VMR; or Fano factor), variance, and Shannon entropy and saves them as csv files. 

"Generate Simulation Heatmaps.R" reads the simulation results from "Simulate Counts.R" and generates heatmaps of the raw and normalized metrics. This script generated Figure 2 and most of the supplemental figures. It also includes code for visualizing probability distribution functions (Figure 1) and for investigating the paradoxical behavior of the Gini index (Figure 3). 

"scRNA-seq Preprocessing.R" is the script used to preprocess UMI counts from the two publicly available datasets analyzed in the manuscript. It also identifies the cells that are cardiomyocytes, and saves them as a separate Seurat object. 

"MatHG Transcriptional Variability Analysis.R" is the script used to calculate transcriptional variability in the scRNA-seq dataset from Manivannan et al. (2022) (GSE193746). This script calculates the change in VMR between matHG and control conditions, and calculates and visualizes correlations between absolute change in VMR and gene characteristics (Figure 4). It also performs GSEA and TF motif enrichment analysis using [RcisTarget](https://www.bioconductor.org/packages/release/bioc/html/RcisTarget.html). This script generated Figure 5 and S1 Table. 

"T21 Transcriptional Variability Analysis.R" is similar to "MatHG Transcriptional Variability Analysis.R". It is the script that was used to analyze transcriptional variability from Lana-Elola et al. (2024) (GSE196447). It calculates change in transcriptional variability between T21 and control conditions, checks correlations between change in VMR and gene characteristics, and performs GSEA and TF motif enrichment analysis. This script generated S2 Table. 

