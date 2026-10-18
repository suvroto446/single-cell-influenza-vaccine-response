# single-cell-influenza-vaccine-response
Exploratory Scanpy and Seurat analysis of single-cell RNA-seq data (GSE211560) to investigate immune-cell composition and transcriptional responses to MF59-adjuvanted influenza vaccination in older adults.
# Project overview
This notebook documents an exploratory analysis of publicly available single-cell RNA-sequencing data (GEO: GSE211560) from PBMCs of older adults (65–80 years) vaccinated with the MF59-adjuvanted influenza vaccine (Fluad). The dataset comprises six donors (three triple responders and three triple non-responders) sampled at baseline (Day 0) and Day 7 post-vaccination (12 total samples), with PBMCs re-stimulated in vitro with vaccine antigens prior to 10x Genomics 5′ sequencing. Original article: [https://www.nature.com/articles/s41467-022-34487-z](https://doi.org/10.1038/s41467-022-34487-z)

The workflow covers quality control, normalization, dimensionality reduction, clustering, marker-based cell annotation, cell-type composition, and focused B-cell, T-cell, and NK-cell analyses. A separate Seurat workflow was used as an independent check of the broad cellular structure.

# Research questions
1. Immune-cell composition: How do the composition and relative abundance of major PBMC populations differ between vaccine responders and non-responders at baseline and Day 7?
2. B-cell responses: What B-cell subpopulations and transcriptional programs are observed after vaccination, and how do their expression patterns vary by response group and time point?
3. T-cell responses: What transcriptional changes occur in T-cell populations between baseline and Day 7, and which patterns differ between responders and non-responders?
4. NK-cell responses: How do NK-cell transcriptional states, including genes associated with cytotoxicity and interferon responses, vary across time points and response groups?
5. Analytical validation: How consistent is the broad cellular structure identified using Scanpy with that obtained through an independent Seurat workflow?
