# single-cell-influenza-vaccine-response
Exploratory Scanpy and Seurat analysis of single-cell RNA-seq data (GSE211560) to investigate immune-cell composition and transcriptional responses to MF59-adjuvanted influenza vaccination in older adults.
# Project overview
This notebook documents an exploratory analysis of publicly available single-cell RNA-sequencing data (GEO: GSE211560) from PBMCs of older adults (65–80 years) vaccinated with the MF59-adjuvanted influenza vaccine (Fluad). The dataset comprises six donors (three triple responders and three triple non-responders) sampled at baseline (Day 0) and Day 7 post-vaccination (12 total samples), with PBMCs re-stimulated in vitro with vaccine antigens prior to 10x Genomics 5′ sequencing. [Original Article Link]([url](https://www.nature.com/articles/s41467-022-34487-z))

The workflow covers quality control, normalization, dimensionality reduction, clustering, marker-based cell annotation, cell-type composition, and focused B-cell, T-cell, and NK-cell analyses. A separate Seurat workflow was used as an independent check of the broad cellular structure.
