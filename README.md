# Single-Cell RNA-seq Analysis of Influenza Vaccine Responses in Older Adults
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

# Methods and Analysis Workflow

1. Environment setup and pilot quality control: Set up the computational environment and examine preliminary quality-control metrics.

2. Data loading and quality control: Load the 12 PBMC samples and apply cell-level filtering based on detected genes and mitochondrial RNA percentage.

3. Normalization and dimensionality reduction: Normalize gene expression, identify highly variable genes, perform principal component analysis (PCA), generate uniform manifold approximation and projection (UMAP) embeddings, and cluster cells using Scanpy.

4. Cell annotation and composition analysis: Annotate broad immune-cell populations using marker-gene expression and compare their relative abundance across samples, time points, and response groups.

5. B-cell subpopulation analysis: Subcluster B cells and investigate their marker expression, transcriptional characteristics, differential gene expression, and relevant biological pathways.

6. T-cell and NK-cell analyses: Investigate cell states, differential gene expression, and gene-set enrichment within T-cell and NK-cell populations.

7. Independent validation: Use a separate Seurat workflow to assess the concordance of the broad cellular structure identified by Scanpy.


# Computational Environment

Python: 3.13.16
Operating environment: Linux (Google Colab)
Package dependencies: See `requirements.txt`


# Key Findings

1. Across both Scanpy and Seurat workflows, marker gene analysis resolved 7 broad lymphocyte populations—B cells, Naive/Central Memory CD4 T cells, Resting CD4 T cells, Naive/central memory mixed T cells, Cytotoxic CD8 T cells, NK cells, and NK/T-like cytotoxic cells—with no myeloid cells identified in the dataset.

2. Cellular composition profiling demonstrates that non-responders are defined by a high baseline frequency of NK cells that remains static post-vaccination. Conversely, protective vaccine responsiveness is driven by a higher baseline reservoir of Naïve/Central Memory CD4 T cells alongside a robust expansion of B cells at Day 7 post-vaccination in Responder 3, while Responder 1 and 2 maintained a base line B-cell frequency.

3. B-cell profiling reveals that protective influenza vaccination response in Responders is characterized by dual activation features: (1) a prominent baseline and post-vaccination enrichment of Class-Switched Plasmablasts/ASCs in specific high-responder phenotypes (representing ~42% of B cells in donor R2) heavily biased toward IGHA1 and IGHM expression, and (2) broad, cohort-wide induction of antigen-presentation capacity mediated by HLA-DR/MHC-II upregulation and translational biogenesis machinery at Day 7 post-vaccination. In contrast, Non-responders fail to maintain class-switched plasmablast populations or execute MHC-II-driven antigen presentation routines.

4. Analysis of the T-cell compartment indicates that protective influenza vaccine responses in Responders are mediated by broad transcriptional activation rather than simple frequency expansion. At Day 7 post-vaccination, Responder T cells execute a robust transcriptional program defined by strong Type-I ISG induction (ISG15, IFITM1, IFI6), enhanced MHC-II assembly, and energetic mitochondrial/translational biogenesis (ATP5ME, EEF1A1). This coordinated activation spans from CD4+ naive/memory pools to terminal CD8+ cytotoxic effectors.

5. Functional characterization of the NK cell compartment reveals that protective influenza vaccination triggers a two-pronged innate response in Responders: (1) robust acute induction of Type-I Interferon-stimulated pathways (ISG15, IFI6, HLA-B/C) at Day 7 post-vaccination, and (2) progressive metabolic and cytotoxic maturation driven by cytosolic calcium signaling, mitochondrial energy generation (MT-ATP8), and heightened granulysin (GNLY) expression. In contrast, Non-responders maintain a blunted Type-I IFN signature and fail to execute coordinate metabolic upregulation.
