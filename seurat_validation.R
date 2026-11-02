library(Seurat)

getwd()

data_dir <- "Data/Day0_NR1"

dir.exists(data_dir)

list.files(data_dir)

library(Matrix) #reading the 3 files separately as Read10X needs compressed files (.gz)

counts <- readMM(file.path(data_dir, "matrix.mtx"))

class(counts)

dim(counts)

features <- read.delim(
  file.path(data_dir, "features.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE
)

barcodes <- read.delim(
  file.path(data_dir, "barcodes.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE
)

dim(features)
dim(barcodes)

head(features)
head(barcodes)

#Assign gene names and cell barcodes
rownames(counts) <- features$V2
colnames(counts) <- barcodes$V1

sum(duplicated(rownames(counts))) #to check duplicate gene names

rownames(counts) <- make.unique(rownames(counts))

sum(duplicated(rownames(counts))) #to check duplicate gene names

#Create the Seurat Object
seu <- CreateSeuratObject(
  counts = counts,
  project = "Influenza"
)

seu

#Normalization with same command like python
seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

seu

#Highly Variable Genes
seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000
)

length(VariableFeatures(seu))
head(VariableFeatures(seu))

#Scaling the data before PCA
seu <- ScaleData(
  seu,
  features = VariableFeatures(seu)
)

seu

#PCA
seu <- RunPCA(
  seu,
  features = VariableFeatures(seu),
  npcs = 50
)

seu[["pca"]]

dim(Embeddings(seu,"pca"))

#Inspect PCA variance
ElbowPlot(seu, ndims = 50)

seu <- FindNeighbors(
  seu,
  dims = 1:15,
  k.param = 15
)

names(seu@graphs)

#UMAP
seu <- RunUMAP(
  seu,
  dims = 1:15,
  n.neighbors = 15,
  seed.use = 42
)

DimPlot(
  seu,
  reduction = "umap"
)

###End of pilot run only on one timepoint and patient


###Addition of rest of the data and main analysis

samples <- c(
  "Day0_NR1", "Day0_NR2", "Day0_NR3",
  "Day0_R1",  "Day0_R2",  "Day0_R3",
  "Day7_NR1", "Day7_NR2", "Day7_NR3",
  "Day7_R1",  "Day7_R2",  "Day7_R3"
)

seu_list <- list()

#Load the 12 raw count matrices
for (sample in samples) {
  
  data_dir <- file.path("Data", sample)
  
  counts <- readMM(file.path(data_dir, "matrix.mtx"))
  
  features <- read.delim(
    file.path(data_dir, "features.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  barcodes <- read.delim(
    file.path(data_dir, "barcodes.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  rownames(counts) <- make.unique(features$V2)
  colnames(counts) <- barcodes$V1
  
  seu_list[[sample]] <- CreateSeuratObject(
    counts = counts,
    project = "Influenza"
  )
  
  print(paste(sample, "loaded:", ncol(seu_list[[sample]]), "cells"))
}

#Calculate QC metrics
for (sample in names(seu_list)) {
  
  seu_list[[sample]][["percent.mt"]] <- PercentageFeatureSet(
    seu_list[[sample]],
    pattern = "^MT-"
  )
}

head(seu_list[["Day0_NR1"]]@meta.data)

#Apply same QC filter as python to all 12 samples
for (sample in names(seu_list)) {
  
  seu_list[[sample]] <- subset(
    seu_list[[sample]],
    subset = nFeature_RNA >= 500 & percent.mt <= 20
  )
}

for (sample in names(seu_list)) {
  print(paste(sample, ":", ncol(seu_list[[sample]]), "cells"))
}

#Add biological metadata
for (sample in names(seu_list)) {
  
  seu_list[[sample]]$Sample <- sample
  
  seu_list[[sample]]$Time <- ifelse(
    grepl("^Day0", sample),
    "Day0",
    "Day7"
  )
  
  seu_list[[sample]]$Donor <- sub(
    "^(Day0|Day7)_",
    "",
    sample
  )
  
  seu_list[[sample]]$Response <- ifelse(
    grepl("_R", sample),
    "Responder",
    "Non-responder"
  )
}

head(seu_list[["Day7_R2"]]@meta.data)

#Merge 12 samples
seu <- merge(
  x = seu_list[[1]],
  y = seu_list[2:length(seu_list)],
  add.cell.ids = names(seu_list),
  project = "Influenza"
)


#Normalization
seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

print(seu)

#Highly Variable Genes
seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000
)

print(length(VariableFeatures(seu)))

#Scaling
seu <- ScaleData(
  seu,
  features = VariableFeatures(seu)
)

seu <- RunPCA(
  seu,
  features = VariableFeatures(seu),
  npcs = 50
)

print(seu[["pca"]])

ElbowPlot(
  seu,
  ndims = 50
)

#Using 30 PC to directly compare with scanpy
seu <- FindNeighbors(
  seu,
  dims = 1:30,
  k.param = 15
)
seu <- RunUMAP(
  seu,
  dims = 1:30,
  n.neighbors = 15,
  seed.use = 42
)
DimPlot(
  seu,
  reduction = "umap"
)

#Proceeding with 30 PC as difference with 15PC is minimal


seu <- FindClusters(
  seu,
  resolution = 0.2
)

#Got 12 clusters with 0.2 resoulution. Will proceed to compare with 10 clusters of scanpy with 0.5 resolution

saveRDS(seu, "influenza_seurat_clustered.rds")


#Loading the saved data
seu <- readRDS("influenza_seurat_clustered.rds")

seu <- JoinLayers(seu)

#Cluster Annotation and Marker Validation
marker_genes <- c(
  "CD3D", "CD3E", "TRAC",
  "IL7R", "LTB", "CCR7",
  "CD8A", "CD8B",
  "NKG7", "GNLY", "KLRD1",
  "MS4A1", "CD79A", "CD37",
  "LYZ", "S100A8", "S100A9",
  "FCGR3A", "MS4A7",
  "FCER1A", "CST3"
)

library(ggplot2)

DotPlot(
  seu,
  features = marker_genes
) +
  RotatedAxis() +
  scale_size(
    range = c(1, 8),
    name = "% Expressing"
  )


cell_type_map <- c(
  "0"  = "Naive/Central Memory CD4 T cells",
  "1"  = "Resting CD4 T cells",
  "2"  = "NK cells",
  "3"  = "NK/T-like cytotoxic cells",
  "4"  = "NK cells",
  "5"  = "Cytotoxic CD8 T cells",
  "6"  = "B cells",
  "7"  = "NK/T-like cytotoxic cells",
  "8"  = "B cells",
  "9"  = "Cytotoxic CD8 T cells",
  "10" = "Naive/central memory mixed T cells",
  "11" = "B cells"
)

colnames(seu@meta.data)
table(seu$seurat_clusters)

cell_type_map
seu$cell_type <- unname(
  cell_type_map[as.character(seu$seurat_clusters)]
)

table(seu$cell_type)


library(ggplot2)

p <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "cell_type"
) +
  NoLegend()

LabelClusters(
  plot = p,
  id = "cell_type",
  repel = TRUE,
  size = 4,
  color = "black",
  fontface = "bold",
  box.padding = 1.5,
  point.padding = 0.5
)



