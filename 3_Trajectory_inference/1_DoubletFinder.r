#!/share/apps/R/4.1.3/bin/Rscript 

# Passing the command line parameters
args <- commandArgs()
input_dir <- args[6]  # specifying the directory containing expression matrix
output_fl <- args[7]  # specifying the file storing the predicted result

# Load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)

# Set the working directory
setwd("/path/to/working/directory")

# Read the expression matrix
tumor.data <- Read10X(data.dir = input_dir)

# Initialize the Seurat object with the raw (non-normalized data)
tumor <- CreateSeuratObject(counts = tumor.data, project = "tumor", min.cells = 3, min.features = 200)

# Calculate the percentage of mitochondria-derived genes
tumor[["percent.mt"]] <- PercentageFeatureSet(tumor, pattern = "^MT-")


################################################################
### The following codes are extracted from the authorized 
### website, 'https://github.com/chris-mcginnis-ucsf/DoubletFinder'.
################################################################

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
tumor <- NormalizeData(tumor)
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000)
tumor <- ScaleData(tumor)
tumor <- RunPCA(tumor)
tumor <- RunUMAP(tumor, dims = 1:10)
tumor <- FindNeighbors(tumor, dims = 1:10)  # tumor and tumor_SD
tumor <- FindClusters(tumor, resolution = 0.5)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_tumor <- paramSweep_v3(tumor, PCs = 1:10, sct = FALSE)
sweep.stats_tumor <- summarizeSweep(sweep.res.list_tumor, GT = FALSE)
bcmvn_tumor <- find.pK(sweep.stats_tumor)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(tumor$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round((nrow(tumor@meta.data)/10000)*0.08*nrow(tumor@meta.data))  ## Assuming a dataset-dependent doublet formation rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
tumor <- doubletFinder_v3(tumor, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
tumor <- doubletFinder_v3(tumor, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = names(tumor[[]][7]), sct = FALSE)

# Write the prediction result
write.table(tumor[[]][9], output_fl, col.names=FALSE, sep=',')