#!/share/apps/R/4.1.3/bin/Rscript 

'''
	This script extracts B lineage cells, integrates
	different datasets, re-clusters B cells and annotate
	B cell subpopulations.
'''

# load required packages
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)

###########################################################
#              Extract B lineage cells
###########################################################
bc <- data.combined[,data.combined$cell_group=="B lineage"]


###########################################################
#              Integrate different datasets
###########################################################
# Split Seurat object to generate dataset list for 
# downstream integration
bc.proj.list <- SplitObject(bc, split.by = "Project")

# Normalize and identify variable features for each dataset
# independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across 
# datasets for integration run PCA on each dataset using 
# these features
features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"))
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Identify anchors and use these anchors to integrate the 
# two datasets. We set the BMMC data as reference 
# (`DENG.BM2 has the maximum number of cells`).
data.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, reduction = "rpca", reference=28)

# This command creates an 'integrated' data assay
bc.proj <- IntegrateData(anchorset = data.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"

###########################################################
#              Re-cluster B cells
###########################################################
# Run the standard workflow for visualization and clustering
bc.proj <- ScaleData(bc.proj, verbose = FALSE)
bc.proj <- RunPCA(bc.proj, npcs = 30, verbose = FALSE)
bc.proj <- RunUMAP(bc.proj, reduction = "pca", dims = 1:30)
bc.proj <- FindNeighbors(bc.proj, reduction = "pca", dims = 1:30)
bc.proj <- FindClusters(bc.proj, resolution = 0.8)


###########################################################
#              Annotate B cell subpopulations
###########################################################
Idents(bc.proj) <- "seurat_clusters_0.8"
new.cluster.ids <- c("Naive B", "Naive B", "IgM+ memory B", "Immature B", "Naive B",
                     "Pre B", "LZ GC B", "FOShi naive B", "DZ GC B", "IgM+ memory B",
                     "Classical memory B", "S100A8hi immature B", "FTLhi immature B", "CD27-IgM+IgD+ memory B", "PC",
                     "PC", "Erythrocyte", "Cycling pre B", "HSPA1Ahi naive B", "Erythrocyte",
                     "Pre-pro B", "Mixed", "LZ GC B", "Cycling pro B", "Pro B",
                     "DZ GC B")
names(new.cluster.ids) <- levels(bc.proj$seurat_clusters)
bc.proj@meta.data["CellType"] <- new.cluster.ids[bc.proj$seurat_clusters]
bc.proj@meta.data["CellType"] <- factor(bc.proj$CellType, 
                levels = c("Pre-pro B", "Cycling pro B", "Pro B", 
                               "Cycling pre B", "Pre B", 
                               "Immature B", "FTLhi immature B", "S100A8hi immature B",
                               "Naive B", "FOShi naive B", "HSPA1Ahi naive B",
                               "LZ GC B", "DZ GC B", 
                               "Classical memory B", "IgM+ memory B", "CD27-IgM+IgD+ memory B", "PC",
                               "Erythrocyte", "Mixed"))
cell.list <- c("Erythrocyte", "Mixed")
bc.proj.clean <- bc.proj[,! bc.proj$CellType %in% cell.list]
Idents(bc.proj.clean) <- "CellType"
bc.proj.clean@meta.data["CellType"] <- factor(as.vector(bc.proj.clean$CellType), 
                  levels = c("Pre-pro B", "Cycling pro B", "Pro B", 
                             "Cycling pre B", "Pre B", 
                             "Immature B", "FTLhi immature B", "S100A8hi immature B",
                             "Naive B", "FOShi naive B", "HSPA1Ahi naive B",
                             "LZ GC B", "DZ GC B", 
                             "Classical memory B", "IgM+ memory B", "CD27-IgM+IgD+ memory B", 
                             "PC"))
							 

###########################################################
#        Distinguish plasmablast from plasma cells
###########################################################							 
# Extract PC for re-clustering
pc.subset <- bc.proj.clean[,bc.proj.clean$CellType=="PC"]
# Re-clustering
DefaultAssay(pc.subset) <- "integrated"
pc.subset <- ScaleData(pc.subset, verbose = FALSE)
pc.subset <- RunPCA(pc.subset, npcs = 30, verbose = FALSE)
pc.subset <- RunUMAP(pc.subset, reduction = "pca", dims = 1:30)
pc.subset <- FindNeighbors(pc.subset, reduction = "pca", dims = 1:30)
pc.subset <- FindClusters(pc.subset, resolution = 0.5)
# Annotate the cell clusters
new.cluster.ids <- c("HighNumber", "HighNumber", "LowNumber", "HighNumber", "HighNumber", 
                     "LowNumber", "LowNumber", "LowNumber", "LowNumber", "LowNumber", 
                     "LowNumber", "HighNumber", "LowNumber", "HighNumber", "LowNumber", 
                     "LowNumber", "HighNumber")
names(new.cluster.ids) <- levels(pc.subset$seurat_clusters)
pc.subset@meta.data["CellTypePC"] <- new.cluster.ids[pc.subset$seurat_clusters]
pc.subset@meta.data["CellTypePC"] <- factor(pc.subset$CellTypePC, 
                                         levels = c("HighNumber", "LowNumber"))
DimPlot(pc.subset, reduction = "umap", group.by="CellTypePC", label = FALSE)

# Reassign B cell subpopulation label for `bc.proj.clean` object
pb.barcodes <- colnames(pc.subset)[pc.subset$CellTypePC=="LowNumber"]
pb.anno <- rep("PB", length(pb.barcodes))
names(pb.anno) <- pb.barcodes
other.barcodes <- colnames(bc.proj.clean)[!colnames(bc.proj.clean) %in% pb.barcodes]
other.anno <- as.vector(bc.proj.clean$CellType)[!colnames(bc.proj.clean) %in% pb.barcodes]
names(other.anno) <- other.barcodes
comb.anno <- c(pb.anno, other.anno)

bc.proj.clean$CellTypeGamma <- factor(as.vector(comb.anno[colnames(bc.proj.clean)]), 
			  levels=c("Pre-pro B", "Cycling pro B", "Pro B", "Cycling pre B", "Pre B",
						"Immature B", "FTLhi immature B", "S100A8hi immature B",
						"Naive B", "FOShi naive B", "HSPA1Ahi naive B",
						"LZ GC B", "DZ GC B",
						"Classical memory B", "IgM+ memory B", "CD27-IgM+IgD+ memory B", "PB", "PC"))