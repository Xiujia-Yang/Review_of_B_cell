#!/share/apps/R/4.1.3/bin/Rscript 

'''
	This script implements the integration all of scRNA-seq 
	data included in this study, which are composed of in-house 
	enriched B cells (BM+PB) and BMMCs and external BM early
	B cells, PBMC and tonsil immune cells (GC).
'''


library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)


###########################################################
#  Preparation of in-house scRNA-seq data of 
#  Enriched B cells (BM+PB)
###########################################################
# Load the single B cell transcriptomic dataset
enriched.b.data <- Read10X(data.dir = "/path/to/the/merged/enriched/b/data")
# Initialize the Seurat object with the raw (non-normalized data).
enriched.b <- CreateSeuratObject(counts = enriched.b.data, project = "Enriched B", min.cells = 0, min.features = 0)
# Assign `donor` and `tissue` label
enriched.b$donor <- str_split_fixed(colnames(enriched.b), '-', n=3)[,2]
tissues <- c("BM", "PB")
names(tissues) <- c("BB", "PB")
enriched.b$tissue <- tissues[str_sub(str_split_fixed(colnames(enriched.b), '-', 4)[,3], 1, 2)]
# QC and selecting cells for further analysis
enriched.b[["percent.mt"]] <- PercentageFeatureSet(enriched.b, pattern = "^MT-")
# Remove BCR and TCR VDJ genes
subclasses <- c('IGHM', 'IGHD', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHE')
enriched.b <- enriched.b[(!str_detect(rownames(enriched.b),"^IGH") | rownames(enriched.b)%in%subclasses)
                                           & (!str_detect(rownames(enriched.b),"^IGK") | str_detect(rownames(enriched.b),"^IGKC"))
                                           & (!str_detect(rownames(enriched.b),"^IGL") | str_detect(rownames(enriched.b),"^IGLC"))
                                           & (!str_detect(rownames(enriched.b),"^TRB") | str_detect(rownames(enriched.b),"^TRBC"))
                                           & (!str_detect(rownames(enriched.b),"^TRA") | str_detect(rownames(enriched.b),"^TRAC"))
                                           & (!str_detect(rownames(enriched.b),"^TRD") | str_detect(rownames(enriched.b),"^TRDC"))
                                           & (!str_detect(rownames(enriched.b),"^TRG") | str_detect(rownames(enriched.b),"^TRGC"))]
# Cell and gene QC
selected_c <- WhichCells(enriched.b, expression = nFeature_RNA >= 200 & percent.mt < 10)
selected_f <- rownames(enriched.b)[Matrix::rowSums(enriched.b) >= 3]
enriched.b <- subset(enriched.b, features = selected_f, cells = selected_c)
# Remove doublets 
doublet_dir <- '/path/to/predicted/doublet/result'
doublet_files <- list.files(doublet_dir, '^R-[A-Z]{3,5}-(BB|PB).*_df_prediction.csv')
for(i in 1:length(doublet_files)){
  df_temp <- read.table(paste0(doublet_dir, doublet_files[i]), sep=',', row.names=1)
  sample_name <- str_split(doublet_files[i], '_')[[1]][1]
  barcode <- paste0(sample_name, '-', rownames(df_temp))
  rownames(df_temp) <- barcode
  if(i==1){
    df_merge <- df_temp
  }else{
    df_merge <- rbind(df_merge, df_temp)
  }
}
singlet_bar <- rownames(df_merge)[df_merge$V2=="Singlet"]
enriched.b <- enriched.b[,colnames(enriched.b) %in% singlet_bar]
# Split the sample
enriched.b$sample <- paste0(str_split_fixed(colnames(enriched.b), '-', n=3)[,2], '-', str_split_fixed(colnames(enriched.b), '-', n=4)[,3])
enriched.b.list <- SplitObject(enriched.b, split.by = "sample")


###########################################################
#  Preparation of in-house scRNA-seq data of BMMC
###########################################################
# Load the single B cell transcriptomic dataset
bmmc.data <- Read10X(data.dir = "/path/to/the/merged/BMMC/data")
# Initialize the Seurat object with the raw (non-normalized data).
bmmc <- CreateSeuratObject(counts = bmmc.data, project = "BMMC", min.cells = 0, min.features = 0)
bmmc$donor <- str_split_fixed(colnames(bmmc), '-', 3)[,2]
# QC and selecting cells for further analysis
bmmc[["percent.mt"]] <- PercentageFeatureSet(bmmc, pattern = "^MT-")
# Remove BCR and TCR VDJ genes
subclasses <- c('IGHM', 'IGHD', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHE')
bmmc <- bmmc[(!str_detect(rownames(bmmc),"^IGH") | rownames(bmmc) %in% subclasses)
             & (!str_detect(rownames(bmmc),"^IGK") | str_detect(rownames(bmmc),"^IGKC"))
             & (!str_detect(rownames(bmmc),"^IGL") | str_detect(rownames(bmmc),"^IGLC"))
             & (!str_detect(rownames(bmmc),"^TRB") | str_detect(rownames(bmmc),"^TRBC"))
             & (!str_detect(rownames(bmmc),"^TRA") | str_detect(rownames(bmmc),"^TRAC"))
             & (!str_detect(rownames(bmmc),"^TRD") | str_detect(rownames(bmmc),"^TRDC"))
             & (!str_detect(rownames(bmmc),"^TRG") | str_detect(rownames(bmmc),"^TRGC"))]

selected_c <- WhichCells(bmmc, expression = nFeature_RNA >= 200 & percent.mt < 10)
selected_f <- rownames(bmmc)[Matrix::rowSums(bmmc) >= 3]
bmmc <- subset(bmmc, features = selected_f, cells = selected_c)
# Remove doublets 
doublet_dir <- '/path/to/predicted/doublet/result'
doublet_files <- list.files(doublet_dir, '^R-[A-Z]{3,5}-BM.*_df_prediction.csv')
for(i in 1:length(doublet_files)){
  df_temp <- read.table(paste0(doublet_dir, doublet_files[i]), sep=',', row.names=1)
  sample_name <- str_split(doublet_files[i], '_')[[1]][1]
  barcode <- paste0(sample_name, '-', rownames(df_temp))
  rownames(df_temp) <- barcode
  if(i==1){
    df_merge <- df_temp
  }else{
    df_merge <- rbind(df_merge, df_temp)
  }
}
singlet_bar <- rownames(df_merge)[df_merge$V2=="Singlet"]
bmmc <- bmmc[,colnames(bmmc) %in% singlet_bar]
# Split the sample
sample <- paste0(str_split_fixed(colnames(bmmc), '-', n=3)[,2], '.', str_split_fixed(colnames(bmmc), '-', n=4)[,3])
bmmc$sample <- sample
bmmc.list <- SplitObject(bmmc, split.by = "sample")


###########################################################
#  Preparation of external scRNA-seq data of early B cells
###########################################################
# Load the single B cell transcriptomic dataset
hd1.cd34.data <- Read10X(data.dir = "/path/to/10.1038-s41556-021-00814-7/HRR178260/filtered_feature_bc_matrix")
hd2.cd34.data <- Read10X(data.dir = "/path/to/10.1038-s41556-021-00814-7/HRR178261/filtered_feature_bc_matrix")
hd1.cd19.data <- Read10X(data.dir = "/path/to/10.1038-s41556-021-00814-7/HRR178258/filtered_feature_bc_matrix")
hd2.cd19.data <- Read10X(data.dir = "/path/to/10.1038-s41556-021-00814-7/HRR178259/filtered_feature_bc_matrix")
hd1.cd34 <- CreateSeuratObject(counts = hd1.cd34.data, project = "BM_CD34", min.cells = 0, min.features = 0)
hd2.cd34 <- CreateSeuratObject(counts = hd2.cd34.data, project = "BM_CD34", min.cells = 0, min.features = 0)
hd1.cd19 <- CreateSeuratObject(counts = hd1.cd19.data, project = "BM_CD19", min.cells = 0, min.features = 0)
hd2.cd19 <- CreateSeuratObject(counts = hd2.cd19.data, project = "BM_CD19", min.cells = 0, min.features = 0)
# Merge individual samples
pro_hd <- merge(hd1.cd34, y = c(hd2.cd34, hd1.cd19, hd2.cd19), 
                add.cell.ids = c("hd1.cd34", "hd2.cd34", "hd1.cd19", "hd2.cd19"), project = "Progenitor")
# QC and selecting cells for further analysis
pro_hd[["percent.mt"]] <- PercentageFeatureSet(pro_hd, pattern = "^MT-")
# Remove BCR and TCR VDJ genes
subclasses <- c('IGHM', 'IGHD', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHE')
pro_hd <- pro_hd[(!str_detect(rownames(pro_hd),"^IGH") | rownames(pro_hd) %in% subclasses)
                   & (!str_detect(rownames(pro_hd),"^IGK") | str_detect(rownames(pro_hd),"^IGKC"))
                   & (!str_detect(rownames(pro_hd),"^IGL") | str_detect(rownames(pro_hd),"^IGLC"))
                   & (!str_detect(rownames(pro_hd),"^TRB") | str_detect(rownames(pro_hd),"^TRBC"))
                   & (!str_detect(rownames(pro_hd),"^TRA") | str_detect(rownames(pro_hd),"^TRAC"))
                   & (!str_detect(rownames(pro_hd),"^TRD") | str_detect(rownames(pro_hd),"^TRDC"))
                   & (!str_detect(rownames(pro_hd),"^TRG") | str_detect(rownames(pro_hd),"^TRGC"))]
selected_c <- WhichCells(pro_hd, expression = nFeature_RNA >= 200 & percent.mt < 10)
selected_f <- rownames(pro_hd)[Matrix::rowSums(pro_hd) >= 3]
pro_hd <- subset(pro_hd, features = selected_f, cells = selected_c)
# Remove doublets 
doublet_dir <- '/path/to/predicted/doublet/result'
doublet_files <- list.files(doublet_dir, '*_df_prediction.csv')
sample.ids <- c("hd1.cd34", "hd2.cd34", "hd1.cd19", "hd2.cd19")
names(sample.ids) <- c("HRR178260", "HRR178261", "HRR178258", "HRR178259")
for(i in 1:4){
  df_temp <- read.table(paste0(doublet_dir, doublet_files[i]), sep=',', row.names=1)
  sample_name <- str_split(doublet_files[i], '_')[[1]][1]
  sample_name <- sample.ids[sample_name]
  barcode <- paste0(sample_name, '_', rownames(df_temp))
  rownames(df_temp) <- barcode
  if(i==1){
    df_merge <- df_temp
  }else{
    df_merge <- rbind(df_merge, df_temp)
  }
}
singlet_bar <- rownames(df_merge)[df_merge$V2=="Singlet"]
pro_hd <- pro_hd[,colnames(pro_hd) %in% singlet_bar]

# Split the sample
sample <- str_split_fixed(colnames(pro_hd), '_', n=2)[,1]
pro_hd$sample <- sample
pro_hd.list <- SplitObject(pro_hd, split.by = "sample")


###########################################################
#  Preparation of external scRNA-seq data of tonsil B cells
###########################################################
# Import the raw data
setwd("/path/to/10.1126-sciimmunol.abe629/E-MTAB-9005")
BCP002_Total <- Read10X(data.dir = "BCP002_Total_3GEX_filtered_feature_bc_matrix")
BCP003_Total <- Read10X(data.dir = "BCP003_Total_5GEX_filtered_feature_bc_matrix")
BCP004_Total <- Read10X(data.dir = "BCP004_Total_5GEX_filtered_feature_bc_matrix")
BCP005_Total <- Read10X(data.dir = "BCP005_Total_5GEX_filtered_feature_bc_matrix")
BCP006_Total <- Read10X(data.dir = "BCP006_Total_5GEX_filtered_feature_bc_matrix")
BCP008_Total <- Read10X(data.dir = "BCP008_Total_5GEX_filtered_feature_bc_matrix")
BCP009_Total <- Read10X(data.dir = "BCP009_Total_5GEX_filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
BCP002_total <- CreateSeuratObject(counts = BCP002_Total, project = "GC", min.cells = 0, min.features = 0)
BCP003_total <- CreateSeuratObject(counts = BCP003_Total, project = "GC", min.cells = 0, min.features = 0)
BCP004_total <- CreateSeuratObject(counts = BCP004_Total, project = "GC", min.cells = 0, min.features = 0)
BCP005_total <- CreateSeuratObject(counts = BCP005_Total, project = "GC", min.cells = 0, min.features = 0)
BCP006_total <- CreateSeuratObject(counts = BCP006_Total, project = "GC", min.cells = 0, min.features = 0)
BCP008_total <- CreateSeuratObject(counts = BCP008_Total, project = "GC", min.cells = 0, min.features = 0)
BCP009_total <- CreateSeuratObject(counts = BCP009_Total, project = "GC", min.cells = 0, min.features = 0)
# Merge the all objects and then perform the downstream standard analyses
gc.sample.list <- c(BCP003_total, BCP004_total, BCP005_total,
                     BCP006_total, BCP008_total, BCP009_total)
gc.sample.ids <- c("BCP002_total", "BCP003_total", "BCP004_total", "BCP005_total",
                    "BCP006_total", "BCP008_total", "BCP009_total")
gc.combined <- merge(BCP002_total, y = gc.sample.list, 
                      add.cell.ids = gc.sample.ids, project = "GC project")
# Remove predicted doublets and TCR and BCR-related genes
gc.combined[["percent.mt"]] <- PercentageFeatureSet(gc.combined, pattern = "^MT-")
subclasses <- c('IGHM', 'IGHD', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHE')
gc.combined <- gc.combined[(!str_detect(rownames(gc.combined),"^IGH") | rownames(gc.combined)%in%subclasses)
                             & (!str_detect(rownames(gc.combined),"^IGK") | str_detect(rownames(gc.combined),"^IGKC"))
                             & (!str_detect(rownames(gc.combined),"^IGL") | str_detect(rownames(gc.combined),"^IGLC"))
                             & (!str_detect(rownames(gc.combined),"^TRB") | str_detect(rownames(gc.combined),"^TRBC"))
                             & (!str_detect(rownames(gc.combined),"^TRA") | str_detect(rownames(gc.combined),"^TRAC"))
                             & (!str_detect(rownames(gc.combined),"^TRD") | str_detect(rownames(gc.combined),"^TRDC"))
                             & (!str_detect(rownames(gc.combined),"^TRG") | str_detect(rownames(gc.combined),"^TRGC"))]
selected_c <- WhichCells(gc.combined, expression = nFeature_RNA >= 200 & percent.mt < 10)
selected_f <- rownames(gc.combined)[Matrix::rowSums(gc.combined) >= 3]
gc.combined <- subset(gc.combined, features = selected_f, cells = selected_c)
# Remove doublets 
doublet_dir <- '/path/to/predicted/doublet/result'
doublet_files <- list.files(doublet_dir, '*Total_df_prediction.csv')
names(gc.sample.ids) <- c("BCP002_Total", "BCP003_Total", "BCP004_Total", "BCP005_Total",
                           "BCP006_Total", "BCP008_Total", "BCP009_Total")
for(i in 1:length(doublet_files)){
  df_temp <- read.table(paste0(doublet_dir, doublet_files[i]), sep=',', row.names=1)
  sample_name <- paste0(str_split(doublet_files[i], '_')[[1]][1], "_", str_split(doublet_files[i], '_')[[1]][2])
  sample_name <- gc.sample.ids[sample_name]
  barcode <- paste0(sample_name, '_', rownames(df_temp))
  rownames(df_temp) <- barcode
  if(i==1){
    df_merge <- df_temp
  }else{
    df_merge <- rbind(df_merge, df_temp)
  }
}
singlet_bar <- rownames(df_merge)[df_merge$V2=="Singlet"]
gc.combined <- gc.combined[,colnames(gc.combined) %in% singlet_bar]
# Split the sample
sample <- str_split_fixed(colnames(gc.combined), '_', n=2)[,1]
gc.combined$sample <- sample
gc.combined.list <- SplitObject(gc.combined, split.by = "sample")


###########################################################
#  Preparation of external scRNA-seq data of PBMC dataset
###########################################################
pbmc <- readRDS("/path/to/10.1016-j.immuni.2020.07.009/Final_nCoV_0716_upload.RDS")
pbmc.3ctrl <- pbmc[,pbmc$batch %in% c("Ctrl-1", "Ctrl-2", "Ctrl-3")]
DefaultAssay(pbmc.3ctrl) <- "RNA"


###########################################################
#      Merge individual datasets into a single object
###########################################################
data.list <- append(append(append(append(pro_hd.list, bmmc.list), list(pbmc=pbmc.3ctrl)), gc.combined.list), enriched.b.list)


###########################################################
#      Perform dataset integration
###########################################################

# Normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# ScaleDataelect features that are repeatedly variable across datasets 
# for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Identify anchors and use these anchors to integrate the two datasets
# we set the BMMC data as reference (`DENG.BM2 has the maximum number of cells`)
data.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, reduction = "rpca", reference=8)

# This command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = data.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.8)

# Then cell lineages were mannually annotated with canonical marker genes