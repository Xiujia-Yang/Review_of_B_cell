#!/share/apps/R/4.1.3/bin/Rscript 

'''
	This script constructs single-cell trajectory for naive
	B cell subpopulations (as an example).
	
	To avoid the noise brought with by the imbalance in 
	the number of cells in different cell subsets, 200 
	cells were sampled from each cluster.
'''

# Load required packages
library(monocle)

###############################################
# Step 1: Create CellDataSet
###############################################
set.seed(1)  # specify the seed for result reproduction
selected.bar <- c()

temp <- sample(colnames(bc.proj.clean)[bc.proj.clean$CellType=="Naive B"], 
               min(200, sum(bc.proj.clean$CellType=="Naive B")))
selected.bar <- c(selected.bar, temp)
temp <- sample(colnames(bc.proj.clean)[bc.proj.clean$CellType=="FOShi naive B"], 
               min(200, sum(bc.proj.clean$CellType=="FOShi naive B")))
selected.bar <- c(selected.bar, temp)
temp <- sample(colnames(bc.proj.clean)[bc.proj.clean$CellType=="HSPA1Ahi naive B"], 
               min(200, sum(bc.proj.clean$CellType=="HSPA1Ahi naive B")))
selected.bar <- c(selected.bar, temp)

bc.proj.clean.naive.select <- bc.proj.clean[, colnames(bc.proj.clean) %in% selected.bar]
#Extract data, phenotype data, and feature data from the SeuratObject
# refer to https://github.com/cole-trapnell-lab/monocle-release/issues/262
data <- as(as.matrix(bc.proj.clean.naive.select@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = bc.proj.clean.naive.select@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

###############################################
# Step 2: Estimate size factors and dispersions
###############################################
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

###############################################
# Step 3: Filtering low-quality cells
###############################################
# This step is not needed here, since both cells
# and genes have been subjected to quality controls.


###############################################
# Step 4: Classifying and Counting Cells
###############################################
# This step is not needed here, since we have their
# identifies beforehand.

###############################################
# Step 5: Constructing Single Cell Trajectories
###############################################
# Trajectory step 1: choose genes that define a
# cell's progress
expressed_genes <- rownames(data)[rowSums(data>0) >= 3]
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~CellType")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)

# Trajectory step 2: reduce data dimensionality
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                               method = 'DDRTree')

# Trajectory step 3: order cells along the trajectory
#monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, reverse=TRUE)  # for naive cells

plot_cell_trajectory(monocle_cds, color_by = "CellType", 
                     show_branch_points = F, show_state_number=F,
                     ) +
  scale_color_manual(breaks = c("Naive B", "FOShi naive B", "HSPA1Ahi naive B"), 
                     values=c("#14c38e", "#99f3bd", "#e3fcbf"))