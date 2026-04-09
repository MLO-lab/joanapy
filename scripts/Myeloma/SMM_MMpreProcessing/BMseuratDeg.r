setwd("pathto/rawData/single_cell/SMM_MM/")
library(Seurat)
#library(dplyr)

print(getwd)

matSMMBM=read.csv("./counts_SMM_BMPC#.csv",header=T)
rownames(matSMMBM)=matSMMBM[,1]
matSMMBM=matSMMBM[,-1]
matSMMBM=as.matrix(matSMMBM)



matMMBM=read.csv("./counts_MM_BMPC#.csv",header=T)
rownames(matMMBM)=matMMBM[,1]
matMMBM=matMMBM[,-1]
matMMBM=as.matrix(matMMBM)


dim(matSMMBM)
dim(matMMBM)
rownames(matSMMBM)[1:10]
rownames(matMMBM)[1:10]
print(identical(rownames(matSMMBM),rownames(matMMBM)))





# Create Seurat objects for each condition
seurat_condition1 <- CreateSeuratObject(counts = matSMMBM, project = "Condition1")
seurat_condition2 <- CreateSeuratObject(counts = matMMBM, project = "Condition2")

# Add metadata to indicate the condition
seurat_condition1$condition <- "Condition1"
seurat_condition2$condition <- "Condition2"

# Merge the two Seurat objects
combined_seurat <- merge(seurat_condition1, y = seurat_condition2, add.cell.ids = c("Condition1", "Condition2"))

# Normalize the data
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)

combined_seurat <- JoinLayers(combined_seurat)

# Identify clusters
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Set the identities to the condition labels
Idents(combined_seurat) <- combined_seurat$condition

# Perform differential expression analysis using Wilcoxon test with loose thresholds
deg_results <- FindMarkers(combined_seurat, ident.1 = "Condition1", ident.2 = "Condition2", 
                           min.pct = 0, logfc.threshold = 0)

write.table(deg_results, file = "./SMM_MMDegBM.txt", col.names=T,row.names=T)

                           
