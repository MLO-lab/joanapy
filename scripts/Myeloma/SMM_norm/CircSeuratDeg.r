setwd("pathto/rawData/single_cell/SMM_norm/")
library(Seurat)

matSMMCirc=read.csv("/home/samerifar/joana/Myeloma/counts/counts_SMM_circPC#.csv",header=T)
rownames(matSMMCirc)=matSMMCirc[,1]
matSMMCirc=matSMMCirc[,-1]
matSMMCirc=as.matrix(matSMMCirc)

matNormCirc=read.table("./countsNormCirc.txt")
dim(matSMMCirc)
dim(matNormCirc)
rownames(matSMMCirc)[1:10]
rownames(matNormCirc)[1:10]
print(identical(rownames(matSMMCirc),rownames(matNormCirc)))





# Create Seurat objects for each condition
seurat_condition1 <- CreateSeuratObject(counts = matSMMCirc, project = "Condition1")
seurat_condition2 <- CreateSeuratObject(counts = matNormCirc, project = "Condition2")

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

write.table(deg_results, file = "./SMM_normDegCirc.txt", col.names=T,row.names=T)

                           
