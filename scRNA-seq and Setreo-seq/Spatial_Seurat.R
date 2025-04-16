library(Seurat)
library(ggplot2)

ifnb <- readRDS("MmuRbMerge.rds")
# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "batch")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

seurat_spatialObj <- immune.combined

seurat_spatialObj <- ScaleData(seurat_spatialObj, verbose = FALSE)

####----2. cluster----####

seurat_spatialObj <- RunPCA(seurat_spatialObj,verbose = T)
seurat_spatialObj <- FindNeighbors(seurat_spatialObj, reduction = "pca", dims = 1:20)
seurat_spatialObj <- RunUMAP(seurat_spatialObj, reduction = "pca", dims = 1:20)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = T,resolution = 0.5)

pdf("MmuRb_merge_UMAP.pdf", width = 5, height = 5)
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE)
dev.off()

pdf("MmuRb_merge_UMAP_split.pdf", width = 15, height = 5)
DimPlot(seurat_spatialObj, reduction = "umap", label = T, split.by = "batch")
dev.off()

pdf("MmuRb_merge_UMAP_group.pdf", width = 5, height = 5)
DimPlot(seurat_spatialObj, reduction = "umap", label = T, group.by = "batch")
dev.off()

#pdf("merge_SpatialDimPlot.pdf", width = 15, height = 5)
#SpatialDimPlot(seurat_spatialObj, label = F, stroke = 0, pt.size.factor = 1.2) + scale_y_reverse()
#dev.off()

#saveRDS(seurat_spatialObj, "MmuRb_merge_cluster.rds")

