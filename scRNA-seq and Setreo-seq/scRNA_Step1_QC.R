#!/usr/bin/Rscript
args=commandArgs(T)

### Seurat clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(future)
library(SeuratObject)
plan()
plan("multiprocess", workers = as.numeric(args[1]))
plan()
options(future.globals.maxSize= 429496729600)

organ <- args[2] 
SampleID <- args[3] 
minUMIs <- as.numeric(args[4]) #1000
minGenes <- as.numeric(args[5]) #400
maxPercent.mt <- as.numeric(args[6]) #10
dim.usage <- as.numeric(args[7]) #20
res.usage <- as.numeric(args[8]) #0.5
doublets.percentage <- as.numeric(args[9]) #0.05

setwd(args[10])
count <- Read10X(SampleID)
colnames(count) <- paste(SampleID,colnames(count),sep=":")


###outdir
dir.create(args[11], recursive = TRUE)
setwd(args[11]) 
dir.create(SampleID, recursive = TRUE)
setwd(SampleID)

#Creat Seurat object
scRNA = CreateSeuratObject(counts = count, min.cells = 1, min.features=1)
cat("1.raw matrix:",dim(scRNA),"\n")
scRNA@meta.data$Batch <- SampleID
cat("Batch:\n")
table(scRNA@meta.data$Batch)


### "MT-"
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
saveRDS(scRNA,"preQC.rds")

### Plot raw QC plot
pdf(paste0(SampleID,"_raw_QCplot_Batch.pdf"),10,5)
p <- VlnPlot(scRNA, group.by="Batch", pt.size=0.01, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
print(p)
dev.off()

pdf(paste0(SampleID,"_raw_scatterplot_Batch.pdf"),12,5)
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", group.by="Batch", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", group.by="Batch", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


### QC
qc_stat <- rbind(summary(scRNA@meta.data$nCount_RNA),summary(scRNA@meta.data$nFeature_RNA),summary(scRNA@meta.data$percent.mt))
rownames(qc_stat) <- c("nUMIs","nGenes","percent.mt")
write.table(qc_stat, file = paste0(SampleID,"_raw_qc_stat.txt"), quote = F, sep = "\t",row.names = T,col.names=NA)

scRNA <- subset(scRNA, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & percent.mt < maxPercent.mt)

scRNA <- NormalizeData(object = scRNA, verbose = FALSE)
scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

pdf(paste0(SampleID,"_HVG.pdf"),6,5)
p <- LabelPoints(plot = VariableFeaturePlot(scRNA), points = (head(VariableFeatures(scRNA), 20)), repel = TRUE)
print(p)
dev.off()

scRNA <- ScaleData(scRNA) 
scRNA <- RunPCA(scRNA, verbose=F)
scRNA <- RunUMAP(scRNA, dims = 1:dim.usage)

### Define Find_doublet function
Find_doublet <- function(data){
	sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn <- find.pK(sweep.stats) ### output plot
	nExp_poi <- round(doublets.percentage*ncol(data))
	p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK)) ### pK Selection
	data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) ### output plot
	colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
	return(data)
}


pdf(paste0(SampleID,"_find_doublet.pdf"),7,7)
scRNA <- Find_doublet(scRNA)
dev.off()
write.table(scRNA@meta.data,paste0(SampleID,"_doublets_info.txt"),quote = F, sep = "\t",row.names = T,col.names=NA)

scRNA <- subset(scRNA,subset=doublet_info=="Singlet")

cat("2.FilterDoublet matrix:",dim(scRNA),"\n")

### Find cluster
scRNA <- FindNeighbors(scRNA, dims = 1:dim.usage)
scRNA <- FindClusters(scRNA, resolution = res.usage)

saveRDS(scRNA, file = paste0(SampleID,"_scRNA.rds"))

### Cluster_Batch
pdf(paste0(SampleID,"_umap_cluster.pdf"),8,7)
DimPlot(object = scRNA, reduction = "umap",pt.size = 0.1,label=T)+ggtitle(label = "umap_Cluster")
DimPlot(object = scRNA, reduction = "umap",group.by = "Batch",pt.size = 0.1,label=F)+ggtitle(label = "umap_Batch")
dev.off()

### umi_gene
pdf(paste0(SampleID,"_umap_gene_umi.pdf"),12,5)
FeaturePlot(scRNA, order = T, features = c("nCount_RNA","nFeature_RNA"),reduction = "umap")
FeaturePlot(scRNA, order = T, features = c("nCount_RNA","nFeature_RNA"),reduction = "umap",max.cutoff = "q95")
dev.off()

### count matrix
write.table(as.matrix(scRNA@assays$RNA@counts),paste0(SampleID,"_count_mat_FilterDoublet.txt"),sep="\t",quote = FALSE,row.names = T,col.names=NA)


