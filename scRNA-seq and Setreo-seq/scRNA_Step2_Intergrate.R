#!/usr/bin/Rscript
args=commandArgs(T)

if(length(args) != 7){
	print("this script need args")
	print("Rscript *R 4 organ dim res countfilelist outdir homologousgenefile")
	quit()
}

### Seurat clustering
library(Seurat)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(cowplot)
library(future)
library(harmony)

plan()
plan("multiprocess", workers = as.numeric(args[1])) #4
plan()
options(future.globals.maxSize= 429496729600)

organ <- args[2] # keolid
dim.usage <- as.numeric(args[3]) #20
res.usage <- as.numeric(args[4]) #0.6

### 1. load data
file_list <- read.table(args[5],sep="\t") 
file_list <- file_list$V1
cat(organ,"has",length(file_list),"samples.\n")



if (length(file_list) > 1){
  cat(organ,"has",length(file_list),"qualified samples to merge.\n")

  count = list()
  for(i in file_list){
    count[[i]] <- as.data.frame(fread(i))
    colnames(count[[i]])[1] <- "ID"
    SampleID <- gsub("_count_mat_FilterDoublet.txt","",i)
    cat(SampleID,"\n")
  }

  counts <- Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="ID"),count)
  rownames(counts) = counts[,which(colnames(counts) %in% "ID")]
  counts[is.na(counts)] <- 0
  counts = counts[, -which(colnames(counts) %in% "ID")]
  cat(organ,"raw count matrix:",dim(counts),"\n")


  ### 
  dir.create(args[6], recursive = TRUE)#keloid_result
  setwd(args[6])
  homologousgenefile=args[7]
#
  ### Creat Seurat object
  scRNA = CreateSeuratObject(counts = counts, min.cells = 1, min.features=1)
  ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
  scRNA@meta.data$organ <- organ
  scRNA@meta.data$batch <- unlist(lapply(strsplit(as.character(colnames(scRNA)),":"),"[",1))
#  print(paste0("batch"),unlist(lapply(strsplit(as.character(colnames(scRNA)),":"),"[",1)))
  scRNA$batch=factor(scRNA$batch,levels=c("Rb_D0_1","Rb_D0_2","Rb_D5","Rb_D10_1","Rb_D10_2","Mmu_D0","Mmu_D3","Mmu_D10"))
  cat("\n organ")
  print(table(scRNA@meta.data$organ))
  cat("\n batch")
  print(table(scRNA@meta.data$batch))

  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
  scRNA$sptime="N"
  scRNA$sptime[grep("Rb_D0",scRNA$batch)] = "RbD0"
  scRNA$sptime[grep("Rb_D5",scRNA$batch)] = "RbD5"
  scRNA$sptime[grep("Rb_D10",scRNA$batch)] = "RbD10"
  scRNA$sptime[grep("Mmu_D0",scRNA$batch)] = "MmuD0"
  scRNA$sptime[grep("Mmu_D3",scRNA$batch)] = "MmuD3"
  scRNA$sptime[grep("Mmu_D10",scRNA$batch)] = "MmuD10"
  scRNA$species="N"
  scRNA$species[grep("Rb",scRNA$sptime)]="Rabbit"
  scRNA$species[grep("Mmu",scRNA$sptime)]="Mice"
  scRNA$time="NA"
  scRNA$time[grep("D0",scRNA$sptime)] = "D0"
  scRNA$time[grep("D3",scRNA$sptime)] = "RbD5andMiceD3"
  scRNA$time[grep("D5",scRNA$sptime)] = "RbD5andMiceD3"
  scRNA$time[grep("D10",scRNA$sptime)] = "D10"
  scRNA$time=factor(scRNA$time,levels=c("D0","RbD5andMiceD3","D10"))
  scRNA$sptime=factor(scRNA$sptime,levels=c("MmuD0","MmuD3","MmuD10","RbD0","RbD5","RbD10"))

  ### Plot raw QC plot
  pdf(paste0(organ,"_filter_QCplot.pdf"), 10, 5)
  p <- VlnPlot(scRNA, group.by="organ", pt.size=0.01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  p <- VlnPlot(scRNA, group.by="batch", pt.size=0.01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  dev.off()

  pdf(paste0(organ,"_filter_scatterplot.pdf"), 12, 5)
  plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", group.by="organ", feature2 = "percent.mt")
  plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", group.by="organ", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", group.by="batch", feature2 = "percent.mt")
  plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", group.by="batch", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
  
  saveRDS(scRNA, file = paste0(organ,"_scRNA.pre.rds"))

################################################## 1. merge ##################################################
  scRNA <- NormalizeData(object = scRNA, verbose = FALSE)
  scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  scRNA <- ScaleData(scRNA) 
  scRNA <- RunPCA(scRNA,verbose=F)

  ### Find cluster
  scRNA <- FindNeighbors(scRNA, dims = 1:dim.usage)
  scRNA <- FindClusters(scRNA, resolution = res.usage)
  scRNA <- RunUMAP(scRNA, dims = 1:dim.usage)

  saveRDS(scRNA, file = paste0(organ,"_scRNA.rds"))


  ### Cluster_Batch
  pdf(paste0(organ,"_Cluster.pdf"),8,7)
  p <- DimPlot(object = scRNA, reduction = "umap",pt.size = 0,label=T)+ggtitle(label = "umap_Cluster")
  print(p)
  p <- DimPlot(object = scRNA, reduction = "umap",group.by = "organ",pt.size = 0,label=F)+ggtitle(label = "umap_organ")
  print(p)
  p <- DimPlot(object = scRNA, reduction = "umap",group.by = "batch",pt.size = 0,label=F)+ggtitle(label = "umap_batch")
  print(p)
  dev.off()

  pdf(paste0(organ,"_Cluster_batchSplit.pdf"),length(unique(scRNA@meta.data$batch))*5+1,6)
  p <- DimPlot(scRNA, reduction = "umap", split.by = "batch",pt.size = 0.0,label=F)+ggtitle(label = "umap_batchSplit")
  print(p)
  dev.off()

  ### umi_gene
  pdf(paste0(organ,"_UMAP_gene_umi.pdf"),12,5)
  p <- FeaturePlot(scRNA, order = T, features = c("nCount_RNA","nFeature_RNA"),reduction = "umap")
  print(p)
  p <- FeaturePlot(scRNA, order = T, features = c("nCount_RNA","nFeature_RNA"),reduction = "umap",max.cutoff="q95")
  print(p)
  dev.off()


  ### meta data
  write.table(scRNA@meta.data,file = paste0(organ,"_meta_data.txt"),sep="\t",quote=FALSE)


################################################## 2. integrated ##################################################
  if (length(unique(scRNA@meta.data$batch)) > 1){
    cat(organ,"has",length(file_list),"(",file_list,")","qualified samples from",length(unique(scRNA@meta.data$batch)),"(",unique(scRNA@meta.data$batch),") batchs to integrate.\n")
    
    scRNA.list <- SplitObject(object = scRNA, split.by = "batch")

    for (i in 1:length(x = scRNA.list)) {
        scRNA.list[[i]] <- NormalizeData(object = scRNA.list[[i]], verbose = FALSE)
        scRNA.list[[i]] <- FindVariableFeatures(object = scRNA.list[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    }

    reference.list <- scRNA.list
    features <- SelectIntegrationFeatures(object.list = reference.list,nfeatures = 4000)
    homologousgene=read.table(homologousgenefile,header=T)
    interfeatures = intersect(features,homologousgene$mmu)
    write.table(interfeatures,paste0(organ,"interfeature.txt"),quote=F)
    scRNA.anchors <- FindIntegrationAnchors(object.list = reference.list,anchor.features=interfeatures, dims = 1:dim.usage)
    scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors, dims = 1:dim.usage)

    DefaultAssay(object = scRNA.integrated) <- "integrated"
    scRNA.integrated <- ScaleData(object = scRNA.integrated, verbose = FALSE)
    scRNA.integrated <- RunPCA(object = scRNA.integrated, npcs = dim.usage, verbose = FALSE)
    scRNA.integrated <- RunUMAP(object = scRNA.integrated, reduction = "pca", dims = 1:dim.usage)

    scRNA.integrated <- FindNeighbors(object = scRNA.integrated,dims = 1:dim.usage)
    scRNA.integrated.Find <- FindClusters(object = scRNA.integrated,resolution = res.usage)
    scRNA.integrated.Find <- RunTSNE(object = scRNA.integrated.Find,dims = 1:dim.usage)
  scRNA.integrated.Find=scRNA
    DefaultAssay(scRNA.integrated.Find)="RNA"
  scRNA.integrated.Find$time=as.factor(scRNA.integrated.Find$time)
  scRNA.integrated.Find$time=factor(scRNA.integrated.Find$time,levels=c("D0","RbD5andMiceD3","D10"))
  scRNA.integrated.Find$sptime=as.factor(scRNA.integrated.Find$sptime)
  scRNA.integrated.Find$sptime=factor(scRNA.integrated.Find$sptime,levels=c("MmuD0","MmuD3","MmuD10","RbD0","RbD5","RbD10"))
  scRNA.integrated.Find$batch=as.factor(scRNA.integrated.Find$batch)
  scRNA.integrated.Find$batch=factor(scRNA.integrated.Find$batch,levels=c("Rb_D0_1","Rb_D0_2","Rb_D5","Rb_D10_1","Rb_D10_2","Mmu_D0","Mmu_D3","Mmu_D10"))

#    saveRDS(scRNA.integrated.Find, file = paste0(organ,"_scRNA.integrated.Find.rds"))
      
    ### Cluster_Batch
    pdf(paste0(organ,"_Cluster_integrated.pdf"),8,7)
    p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",pt.size = 0,label=T,cols=allcolour)+ggtitle(label = "umap_Cluster")
    print(p)
    p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",group.by = "organ",pt.size = 0,label=F,cols=allcolour)+ggtitle(label = "umap_organ")
    print(p)
    p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",group.by = "batch",pt.size = 0,label=F)+ggtitle(label = "umap_batch")
    print(p)
    p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",group.by = "sptime",pt.size = 0,label=F)+ggtitle(label = "umap_batch")
    print(p)
    dev.off()

    pdf(paste0(organ,"_Cluster_batchSplit_integrated.pdf"),length(unique(scRNA@meta.data$batch))*5+1,6)
    p <- DimPlot(scRNA.integrated.Find, reduction = "umap", split.by = "batch",pt.size = 0,label=F,cols=allcolour)+ggtitle(label = "umap_batchSplit")
    print(p)
    p <- DimPlot(scRNA.integrated.Find, reduction = "umap", split.by = "sptime",pt.size = 0,label=F,cols=allcolour)+ggtitle(label = "umap_batchSplit")
    print(p)
    dev.off()

    ### umi_gene
    pdf(paste0(organ,"_UMAP_gene_umi_integrated.pdf"),12,5)
    p <- FeaturePlot(scRNA.integrated.Find, order = T, features = c("nCount_RNA","nFeature_RNA"),reduction = "umap")
    print(p)
    p <- FeaturePlot(scRNA.integrated.Find, order = T, features = c("nCount_RNA","nFeature_RNA"),reduction = "umap",max.cutoff="q95")
    print(p)
    dev.off()



    ### meta data
    write.table(scRNA.integrated.Find@meta.data,file = paste0(organ,"_meta_data_integrated.txt"),sep="\t",quote=FALSE)

  }else{
    cat(organ,"has",length(file_list),"(",file_list,")","qualified samples from",length(unique(scRNA@meta.data$batch)),"(",unique(scRNA@meta.data$batch),") batchs.\n")
  }

} else {
  cat(organ,"has",length(file_list),"qualified samples:",file_list,".\n")
}

##
################### StackedPlot
scRNA=scRNA.integrated.Find
subcluster=c(4,6,9,5,22,17,7,25)
scRNA=scRNA[,scRNA$seurat_clusters %in% subcluster]
cellratio <- prop.table(table(Idents(scRNA),scRNA$sptime),margin=2)
cellratio <- as.data.frame(cellratio)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
pdf(paste0(organ,"StackeRatio.pdf"))
p <- ggplot(cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity") + theme_classic() + labs(x='Time',y = 'Ratio') + scale_fill_manual(values = allcolour) + theme(panel.border = element_rect(fill=NA,color="black", size=0, linetype="blank"))
print(p)
p <- ggplot(cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity") + theme_classic() + labs(x='Time',y = 'Ratio') + scale_fill_manual(values = c('4'="#9370DB",'6'="#F08080",'9'="#FFFF00",'5'="#98FB98",'22'="#FF1493",'17'="#D2B48C",'7'="#1E90FF",'25'="#FFE4B5")) + theme(panel.border = element_rect(fill=NA,color="black", size=0, linetype="blank"))
print(p)
dev.off()

################### clustree
library(clustree)
DefaultAssay(object = scRNA) <- "integrated"
scRNA <- FindClusters(scRNA,resolution = c(seq(0.6,1.2,.1)))
pdf(paste0(organ,"clustree.pdf"))
p=clustree(scRNA@meta.data, prefix = "integrated_snn_res.")
print(p)
dev.off()

