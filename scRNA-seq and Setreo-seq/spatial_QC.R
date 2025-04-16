args=commandArgs(T)
if(length(args) != 5){
	print("this script need args")
	print("Rscript *R infile.gem filtgene cellmask outdir ")
	quit()
}

library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
library(ggsci)

library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ggpubr)
library(stringr)


infile = args[1] ##gem file
flitG = as.numeric(args[2]) ## 0
cellmask = args[3] ## cellbin maskfile
res=as.numeric(args[4]) ## 0.3
outdir = args[5]
bs="_singlecell"
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}

pro=tail(unlist(strsplit(infile,"/")),1)
pro=gsub("mat.|.txt|.tsv|.gz|_filterd|.gem","",pro)
options(future.globals.maxSize = 100000000 * 1024^2,scipen = 20000000)


dir.create(outdir,recursive = TRUE)
setwd(outdir)
dir.create("clusters",recursive = TRUE)
############################## 1. bin data  ##############################
dat <- fread(file = infile,sep="\t")
dat1 <- dat[,1:4];dat2=dat[,6];dat=cbind(dat1,dat2)
colnames(dat) <- c("geneID","x","y","UMICount","label")
out <- as.data.frame(dat)

dat.x<-as.data.frame(dat[, ceiling(median(x)), by = .(label)])
hash.x<-data.frame(row.names = dat.x$label, values = dat.x$V1)
dat.y<-as.data.frame(dat[, ceiling(median(y)), by = .(label)])
hash.y<-data.frame(row.names = dat.y$label, values = dat.y$V1)
dat$x<-hash.x[sprintf('%d', dat$label), 'values']
dat$y<-hash.y[sprintf('%d', dat$label), 'values']

#check duplicated data
te <- data.frame(unique(paste(dat$x,dat$y,dat$label,sep="_")))
colnames(te) <- "xyb"
split_b<-str_split(te$xyb,"_")
te$x <- as.vector(sapply(split_b,"[",1))
te$y <- as.vector(sapply(split_b,"[",2))
te$b <- as.vector(sapply(split_b,"[",3))
te$xy <- paste(te$x,te$y,sep="_")
dp <- which(duplicated(te$xy))
dplen <- length(dp)

for (i in dplen) {
  tmp <- dp[i]
  xy <- te[tmp,]$xy
  cxy <- which(te$xy %in% xy)
  cxylen <- length(cxy)
  keep <- te[cxy[1],]$b
  for (m in 2:cxylen) {
    c <- which(dat$label %in% te[cxy[m],]$b)
    #dat$label[grep(as.character(te[cxy[m],]$b),dat$label)] <- keep
    dat$label[c] <- as.numeric(keep)
  }
}

dat$label <- as.numeric(dat$label)

infoout <- cbind(dat$y,dat$x,out)
colnames(infoout)[1:2] <- c(paste0("cell",".y"),paste0("cell",".x"))
fwrite(infoout,paste0(pro,"_cell","_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)


dat <- dat[, sum(UMICount), by = .(geneID, x, y,label)]
#dat$bin_ID <- dat$label
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
bin.coor <- dat[, sum(V1), by = .(x, y)]

out <- as.data.frame(cbind(paste0("CELL.",unique(dat$bin_ID)),bin.coor$y,bin.coor$x))
colnames(out) <- c(paste0("CELL.",bs),paste0("bin",bs,".y"),paste0("bin",bs,".x"))
rownames(out) <- out[,1]
fwrite(out,paste0(pro,"_bin",bs,"_position.txt"),col.names=T,row.names=F,sep="\t",quote=F)
bkpos <- out

##
geneID <- seq(length(unique(dat$geneID)))
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, 'values']


##
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']


##
cnt <- dat$V1


##
rm(dat)
gc()


##
tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))

tissue_positions_list <- data.frame(row.names = paste('CELL', rownames(hash.B), sep = '.'),
                                    tissue = 1, 
                                    row = bin.coor$y, col = bin.coor$x,
                                    imagerow = bin.coor$y, imagecol = bin.coor$x)

scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))


##
mat <- sparseMatrix(i = gen, j = bin, x = cnt)

rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste('CELL', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')
############################ 2. creat Spatial Object  #####################
seurat_spatialObj <- CreateSeuratObject(mat, project = pro, assay = 'Spatial',min.cells=1, min.features=1)

##
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) 
{
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  
  spot.radius <- unnormalized.radius / max(dim(image))
  
  return(new(Class = 'VisiumV1', 
             image = image, 
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                          fiducial = scale.factors$fiducial_diameter_fullres, 
                                          hires = scale.factors$tissue_hires_scalef, 
                                          lowres = scale.factors$tissue_lowres_scalef), 
             coordinates = tissue.positions, 
             spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                  scale.factors = fromJSON(scalefactors_json), 
                                  tissue.positions = tissue_positions_list)


##
spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'
#DefaultAssay(seurat_spatialObj) <- 'Spatial'

seurat_spatialObj[['slice1']] <- spatialObj


rm(mat)
rm(bin.coor)
rm(hash.G)
rm(hash.B)
rm(bin)
rm(gen)
rm(cnt)

dir.create("figures")
##############################  3. Spatial Analyse  ##############################
### for mouse
seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, pattern = "^mt-") 

Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3+1.5*(Q3-Q1))
lower <- as.numeric(Q1-1.5*(Q3-Q1))

mean_gene<-as.integer(mean(seurat_spatialObj$nFeature_Spatial))
median_gene<-as.integer(median(seurat_spatialObj$nFeature_Spatial))
mean_UMI<-as.integer(mean(seurat_spatialObj$nCount_Spatial))
median_UMI<-as.integer(median(seurat_spatialObj$nCount_Spatial))
mean_mito<-as.integer(mean(seurat_spatialObj$percent.mt))
median_mito<-as.integer(median(seurat_spatialObj$percent.mt))
bin_number<-dim(seurat_spatialObj)[2]

pdf(paste0("figures/",pro,"_bin",bs,"_preQC.pdf"),width=8,height=8)
p1 <- VlnPlot(seurat_spatialObj, features=c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3, pt.size=0)
print(p1)
p2 <- ggplot(seurat_spatialObj@meta.data,aes(x=nFeature_Spatial)) +geom_density(colour="black") + theme_classic() + theme(plot.title=element_text(hjust=0.5,size=18, face="bold.italic"), legend.position="none",axis.title=element_text(size=15, face="bold.italic"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank()) + 
  geom_vline(aes(xintercept=flitG,colour="#999999",linetype="twodash")) + geom_vline(aes(xintercept=lower, colour="#377EB8",linetype="twodash")) + geom_vline(aes(xintercept=upper, colour="#E41A1C", linetype="twodash"))+xlim(min(seurat_spatialObj@meta.data$nFeature_Spatial),max(seurat_spatialObj@meta.data$nFeature_Spatial)) + 
  ggtitle(paste0(pro,".BIN_",bs,":","\n","nGene:",dim(seurat_spatialObj)[1],"; ","nBIN:",dim(seurat_spatialObj)[2],"\n","Gene:",mean_gene," ",median_gene,"\n","UMI:",mean_UMI," ",median_UMI,"\n","Mt ratio:",mean_mito," ",median_mito))
print(p2)
p3 <- SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial") + theme(legend.position = "right") + scale_y_reverse()
print(p3)
p4 <- SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial") + theme(legend.position = "right") + scale_y_reverse()
print(p4)
dev.off()

png(file=paste0("figures/",pro,"_bin",bs,"_preQC.png"), width=1000,height=1000,res=80)
print(cowplot::plot_grid(p1, p2, p3, p4, rel_widths=c(1,1,1,1), ncol=2))
dev.off()

saveRDS(seurat_spatialObj,file=paste0(pro,"_bin",bs,"_preQC.rds"))

###################################  4. quality check  #################################
seurat_spatialObj <- readRDS(paste0(pro,"_bin",bs,"_preQC.rds"))
seurat_spatialObj <- subset(seurat_spatialObj,nFeature_Spatial>flitG)
dir.create("figures_afterqc")

Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3+1.5*(Q3-Q1))
lower <- as.numeric(Q1-1.5*(Q3-Q1))

mean_gene<-as.integer(mean(seurat_spatialObj$nFeature_Spatial))
median_gene<-as.integer(median(seurat_spatialObj$nFeature_Spatial))
mean_UMI<-as.integer(mean(seurat_spatialObj$nCount_Spatial))
median_UMI<-as.integer(median(seurat_spatialObj$nCount_Spatial))
mean_mito<-as.integer(mean(seurat_spatialObj$percent.mt))
median_mito<-as.integer(median(seurat_spatialObj$percent.mt))
bin_number<-dim(seurat_spatialObj)[2]

pdf(paste0("figures_afterqc/",pro,"_bin",bs,"_afterQC.pdf"),width=8,height=8)
p1 <- VlnPlot(seurat_spatialObj, features=c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3, pt.size=0)
print(p1)
p2 <- ggplot(seurat_spatialObj@meta.data,aes(x=nFeature_Spatial)) +geom_density(colour="black") + theme_classic() + theme(plot.title=element_text(hjust=0.5,size=18, face="bold.italic"), legend.position="none",axis.title=element_text(size=15, face="bold.italic"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank()) + 
  geom_vline(aes(xintercept=flitG,colour="#999999",linetype="twodash")) + geom_vline(aes(xintercept=lower, colour="#377EB8",linetype="twodash")) + geom_vline(aes(xintercept=upper, colour="#E41A1C", linetype="twodash"))+xlim(min(seurat_spatialObj@meta.data$nFeature_Spatial),max(seurat_spatialObj@meta.data$nFeature_Spatial)) + 
  ggtitle(paste0(pro,".BIN_",bs,":","\n","nGene:",dim(seurat_spatialObj)[1],"; ","nBIN:",dim(seurat_spatialObj)[2],"\n","Gene:",mean_gene," ",median_gene,"\n","UMI:",mean_UMI," ",median_UMI,"\n","Mt ratio:",mean_mito," ",median_mito))
print(p2)
p3 <- SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial") + theme(legend.position = "right") + scale_y_reverse()
print(p3)
p4 <- SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial") + theme(legend.position = "right") + scale_y_reverse()
print(p4)
dev.off()

png(file=paste0("figures_afterqc/",pro,"_bin",bs,"_afterQC.png"), width=1000,height=1000,res=80)
print(cowplot::plot_grid(p1, p2, p3, p4, rel_widths=c(1,1,1,1), ncol=2))
dev.off()
################################  4. Data SCTransform  Clusters ##############################
setwd("clusters")
dir.create("figures",recursive = TRUE)

set.seed(6)

seurat_spatialObj <- SCTransform(seurat_spatialObj, assay = "Spatial", verbose = FALSE, variable.features.n=2000)
seurat_spatialObj <- RunPCA(seurat_spatialObj,assay = "SCT",verbose = F)
seurat_spatialObj <- FindNeighbors(seurat_spatialObj, reduction = "pca", dims = 1:20)
seurat_spatialObj <- RunUMAP(seurat_spatialObj, reduction = "pca", dims = 1:20)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE,resolution = res)

col <- c(brewer.pal(9,"Set1")[c(1:5,7,6,8,9)],brewer.pal(9,"Pastel1")[c(7,2,3,5,1)],brewer.pal(12,"Set3")[c(-2,-4,-6,-9,-11)],brewer.pal(8,"Set2")[c(3,4,7)],brewer.pal(8,"Pastel2")[8],brewer.pal(12,"Paired")[c(5,6,9,10,12)])
if(length(levels(seurat_spatialObj))>length(col)){
	  col <- colorRampPalette(col)(length(levels(seurat_spatialObj)))}else{
		    col <- col[1:length(levels(seurat_spatialObj))]
}

plot2 <- ElbowPlot(seurat_spatialObj, ndims=50, reduction="pca")
pdf(paste0("figures/",pro,"_bin",bs,"_SCT_UMAP.pdf"))
print(plot2)
DimPlot(seurat_spatialObj, reduction="umap", label=TRUE,cols = col)
#SpatialDimPlot(seurat_spatialObj, label = TRUE, label.size=3, cols=col, stroke=0, pt.size.factor=pt.s) + scale_y_reverse()
dev.off()
pdf(paste0("figures/",pro,"_bin",bs,"_SCT_UMAP_split.pdf"), length(levels(seurat_spatialObj))*5, 6)
DimPlot(seurat_spatialObj, reduction="umap", label=TRUE, split.by = "seurat_clusters", cols = col)
dev.off()

seurat_spatialObj <- NormalizeData(seurat_spatialObj, normalization.method = "LogNormalize", assay="Spatial")
seurat_spatialObj <- ScaleData(seurat_spatialObj, assay="Spatial", features =rownames(seurat_spatialObj))

######  DEG of patial.assay
de_markers <- FindAllMarkers(seurat_spatialObj,only.pos = FALSE, assay="Spatial", min.pct = 0.1, logfc.threshold = 0.25)
de_markers <- de_markers[de_markers$p_val_adj<0.05,]
de_markers <- de_markers[order(de_markers$cluster,-de_markers$avg_log2FC),]
de_markers$avg_FC <- 2^de_markers$avg_log2FC

top10 <- de_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)

pdf(paste0("figures/",pro,"_bin",bs,"_Cluster_top10Heat.pdf"))
DoHeatmap(seurat_spatialObj, assay="Spatial", features=top10$gene) + NoLegend()
#DoHeatmap(seurat_spatialObj, assay="Spatial", features=genes) + NoLegend()
dev.off()
write.table(de_markers, file=paste0(pro,"_bin",bs,"_DEmarker.txt"),sep = "\t",quote = F)
write.table(top10, file=paste0(pro,"_bin",bs,"_top10DEmarker.txt"),sep = "\t",quote = F)

#seurat_spatialObj@assays$Spatial@scale.data <- seurat_spatialObj@assays$Spatial@scale.data[1:2,1:2]
saveRDS(seurat_spatialObj,file=paste0(pro,"_bin",bs,"_Cluster.rds")) # zhangpanyu

### save info for customizing SpatialDimPlot 
SpatialPlotInfo <- cbind(seurat_spatialObj@meta.data,seurat_spatialObj@images$slice1@coordinates)
write.table(SpatialPlotInfo,paste0(pro,"_bin",bs,"_SpatialPlotInfo.txt"),quote=F,sep="\t",row.names=T,col.names=NA)
write.table(col,paste0(pro,"_bin",bs,"_col.txt"),quote=F,sep="\t",row.names=F,col.names=F)

###
pdf(paste0("figures/",pro,"_bin",bs,"_cluster_vln_QC.pdf"),6,6)
plot_grid(VlnPlot(seurat_spatialObj, features = 'nFeature_Spatial', pt.size = 0, col=col)+xlab("")+NoLegend(), VlnPlot(seurat_spatialObj, features = 'nCount_Spatial', pt.size = 0, col=col)+xlab("")+ NoLegend(), VlnPlot(seurat_spatialObj, features = 'percent.mt', pt.size = 0, col=col)+xlab("")+ NoLegend(), ncol=2)
dev.off()

############################################# cellbinPlot
#colnames(infoout)=c("cell.x","cell.y","geneID","bin1x","bin1y","UMICount","label")
#colnames(SpatialPlotInfo)[13:14]=c("cell.x","cell.y")
#merge_table=merge(infoout,SpatialPlotInfo,by = .(cell.x,cell.y),all=T)
#data=cbind(merge_table$bin1x,merge_table$bin1y,merge_table$seurat_clusters)
#
##clusters=as.character((seurat_spatialObj$seurat_clusters))
##data=cbind(seurat_spatialObj@images$slice1@coordinates$imagerow,seurat_spatialObj@images$slice1@coordinates$imagecol,clusters)
#data=unique(data)
#colnames(data)=c("bin1x","bin1y","cluster")
#write.table(data,"bin1clu.txt",sep="\t",quote=F,row.names=F,col.names=T)

collist=read.csv(paste0(pro,"_bin",bs,"_col.txt"),header=F)
colnames(collist)=c("col")
collist$clu=as.numeric(rownames(collist))-1
collist=cbind(collist$clu,collist$col)
write.table(collist,"col.list",sep="\t",quote=F,row.names=F,col.names=F)

#shell_cmd<-paste0("/data/chenxi/local/bin/cell_bin_plot bin1clu.txt ",cellmask," col.list plot.tif")
##grep_out<-system(shell_cmd, intern = TRUE)
#cat(grep_out)




###################################  5. data preprocessing  ##################################
#seurat_spatialObj@meta.data$sample <- seurat_spatialObj@meta.data$orig.ident
#seurat_spatialObj@meta.data$age <- age
#seurat_spatialObj@meta.data$tissue <- tissue
#names(seurat_spatialObj@images) <- pro
#seurat_spatialObj <- RenameCells(seurat_spatialObj, new.names=paste0(colnames(seurat_spatialObj),"_",pro))
#saveRDS(seurat_spatialObj,file=paste0(pro,"_bin",bs,"_afterpreprocessing.rds"))
