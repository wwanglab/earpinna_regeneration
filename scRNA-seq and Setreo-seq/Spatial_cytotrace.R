library(CytoTRACE)
library(Matrix)
library(Seurat)
library(RColorBrewer)

args <- commandArgs(T)
countrds <- args[1]
seuratrds = args[2]
pre=args[3]
outdir <- args[4]
mask=args[5]

outdir=paste0(outdir,"/",pre,"/")
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)

## seurat celltype
zhp=readRDS(seuratrds)
#zhp=zhp[,zhp$batch %in% pre]
saveRDS(zhp,file=paste0(outdir,"/test.rds"))
emb <- as.data.frame(cbind(zhp$imagecol,zhp$imagerow))
embumap <- as.data.frame(zhp[["umap"]]@cell.embeddings)
write.table(embumap,file=paste0(outdir,"/embumap.txt"),col.name=T,quote=F,sep="\t")
metadata <- zhp@meta.data
write.table(metadata,file=paste0(outdir,"/meta.txt"),quote=F)
phe <- as.vector(metadata$leiden_poly)
#phe <- as.vector(metadata$celltype)
write.table(phe,file=paste0(outdir,"/phe.txt"),quote=F)
names(phe) <- rownames(zhp@meta.data)
write.table(names(phe),file=paste0(outdir,"/phenames.txt"),quote=F)

## count rds
count=readRDS(countrds)
colnames(count@assays$Spatial@counts)=rownames(zhp@meta.data)
data <- GetAssayData(count,slot = "counts")
data <- as.matrix(data) #need about 60GB memory for 115695 cells
colnames(data)=rownames(zhp@meta.data)
write.table(rownames(zhp@meta.data),file=paste0(outdir,"/cellid.txt"),quote=F,sep="\t")
results <- CytoTRACE(mat = data)
## result
cytotracere=cbind(emb,ceiling(results$CytoTRACE*100))
cytotracere=na.omit(cytotracere)
colnames(cytotracere)=c("bin1x","bin1y","score")
write.table(cytotracere,file=paste0(outdir,"/cytotracere.txt"),row.names = F, col.names = T,quote=F,sep="\t")

## Plot
plotCytoGenes(results, numOfGenes = 10,outputDir = outdir)
plotCytoTRACE(results,phenotype = phe,outputDir = outdir,emb = embumap)
plotCytoTRACE(results,outputDir = paste0(outdir,"Spatial"),emb = emb)

#oldname=paste0(outdir,"../../","CytoTRACE_plot.pdf")
#newname=paste0(outdir,pre,"SpatialCytoTRACE_plot.pdf")
#file.rename(oldname,newname)
#file.remove(oldname)


## plot to mask
cytotracere=read.table(paste0(outdir,"/cytotracere.txt"),header = T,sep="\t")
score_quantile <- quantile(cytotracere$score,c(0.95))
cytotracere$score[cytotracere$score>score_quantile[[1]]]=score_quantile[[1]]
genefilename=paste0(outdir,"/",pre,"_score.txt")
write.table(cytotracere,genefilename,sep="\t",quote=F,row.names = F)
collen=length(unique(cytotracere$score))
pal<-rev(brewer.pal(11,'RdYlBu'))
colors=colorRampPalette(pal)(collen)
sortexp=sort(as.numeric(unique(cytotracere$score)))
collist=cbind(sortexp,colors)
colfilename=paste0(outdir,"/",pre,"_col.list")
write.table(collist,file=colfilename,quote=F,sep="\t",row.names=F,col.names=F)
tifname=paste0(outdir,"/",pre,"_score.tif")
shell=paste0("cellbin_plot  ",genefilename," ",mask," ",colfilename," ",tifname)
grep_out<-system(shell, intern = TRUE)
shell


