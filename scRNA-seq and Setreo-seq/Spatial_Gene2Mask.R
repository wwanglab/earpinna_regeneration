args=commandArgs(T)
if(length(args) != 5){ 
	print("this script need args")
	print("Rscript *R rds mask genelist outdir key(rds keywords)") 
	quit()
}

library(Seurat)
library(stringr)
indir=args[1]
mask=args[2]
genelist=args[3]
rawoutdir=args[4]
key=args[5]
if(!dir.exists(rawoutdir)){dir.create(rawoutdir,recursive=T)}

rdsfile=list.files(path=indir,pattern=paste0(key,"_*"))
gene=read.table(genelist)
gene=gene$V1;gene=tolower(gene);gene=str_to_title(gene)
for (rds in rdsfile){
	slicepre=gsub("_imputated_SpaGCN.h5ad.rds","",rds)
	outdir=paste0(rawoutdir,"/",slicepre)
	if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}
	sc=readRDS(paste0(indir,"/",rds))
	DefaultAssay(sc)="RNA"
	sccount=sc@assays$RNA@counts
	allgene=rownames(sc@assays$RNA@counts)
	intergene=intersect(x=gene,y=allgene)
	setdiff(gene,allgene)
	for (g in intergene){
		exp=ceiling((as.data.frame(sccount[g,]))*1000000)
		test=cbind(sc$imagecol,sc$imagerow,exp)
		colnames(test)=c("bin1x","bin1y","exp")
		## 去掉极值
		score_quantile <- quantile(test$exp,c(0.95))
		test$exp[test$exp>score_quantile[[1]]]=score_quantile[[1]]
		genefilename=paste0(outdir,"/",slicepre,"_",g,"_exp.txt")
		write.table(test,genefilename,sep="\t",quote=F,row.names = F)
		
		collen=length(unique(test$exp))
		colors=colorRampPalette(c("#4B0082","#9E0142","#FFFF00"))(collen)
		sortexp=sort(as.numeric(unique(test$exp)))
		collist=cbind(sortexp,colors)
		colfilename=paste0(outdir,"/",slicepre,"_",g,"col.list")
		write.table(collist,file=colfilename,quote=F,sep="\t",row.names=F,col.names=F)
		
		tifname=paste0(outdir,"/",slicepre,"_",g,"_exp.tif")
		shell=paste0("cellbin_plot  ",genefilename," ",mask," ",colfilename," ",tifname)
		grep_out<-system(shell, intern = TRUE)

		## 每张图生成一个图例
#		col_fun = colorRamp2(c(-1, 0, 1), c("#4B0082", "#9E0142", "#FFFF00"))
#		lgd <- Legend(col_fun = col_fun, title = "foo", at = c(min(sortexp),ceilling((max(sortexp)-min(sortexp))/2),max(sortexp)))
#		pdf(paste0(outdir,"/",slicepre,"_",g,"legend"))
#		daw(lgd)
#		dev.off()
	}
}

### 图例
#library(ComplexHeatmap)
#library(circlize)
#col_fun = colorRamp2(c(-1, 0, 1), c("#4B0082", "#9E0142", "#FFFF00"))
#lgd <- Legend(col_fun = col_fun, title = "exp",at = c(-1, 1),labels = c("low", "high"))
#pdf(paste0(outdir,"/legend.pdf"))
#draw(lgd)
#dev.off()
#lgd <- Legend(col_fun = col_fun, title = "exp",at = c(-1, 1),labels = c("low", "high"),direction = "horizontal")
#pdf(paste0(outdir,"/legendhorizontal.pdf"))
#draw(lgd)
#dev.off()
#
