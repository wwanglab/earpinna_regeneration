# "gtf_path","g",2,"character","path of gtf file",
# "id_path","i",1,"character","path of biomart id file",
# "experimental_design","e",2,"character","path of experimental design"


library(Rsamtools)
library(parallel)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)
library(biomaRt)
library(ggplot2)
library(pander)
library(edgeR)
library(plotly)
library(grid)
library(GenomicRanges)
library(reshape)
library(DT)
library(dplyr)
library(getopt)


spec = matrix(
  c(
    "gtf_path","g",2,"character","path of gtf file",
    "id_path","i",1,"character","path of biomart id file",
    "experimental_design","e",2,"character","path of experimental design"
  ),
  byrow = T, ncol = 5
)
option = getopt(spec = spec)
gtf_path = option$gtf_path
txdb_name = unlist(strsplit(gtf_path, split = "/"))
txdb_name = paste(txdb_name[length(txdb_name)],".txdb", collapse = "", sep = "")

if(file.exists(txdb_name)){
  print("txdb file exist")
  txdb = loadDb(txdb_name)
  # you need to load txdb file
} else {
  txdb = makeTxDbFromGFF(file = gtf_path)
  saveDb(txdb, txdb_name)
}

exons.per.gene = exonsBy(txdb, by = "gene")
gene.lengths = as.matrix(sum(width(reduce(exons.per.gene))))

counts = read.csv("./output_read_counts.wanglab.csv")
rownames(counts) = counts$GeneID
counts = counts[,c(2:length(counts))]

targets = read.csv(option$experimental_design)


dgelist = DGEList(counts = counts, group = targets$group)
keep_cutoff_1 = rowSums(cpm(dgelist)>3) >=1
keep_cutoff_2 = rowSums(cpm(dgelist)>3) >=2
keep_cutoff_3 = rowSums(cpm(dgelist)>3) >=3

table(keep_cutoff_1)
table(keep_cutoff_2)
table(keep_cutoff_3)

dgelist_cutoff_1 = dgelist[keep_cutoff_1,,keep.lib.sizes=FALSE]
dgelist_cutoff_1 = calcNormFactors(dgelist_cutoff_1)
dgelist_cutoff_1 <- estimateCommonDisp(dgelist_cutoff_1)
dgelist_cutoff_1 <- estimateTagwiseDisp(dgelist_cutoff_1)


dgelist_cutoff_2 = dgelist[keep_cutoff_2,,keep.lib.sizes=FALSE]
dgelist_cutoff_2 = calcNormFactors(dgelist_cutoff_2)
dgelist_cutoff_2 <- estimateCommonDisp(dgelist_cutoff_2)
dgelist_cutoff_2 <- estimateTagwiseDisp(dgelist_cutoff_2)

dgelist_cutoff_3 = dgelist[keep_cutoff_3,,keep.lib.sizes=FALSE]
dgelist_cutoff_3 = calcNormFactors(dgelist_cutoff_3)
dgelist_cutoff_3 <- estimateCommonDisp(dgelist_cutoff_3)
dgelist_cutoff_3 <- estimateTagwiseDisp(dgelist_cutoff_3)

un.group = unique(as.character(targets$group))
pair = list()
for (i in 2:length(un.group)){
  pair[i-1] = list(c(un.group[1], un.group[i]))
}


#####################
run_eTest_cutoff_1 = lapply(pair, function(x){
  exactTest(dgelist_cutoff_1, pair = x)
})
run_eTest_cutoff_1 <- lapply(run_eTest_cutoff_1, function(x){ within(x$table, {
  padj = p.adjust(x$table$PValue,method="BH")})
})
df_gene = list()
df_summary = list()

# functions
minpositive <- function(x)min(x[x > 0])
maxnegative <- function(x)max(x[x < 0])

for(i in 1:length(run_eTest_cutoff_1)){
  cat("\n\n")
  pandoc.header(paste(un.group[1],"vs.", un.group[i+1], sep = ''), level = 4)
  padj.up = run_eTest_cutoff_1[[i]][(run_eTest_cutoff_1[[i]]$padj < 0.01 & run_eTest_cutoff_1[[i]]$logFC > log2(2)),]
  padj.down = run_eTest_cutoff_1[[i]][(run_eTest_cutoff_1[[i]]$padj < 0.01 & run_eTest_cutoff_1[[i]]$logFC < log2(2)),]
  df_gene[[i]] = rbind(padj.up,padj.down)
  df_summary = c(nrow(padj.up),round(minpositive(padj.up$logFC),3),
                 round(max(padj.up$logFC),3),
                 nrow(padj.down),round(minpositive(padj.down$logFC),3),
                 round(max(padj.down$logFC),3)
  )
  pdf(paste(un.group[1]," vs. ", un.group[i+1],"cutoff_1_MAplot.pdf"),paper='special')
  p <- ggplot(data = run_eTest_cutoff_1[[i]], aes(x=run_eTest_cutoff_1[[i]]$logCPM, y=run_eTest_cutoff_1[[i]]$logFC)) +
    theme_bw() + geom_point(size=1) + ylab("log2(FC)") + xlab("log2(CPM)") +
    #geom_point(data = padj.up, aes(x=padj.up$logCPM,y=padj.up$logFC), color="red", size=1) +
    #geom_point(data = padj.down, aes(x=padj.down$logCPM,y=padj.down$logFC), color="blue", size=1) +
    labs(title=paste(un.group[1]," vs. ", un.group[i+1], ": pval < 0.01 & log2(FC) > log2(2)",sep=''))
  print(p)
  dev.off()
}
#####################

#####################
run_eTest_cutoff_2 = lapply(pair, function(x){
  exactTest(dgelist_cutoff_2, pair = x)
})
run_eTest_cutoff_2 <- lapply(run_eTest_cutoff_2, function(x){ within(x$table, {
  padj = p.adjust(x$table$PValue,method="BH")})
})
df_gene = list()
df_summary = list()

# functions
minpositive <- function(x)min(x[x > 0])
maxnegative <- function(x)max(x[x < 0])

for(i in 1:length(run_eTest_cutoff_2)){
  cat("\n\n")
  pandoc.header(paste(un.group[1],"vs.", un.group[i+1], sep = ''), level = 4)
  padj.up = run_eTest_cutoff_2[[i]][(run_eTest_cutoff_2[[i]]$padj < 0.01 & run_eTest_cutoff_2[[i]]$logFC > log2(2)),]
  padj.down = run_eTest_cutoff_2[[i]][(run_eTest_cutoff_2[[i]]$padj < 0.01 & run_eTest_cutoff_2[[i]]$logFC < log2(2)),]
  df_gene[[i]] = rbind(padj.up,padj.down)
  df_summary = c(nrow(padj.up),round(minpositive(padj.up$logFC),3),
                 round(max(padj.up$logFC),3),
                 nrow(padj.down),round(minpositive(padj.down$logFC),3),
                 round(max(padj.down$logFC),3)
  )
  pdf(paste(un.group[1]," vs. ", un.group[i+1],"cutoff_2_MAplot.pdf"),paper='special')
  p <- ggplot(data = run_eTest_cutoff_2[[i]], aes(x=run_eTest_cutoff_2[[i]]$logCPM, y=run_eTest_cutoff_2[[i]]$logFC)) +
    theme_bw() + geom_point(size=1) + ylab("log2(FC)") + xlab("log2(CPM)") +
    #geom_point(data = padj.up, aes(x=padj.up$logCPM,y=padj.up$logFC), color="red", size=1) +
    #geom_point(data = padj.down, aes(x=padj.down$logCPM,y=padj.down$logFC), color="blue", size=1) +
    labs(title=paste(un.group[1]," vs. ", un.group[i+1], ": pval < 0.01 & log2(FC) > log2(2)",sep=''))
  print(p)
  dev.off()
}
#####################

#####################
run_eTest_cutoff_3 = lapply(pair, function(x){
  exactTest(dgelist_cutoff_3, pair = x)
})
run_eTest_cutoff_3 <- lapply(run_eTest_cutoff_3, function(x){ within(x$table, {
  padj = p.adjust(x$table$PValue,method="BH")})
})
df_gene = list()
df_summary = list()

# functions
minpositive <- function(x)min(x[x > 0])
maxnegative <- function(x)max(x[x < 0])

for(i in 1:length(run_eTest_cutoff_3)){
  cat("\n\n")
  pandoc.header(paste(un.group[1],"vs.", un.group[i+1], sep = ''), level = 4)
  padj.up = run_eTest_cutoff_3[[i]][(run_eTest_cutoff_3[[i]]$padj < 0.01 & run_eTest_cutoff_3[[i]]$logFC > log2(2)),]
  padj.down = run_eTest_cutoff_3[[i]][(run_eTest_cutoff_3[[i]]$padj < 0.01 & run_eTest_cutoff_3[[i]]$logFC < log2(2)),]
  df_gene[[i]] = rbind(padj.up,padj.down)
  df_summary = c(nrow(padj.up),round(minpositive(padj.up$logFC),3),
                 round(max(padj.up$logFC),3),
                 nrow(padj.down),round(minpositive(padj.down$logFC),3),
                 round(max(padj.down$logFC),3)
  )
  pdf(paste(un.group[1]," vs. ", un.group[i+1],"cutoff_3_MAplot.pdf"),paper='special')
  p <- ggplot(data = run_eTest_cutoff_3[[i]], aes(x=run_eTest_cutoff_3[[i]]$logCPM, y=run_eTest_cutoff_3[[i]]$logFC)) +
    theme_bw() + geom_point(size=1) + ylab("log2(FC)") + xlab("log2(CPM)") +
    #geom_point(data = padj.up, aes(x=padj.up$logCPM,y=padj.up$logFC), color="red", size=1) +
    #geom_point(data = padj.down, aes(x=padj.down$logCPM,y=padj.down$logFC), color="blue", size=1) +
    labs(title=paste(un.group[1]," vs. ", un.group[i+1], ": pval < 0.01 & log2(FC) > log2(2)",sep=''))
  print(p)
  dev.off()
}
#####################
IDs = read.csv(option$id_path)
# # this need to modified

df_summary_final_cutoff_1 = do.call("cbind", run_eTest_cutoff_1)
colnames(df_summary_final_cutoff_1)<-c(paste(un.group[1],".vs.", do.call(c, lapply(2:length(un.group), function(x) rep(un.group[x], 4))), "_", colnames(df_summary_final_cutoff_1)[1:4], sep=''))
df_summary_final_cutoff_1$GeneID = rownames(df_summary_final_cutoff_1)
df_summary_final_cutoff_1 = inner_join(df_summary_final_cutoff_1, IDs, by = "GeneID")
write.csv(df_summary_final_cutoff_1, file = "output_DEG_cutoff_1.wanglab.csv", row.names = F)

df_summary_final_cutoff_2 = do.call("cbind", run_eTest_cutoff_2)
colnames(df_summary_final_cutoff_2)<-c(paste(un.group[1],".vs.", do.call(c, lapply(2:length(un.group), function(x) rep(un.group[x], 4))), "_", colnames(df_summary_final_cutoff_2)[1:4], sep=''))
df_summary_final_cutoff_2$GeneID = rownames(df_summary_final_cutoff_2)
df_summary_final_cutoff_2 = inner_join(df_summary_final_cutoff_2, IDs, by = "GeneID")
write.csv(df_summary_final_cutoff_2, file = "output_DEG_cutoff_2.wanglab.csv", row.names = F)

df_summary_final_cutoff_3 = do.call("cbind", run_eTest_cutoff_3)
colnames(df_summary_final_cutoff_3)<-c(paste(un.group[1],".vs.", do.call(c, lapply(2:length(un.group), function(x) rep(un.group[x], 4))), "_", colnames(df_summary_final_cutoff_3)[1:4], sep=''))
df_summary_final_cutoff_3$GeneID = rownames(df_summary_final_cutoff_3)
df_summary_final_cutoff_3 = inner_join(df_summary_final_cutoff_3, IDs, by = "GeneID")
write.csv(df_summary_final_cutoff_3, file = "output_DEG_cutoff_3.wanglab.csv", row.names = F)

#Multi-dimensional scaling (MDS) plots
pdf("sample_MDS_plot.pdf", paper='special')
plotMDS(dgelist)
plotMDS(dgelist, cex = 5, labels = ".")
dev.off()

pdf("sample_PCA_plot.pdf", paper='special')
plotMDS(dgelist, gene.selection = "common")
dev.off()

save.image(file = "./output.RData")

##################################
# tpm: all genes    modified by weifeng
counts.new = as.data.frame(counts)
counts.new$GeneID = rownames(counts.new)
gene.lengths.new = as.data.frame(gene.lengths)
gene.lengths.new$GeneID = rownames(gene.lengths.new)
merge_count_gene.length = inner_join(counts.new, gene.lengths.new, by = "GeneID")

i = length(colnames(merge_count_gene.length))
i = i -2
input.counts = merge_count_gene.length[,1:i]
input.gene.length = as.matrix(merge_count_gene.length[,(i+2)])
rpk.new = apply(input.counts, 2, function(x){x/(input.gene.length/1000)})
tpm.new = apply(rpk.new, 2, function(x){x/(sum(x)/1e6)})
tpm.new = as.data.frame(tpm.new)
tpm.new$GeneID <- merge_count_gene.length$GeneID
tpm.final <- merge(tpm.new, IDs, by = "GeneID")
write.csv(tpm.final,"output_wanglab_tpm_allgenes.csv", row.names = FALSE)
