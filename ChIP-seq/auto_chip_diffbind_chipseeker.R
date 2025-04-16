library(DiffBind)
library(GenomicFeatures)
library(ChIPseeker)
library(dplyr)
library(genomation)

##### load all samples and generate a correlation heatmap ########## 
sample.sheet = read.csv("sample_sheet.csv")
MyChip.DBA = dba(sampleSheet = sample.sheet)

MyChip.DBA

# draw all sample correlation plot
pdf("Allsample_correlation_plot.pdf")
# par: down, left, up, right
par(oma=c(2,0,0,2))
plot(MyChip.DBA)
dev.off()

# pick out good samples and make a correlation heatmap
sample.sheet.final = sample.sheet

# function 1: count, normalize, contrast, analyze and return dba.ana object, draw correlation plot
pre.analysis = function(input.sample.sheet, output.prefix){
  MyChip.DBA = dba(sampleSheet = input.sample.sheet)
  MyChip.DBA.counts = dba.count(MyChip.DBA, summits = 300, bUseSummarizeOverlaps = T)
  # I change summits from 200 to 300
  MyChip.DBA.counts.nor = dba.normalize(MyChip.DBA.counts, method=DBA_DESEQ2, normalize=DBA_NORM_LIB)
  MyChip.DBA.counts.con = dba.contrast(MyChip.DBA.counts.nor, categories = DBA_CONDITION, minMembers = 2)
  MyChip.DBA.counts.ana = dba.analyze(MyChip.DBA.counts.con, method = DBA_DESEQ2, bBlacklist = F, bGreylist = F)
  dba.show(MyChip.DBA.counts.ana, bContrasts = T)
  
  # calculate the peakreads and save as csv files
  info = dba.show(MyChip.DBA.counts)
  libsizes = cbind(LibReads= info$Reads, FRip = info$FRiP, PeakReads = round(info$Reads * info$FRiP))
  rownames(libsizes) = info$ID
  write.csv(info, file = paste(output.prefix, "_DBA.counts.libsizes.peakreads.summary.csv", sep = ""), row.names = T)

  # calculate the normalized library size and save as csv files
  norm <- dba.normalize(MyChip.DBA.counts.nor, bRetrieve=TRUE)
  normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors))
  rownames(normlibs) <- info$ID
  write.csv(normlibs, file = paste(output.prefix, "_NormLibSize.csv", sep = ""), row.names = T)
  
  # save a correlation heatmap base on counts as pdf file
  pdf(paste(output.prefix, "_correlation_plot.pdf", sep = ""), width = 10, height = 10)
  par(oma=c(4,0,0,4))
  plot(MyChip.DBA.counts)
  title(main = "Correlation based on counts", cex.main = 0.5)
  plot(MyChip.DBA.counts.ana)
  title(main = "Correlation based on differential bound sites", cex.main = 0.5)
  dev.off()
  
  return(MyChip.DBA.counts.ana)
  
}

# function 2: generate PCA plot, correaltion heatmap and others
auto.ploting = function(MyChip.DBA.counts.ana, output.prefix, threshold){
  # save PCA plot, MA plot, Venn plot as a pdf file
  pdf(paste(output.prefix,"_output_plot_PCA_MA_volcano_heatmap_plot.pdf", sep = ""), width = 10, height = 8)
  dba.plotVenn(MyChip.DBA.counts.ana, contrast = 1, bDB = T, bGain = T, bLoss = T, bAll = F, th = threshold)
  dba.plotPCA(MyChip.DBA.counts.ana, attributes = DBA_CONDITION, label = DBA_CONDITION)
  
  # PCA plot with differential bound site
  dba.plotPCA(MyChip.DBA.counts.ana, contrast = 1, attributes = DBA_CONDITION, label = DBA_CONDITION)
  title(main = "PCA plot with differential peaks")
  dba.plotMA(MyChip.DBA.counts.ana, th = threshold)
  dba.plotVolcano(MyChip.DBA.counts.ana, th = threshold)
  pvals = dba.plotBox(MyChip.DBA.counts.ana, th = 0.01)
  hmap = colorRampPalette(c("blue","white","red"))(n = 13)
  readscores = dba.plotHeatmap(MyChip.DBA.counts.ana, contrast = 1, correlations = F, scale = "row", colScheme = hmap)
  dev.off()
  
  # write the p.value of Boxplot in the previous figure
  write.csv(pvals, file = paste(output.prefix, "_DBA.plotBox.p.value.csv", sep = ""), row.names = T)
}

# function 3: peak report and annotation
peak.annotation = function(MyChip.DBA.counts.ana, txdb, output.prefix, threshold, gene_ID){
  MyChip.DBA.report = dba.report(MyChip.DBA.counts.ana, th = threshold)
  MyChip.DBA.annotated = annotatePeak(MyChip.DBA.report, TxDb = txdb, tssRegion = c(-1000,1000))
  
  # add annotation to DBA.report object
  MyChip.DBA.annotated.df = data.frame(MyChip.DBA.annotated)
  colnames(MyChip.DBA.annotated.df)[18] = "GeneID"
  MyChip.DBA.annotated.df.final = left_join(MyChip.DBA.annotated.df, gene_ID, by = "GeneID")
  
  write.csv(MyChip.DBA.annotated.df.final, file = paste("Diffbind.",output.prefix,"_threshold_", threshold,"_differential_peaks.annotated.csv", sep = ""), row.names = F)
  
  p2 = plotAnnoBar(MyChip.DBA.annotated)
  p3 = plotDistToTSS(MyChip.DBA.annotated, title="Distribution of transcription factor-binding loci\nrelative to TSS")
  p5 = upsetplot(MyChip.DBA.annotated, vennpie=TRUE)
  pdf(paste(output.prefix,"_threshold_",threshold , "_peak_annotation.pdf"), width = 10, height = 8)
  plotAnnoPie(MyChip.DBA.annotated)
  print(p2)
  print(p3)
  vennpie(MyChip.DBA.annotated)
  print(p5)
  dev.off()

  
  promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
  tagMatrix <- getTagMatrix(MyChip.DBA.report, windows=promoter)
  p1 = plotAvgProf(tagMatrix, xlim=c(-2000, 2000),
                   xlab="Genomic Region (5'->3')", 
                   ylab = "Read Count Frequency")
  pdf(paste(output.prefix,"_threshold_",threshold, "_genomic-wide_distribution_for_TSS.pdf", sep = ""), width = 8, height = 6)
  print(p1)
  dev.off()
  
  return(MyChip.DBA.annotated)
}

# function 4: draw peak heatmap
peak.heatmap = function(MyChip.DBA.counts.ana, input.sample.sheet, output.prefix, threshold){
  MyChip.DBA.report = dba.report(MyChip.DBA.counts.ana, th = threshold)
  MyChip.peaks.up <- resize(MyChip.DBA.report, width = 2000, fix = "center")
  bam.files <- input.sample.sheet$bamReads
  sml.3 <- ScoreMatrixList(bam.files, MyChip.peaks.up, bin.num = 80, type = "bam")
  names(sml.3) = input.sample.sheet$SampleID
  
  pdf(paste(output.prefix,"_threshold_", threshold ,".ChIP_heatmap_peaks.pdf",sep = ""), width = 10, height = 8)
  multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE, 
                  col = c("white","purple"))
  multiHeatMatrix(sml.3, xcoords = c(-2000, 2000), order = TRUE, winsorize=c(0,97), common.scale=TRUE)
  
  plotMeta(sml.3, xcoords = c(-2000, 2000), centralTend = "mean", 
           profile.names = input.sample.sheet$SampleID,
           winsorize=c(0,97))
  
  dev.off()
}


MyChip.DBA.counts.ana = pre.analysis(input.sample.sheet = sample.sheet.final, output.prefix = "mouse")
auto.ploting(MyChip.DBA.counts.ana, "mouse", threshold = 0.01)

txdb = makeTxDbFromGFF(file = "/data/Genome/mouse/GRCm39/Mus_musculus.GRCm39.104.gtf")
geneID = read.csv("/data/Genome/mouse/GRCm39/mice_ensemble104_GRCm39_IDs.csv")
MyChip.DBA.annotated.th_0.01 = peak.annotation(MyChip.DBA.counts.ana, txdb, output.prefix = "mouse", threshold = 0.01, gene_ID = geneID)
MyChip.DBA.annotated.all = peak.annotation(MyChip.DBA.counts.ana, txdb, output.prefix = "mouse", threshold = 1, gene_ID = geneID)

peak.heatmap(MyChip.DBA.counts.ana, input.sample.sheet = sample.sheet.final, output.prefix = "mouse", threshold = 0.01)
peak.heatmap(MyChip.DBA.counts.ana, input.sample.sheet = sample.sheet.final, output.prefix = "mouse", threshold = 1)

temp = read.csv("Diffbind.mouse_differential_peaks.annotated.csv")
temp2 = temp[,c(1:3)]
temp2$strands = "+"
write.table(temp2, file = "mouse.for.motif.analysis.txt", sep = "\t", col.names = F, quote = F)
write.table(temp[,c(1:3,5)], file = "mouse.for.motif.analysis_2.txt", sep = "\t", col.names = F, quote = F)

save.image(file = "./output.RData")
