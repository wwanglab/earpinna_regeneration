args=commandArgs(T)

# 加载必要的库
library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(stringr)

process_gene <- function(gene) {
  if (str_starts(gene, "ENSOCUG")) {
      return(gene)
  } else {
      # 将字符串转换为首字母大写，其他字母小写
      return(str_to_title(str_to_lower(gene)))
    }
}


sc=readRDS("EarCelltypescRNAAnnoCluV1.rds")
DefaultAssay(sc)="RNA"
Idents(sc)=sc$papercelltype
seurat_obj=sc

genefile=args[1]
outdir=args[2]
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}
genelist=read.table(genefile,sep="\t",header=F)
genelist=genelist$V1
processed_genelist <- sapply(genelist, process_gene)

##################################### mmu & rb(两个物种一起) ########################
######################################## step1 热图 
# 使用AverageExpression获取每个细胞类型的平均表达
avg_expression <- AverageExpression(seurat_obj, features = genelist, group.by = "ident")
# 提取RNA assay的结果（如果您使用的是默认的RNA assay）
avg_expression <- avg_expression$RNA
# 计算每个细胞类型中三个基因的平均表达
cell_type_avg <- colMeans(avg_expression)
# 创建一个数据框来存储结果
result_df <- data.frame(
  CellType = names(cell_type_avg),
    AvgExpression = cell_type_avg
  )
# 对平均表达值进行排序（可选）
result_df <- result_df[order(result_df$AvgExpression, decreasing = TRUE), ]
# 准备热图数据
heatmap_data <- matrix(result_df$AvgExpression, ncol = 1)
rownames(heatmap_data) <- result_df$CellType
colnames(heatmap_data) <- "AvgExpression"
# 绘制热
#pheatmap(heatmap_data,
#          cluster_rows = FALSE,  # 不对行进行聚类
#           cluster_cols = FALSE,  # 不对列进行聚类
#            show_colnames = FALSE, # 不显示列名
#         main = "Average Expression of Selected Genes Across Cell Types",
#          fontsize = 10,
#           cellwidth = 25,
#            cellheight = 15,
#         color = colorRampPalette(c("blue", "white", "red"))(100))
pdf(paste0(outdir,"/MmuRb_AvgExpression.pdf"))
color_palette <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
p=pheatmap(heatmap_data,
   fontsize = 10,
     fontsize_row = 10,
   fontsize_col = 10,
     cellwidth = 30,
   cellheight = 20,
     color = color_palette,
   cluster_rows = FALSE,
     cluster_cols = FALSE,
   show_colnames = FALSE,
     angle_col = 0,
   main = "",
     border_color = NA,
   legend = TRUE,
     legend_breaks = seq(min(heatmap_data), max(heatmap_data), length.out = 3),
   legend_labels = c("Low", "Medium", "High")
   )
print(p)
dev.off()

########## 显著性统计 
##  计算每个细胞平均表达
genes_of_interest=genelist
genes_of_interest=intersect(genes_of_interest,rownames(seurat_obj@assays$RNA@data))
seurat_obj$avg_expression <- colMeans(seurat_obj@assays$RNA@data[genes_of_interest, ])
data = as.data.frame(cbind(seurat_obj$papercelltype,seurat_obj$avg_expression))
colnames(data)=c("papercelltype","avg_expression")
## (使用Wilcoxon秩和检验)
# 首先，确保你的数据是按照papercelltype分组的
data$papercelltype <- as.factor(data$papercelltype)
# 获取除了"WIF"之外的所有细胞类型
other_cell_types <- setdiff(levels(data$papercelltype), "WIF")
# 创建一个空列表来存储所有的检验结果
test_results <- list()
# 对"WIF"细胞类型与其它每个细胞类型进行Wilcoxon秩和检验
for(ct in other_cell_types) {
  # 提取"WIF"细胞类型的表达数据
  wif_expression <- as.numeric(data$avg_expression[data$papercelltype == "WIF"])
  # 提取当前细胞类型的表达数据
  other_expression <- as.numeric(data$avg_expression[data$papercelltype == ct])
    # 执行Wilcoxon秩和检验
    test_result <- wilcox.test(wif_expression, other_expression, alternative = "two.sided")
    ## 计算 LogFC
    avg_LogFC <- log2(mean(wif_expression / mean(other_expression)))
  # 保存检验结果
  test_results[[paste("WIF_vs_", ct, sep = "")]] <- list(
     statistic = test_result$statistic,
     p.value = test_result$p.value,
     avg_LogFC = avg_LogFC
   )
}
# 将检验结果转换为数据框以便查看
test_results_df <- do.call(rbind, test_results) %>% as.data.frame()
names(test_results_df) <- c("Statistic", "P.Value","avg_LogFC")
# 查看结果
print(test_results_df)
# 手动构建数据框
Statistic <- sapply(test_results, function(x) x$statistic)
P_Value <- sapply(test_results, function(x) x$p.value)
avg_LogFC <- sapply(test_results, function(x) x$avg_LogFC)

test_results_df <- data.frame(
    Cell_Type = names(test_results),
  Statistic = Statistic,
    P.Value = P_Value,
  avg_LogFC = avg_LogFC
)

write.table(test_results_df,paste0(outdir,"/MmuRb_wilcox.txt"),quote=F,row.names=T,col.names=T,sep="\t")
######################### vlnplot
my_comparisons=list(c("WIF","Fibroblast"),c("WIF","Immune"),c("WIF","Melanocytes"),c("WIF","Keratinocyte"),c("WIF","Muscle"),c("WIF","Chondrocyte"),c("WIF","Endothelium"))
plotfile=cbind(seurat_obj$papercelltype,seurat_obj$avg_expression)
colnames(plotfile)=c("celltype","score")
plotfile=as.data.frame(plotfile)
plotfile$score=as.numeric(plotfile$score)
plotfile <- plotfile %>%
  group_by(celltype) %>%
  mutate(mean_score = mean(score, na.rm = TRUE)) %>%
  ungroup()
# 使用reorder()根据平均score对time进行排序，从大到小
plotfile$celltype <- reorder(plotfile$celltype, plotfile$mean_score, FUN = max, na.rm = TRUE)
pdf(paste0(outdir,"/","MmuRb_VlnPlot.pdf"))
p = ggplot(data = plotfile, aes(x = celltype, y = score, color = celltype)) +
  geom_violin() +
    geom_boxplot(width = 0.1) +
  scale_color_brewer(palette = "Set1") + 
    theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = .9, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
  scale_y_continuous(limits = c(0, max(plotfile$score, na.rm = TRUE) + 2)) +
    ylab("AvgExp")
  print(p)
dev.off()



#################################### mmu (only mmu)########################
# 使用AverageExpression获取每个细胞类型的平均表达
seurat_obj=sc[,sc$species=="Mice"]
avg_expression <- AverageExpression(seurat_obj, features = genelist, group.by = "ident")

# 提取RNA assay的结果（如果您使用的是默认的RNA assay）
avg_expression <- avg_expression$RNA

# 计算每个细胞类型中三个基因的平均表达
cell_type_avg <- colMeans(avg_expression)

# 创建一个数据框来存储结果
result_df <- data.frame(
  CellType = names(cell_type_avg),
    AvgExpression = cell_type_avg
  )

# 对平均表达值进行排序（可选）
result_df <- result_df[order(result_df$AvgExpression, decreasing = TRUE), ]

# 准备热图数据
heatmap_data <- matrix(result_df$AvgExpression, ncol = 1)
rownames(heatmap_data) <- result_df$CellType
colnames(heatmap_data) <- "AvgExpression"

# 绘制热图
pdf(paste0(outdir,"/Mmu","_AvgExpression.pdf"))
#pheatmap(heatmap_data,
#          cluster_rows = FALSE,  # 不对行进行聚类
#           cluster_cols = FALSE,  # 不对列进行聚类
#            show_colnames = FALSE, # 不显示列名
#         main = "Average Expression of Selected Genes Across Cell Types",
#          fontsize = 10,
#           cellwidth = 25,
#            cellheight = 15,
#         color = colorRampPalette(c("blue", "white", "red"))(100))
color_palette <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
p=pheatmap(heatmap_data,
   fontsize = 10,
     fontsize_row = 10,
   fontsize_col = 10,
     cellwidth = 30,
   cellheight = 20,
     color = color_palette,
   cluster_rows = FALSE,
     cluster_cols = FALSE,
   show_colnames = FALSE,
     angle_col = 0,
   main = "",
     border_color = NA,
   legend = TRUE,
     legend_breaks = seq(min(heatmap_data), max(heatmap_data), length.out = 3),
   legend_labels = c("Low", "Medium", "High")
   )
print(p)
dev.off()

########## 显著性统计 
genes_of_interest=genelist
genes_of_interest=intersect(genes_of_interest,rownames(seurat_obj@assays$RNA@data))
seurat_obj$avg_expression <- colMeans(seurat_obj@assays$RNA@data[genes_of_interest, ])
## 2. 计算每个细胞平均表达
data = as.data.frame(cbind(seurat_obj$papercelltype,seurat_obj$avg_expression))
colnames(data)=c("papercelltype","avg_expression")
## (使用Wilcoxon秩和检验)
# 首先，确保你的数据是按照papercelltype分组的
data$papercelltype <- as.factor(data$papercelltype)
# 获取除了"WIF"之外的所有细胞类型
other_cell_types <- setdiff(levels(data$papercelltype), "WIF")
# 创建一个空列表来存储所有的检验结果
test_results <- list()
# 对"WIF"细胞类型与其它每个细胞类型进行Wilcoxon秩和检验
for(ct in other_cell_types) {
  # 提取"WIF"细胞类型的表达数据
  wif_expression <- as.numeric(data$avg_expression[data$papercelltype == "WIF"])
  # 提取当前细胞类型的表达数据
  other_expression <- as.numeric(data$avg_expression[data$papercelltype == ct])
    # 执行Wilcoxon秩和检验
    test_result <- wilcox.test(wif_expression, other_expression, alternative = "two.sided")
    ## 计算 LogFC
    avg_LogFC <- log2(mean(wif_expression / mean(other_expression)))
  # 保存检验结果
  test_results[[paste("WIF_vs_", ct, sep = "")]] <- list(
     statistic = test_result$statistic,
     p.value = test_result$p.value,
     avg_LogFC = avg_LogFC
   )
}
# 将检验结果转换为数据框以便查看
test_results_df <- do.call(rbind, test_results) %>% as.data.frame()
names(test_results_df) <- c("Statistic", "P.Value","avg_LogFC")
# 查看结果
print(test_results_df)
# 手动构建数据框
Statistic <- sapply(test_results, function(x) x$statistic)
P_Value <- sapply(test_results, function(x) x$p.value)
avg_LogFC <- sapply(test_results, function(x) x$avg_LogFC)

test_results_df <- data.frame(
    Cell_Type = names(test_results),
  Statistic = Statistic,
    P.Value = P_Value,
  avg_LogFC = avg_LogFC
)
write.table(test_results_df,paste0(outdir,"/Mmu","_wilcox.txt"),quote=F,row.names=T,col.names=T,sep="\t")
######################### vlnplot
my_comparisons=list(c("WIF","Fibroblast"),c("WIF","Immune"),c("WIF","Melanocytes"),c("WIF","Keratinocyte"),c("WIF","Muscle"),c("WIF","Chondrocyte"),c("WIF","Endothelium"))
plotfile=cbind(seurat_obj$papercelltype,seurat_obj$avg_expression)
colnames(plotfile)=c("celltype","score")
plotfile=as.data.frame(plotfile)
plotfile$score=as.numeric(plotfile$score)
plotfile <- plotfile %>%
  group_by(celltype) %>%
  mutate(mean_score = mean(score, na.rm = TRUE)) %>%
  ungroup()
# 使用reorder()根据平均score对time进行排序，从大到小
plotfile$celltype <- reorder(plotfile$celltype, plotfile$mean_score, FUN = max, na.rm = TRUE)
pdf(paste0(outdir,"/","Mmu_VlnPlot.pdf"))
p=ggplot(data = plotfile, aes(x = celltype, y = score, color = celltype)) +
  geom_violin() +
    geom_boxplot(width = 0.1) +
  scale_color_brewer(palette = "Set1") + 
    theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = .9, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
  scale_y_continuous(limits = c(0, max(plotfile$score, na.rm = TRUE) + 2)) +
    ylab("AvgExp")
  print(p)
dev.off()


# 如果您想保存热图
# ggsave("avg_gene_expression_heatmap.png", width = 6, height = 8)

#################################### rb (only rabbit)#######################
# 使用AverageExpression获取每个细胞类型的平均表达
#seurat_obj=sc[,sc$species=="Rabbit"]
seurat_obj=sc[,sc$sptime=="RbD10"]
genelist=intersect(processed_genelist,rownames(seurat_obj@assays$RNA@counts))
avg_expression <- AverageExpression(seurat_obj, features = genelist, group.by = "ident")

# 提取RNA assay的结果（如果您使用的是默认的RNA assay）
avg_expression <- avg_expression$RNA

# 计算每个细胞类型中三个基因的平均表达
cell_type_avg <- colMeans(avg_expression)

# 创建一个数据框来存储结果
result_df <- data.frame(
  CellType = names(cell_type_avg),
    AvgExpression = cell_type_avg
  )

# 对平均表达值进行排序（可选）
result_df <- result_df[order(result_df$AvgExpression, decreasing = TRUE), ]

# 准备热图数据
heatmap_data <- matrix(result_df$AvgExpression, ncol = 1)
rownames(heatmap_data) <- result_df$CellType
colnames(heatmap_data) <- "AvgExpression"

# 绘制热图
pdf(paste0(outdir,"/Rb","_AvgExpression.pdf"))
#pheatmap(heatmap_data,
#          cluster_rows = FALSE,  # 不对行进行聚类
#           cluster_cols = FALSE,  # 不对列进行聚类
#            show_colnames = FALSE, # 不显示列名
#         main = "Average Expression of Selected Genes Across Cell Types",
#          fontsize = 10,
#           cellwidth = 25,
#            cellheight = 15,
#         color = colorRampPalette(c("blue", "white", "red"))(100))
color_palette <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
p=pheatmap(heatmap_data,
   fontsize = 10,
     fontsize_row = 10,
   fontsize_col = 10,
     cellwidth = 30,
   cellheight = 20,
     color = color_palette,
   cluster_rows = FALSE,
     cluster_cols = FALSE,
   show_colnames = FALSE,
     angle_col = 0,
   main = "",
     border_color = NA,
   legend = TRUE,
     legend_breaks = seq(min(heatmap_data), max(heatmap_data), length.out = 3),
   legend_labels = c("Low", "Medium", "High")
   )
print(p)
dev.off()

########## 显著性统计 
genes_of_interest=genelist
genes_of_interest=intersect(genes_of_interest,rownames(seurat_obj@assays$RNA@data))
seurat_obj$avg_expression <- colMeans(seurat_obj@assays$RNA@data[genes_of_interest, ])
## 2. 计算每个细胞平均表达
data = as.data.frame(cbind(seurat_obj$papercelltype,seurat_obj$avg_expression))
colnames(data)=c("papercelltype","avg_expression")
## (使用Wilcoxon秩和检验)
# 首先，确保你的数据是按照papercelltype分组的
data$papercelltype <- as.factor(data$papercelltype)
# 获取除了"WIF"之外的所有细胞类型
other_cell_types <- setdiff(levels(data$papercelltype), "WIF")
# 创建一个空列表来存储所有的检验结果
test_results <- list()
# 对"WIF"细胞类型与其它每个细胞类型进行Wilcoxon秩和检验
for(ct in other_cell_types) {
  # 提取"WIF"细胞类型的表达数据
  wif_expression <- as.numeric(data$avg_expression[data$papercelltype == "WIF"])
  # 提取当前细胞类型的表达数据
  other_expression <- as.numeric(data$avg_expression[data$papercelltype == ct])
    # 执行Wilcoxon秩和检验
    test_result <- wilcox.test(wif_expression, other_expression, alternative = "two.sided")
    ## 计算 LogFC
    avg_LogFC <- log2(mean(wif_expression / mean(other_expression)))
  # 保存检验结果
  test_results[[paste("WIF_vs_", ct, sep = "")]] <- list(
     statistic = test_result$statistic,
     p.value = test_result$p.value,
     avg_LogFC = avg_LogFC
   )
}
# 将检验结果转换为数据框以便查看
test_results_df <- do.call(rbind, test_results) %>% as.data.frame()
names(test_results_df) <- c("Statistic", "P.Value","avg_LogFC")
# 查看结果
print(test_results_df)
# 手动构建数据框
Statistic <- sapply(test_results, function(x) x$statistic)
P_Value <- sapply(test_results, function(x) x$p.value)
avg_LogFC <- sapply(test_results, function(x) x$avg_LogFC)

test_results_df <- data.frame(
    Cell_Type = names(test_results),
  Statistic = Statistic,
    P.Value = P_Value,
  avg_LogFC = avg_LogFC
)
write.table(test_results_df,paste0(outdir,"/Rb","_wilcox.txt"),quote=F,row.names=T,col.names=T,sep="\t")
######################### vlnplot
my_comparisons=list(c("WIF","Fibroblast"),c("WIF","Immune"),c("WIF","Melanocytes"),c("WIF","Keratinocyte"),c("WIF","Muscle"),c("WIF","Chondrocyte"),c("WIF","Endothelium"))
plotfile=cbind(seurat_obj$papercelltype,seurat_obj$avg_expression)
colnames(plotfile)=c("celltype","score")
plotfile=as.data.frame(plotfile)
plotfile$score=as.numeric(plotfile$score)
plotfile <- plotfile %>%
  group_by(celltype) %>%
  mutate(mean_score = mean(score, na.rm = TRUE)) %>%
  ungroup()
# 使用reorder()根据平均score对time进行排序，从大到小
plotfile$celltype <- reorder(plotfile$celltype, plotfile$mean_score, FUN = max, na.rm = TRUE)
pdf(paste0(outdir,"/","Rb_VlnPlot.pdf"))
p=ggplot(data = plotfile, aes(x = celltype, y = score, color = celltype)) +
  geom_violin() +
    geom_boxplot(width = 0.1) +
  scale_color_brewer(palette = "Set1") + 
    theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = .9, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
  scale_y_continuous(limits = c(0, max(plotfile$score, na.rm = TRUE) + 1.5)) +
    ylab("AvgExp")
print(p)
dev.off()


# 如果您想保存热图
# ggsave("avg_gene_expression_heatmap.png", width = 6, height = 8)


