library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)

bf_folder = "./output"

countFiles = BamFileList(dir(bf_folder, pattern = "ReadsPerGene.out.tab$", full = TRUE))
toProcess = data.frame(name=names(countFiles), file=path(countFiles), stringsAsFactors = FALSE)
toProcess = mutate(toProcess, status=gsub("ReadsPerGene.out.tab","", toProcess$name) )

data = toProcess$file
sample_name = toProcess$status

merge.data = read.table(data[1], header = FALSE, sep = '\t')
merge.data = tail(merge.data, -4)
merge.data = merge.data[,c(1,2)]

colnames(merge.data) = c("GeneID", sample_name[1])

for(i in 2:length(sample_name)){
  new_data = read.table(data[i], header = FALSE, sep = '\t')
  new_data = tail(new_data, -4)
  new_data = new_data[,c(1,2)]
  colnames(new_data) = c("GeneID", sample_name[i])
  merge.data = full_join(merge.data, new_data, by = "GeneID")
}

write.csv(merge.data, file = "output_read_counts.wanglab.csv", row.names = F)