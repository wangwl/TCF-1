data=read.table("data/NIH3T3.DEseq2.csv", header=F)
rownames(data) = data$V1
data=data[,-1]
pheatmap(data, show_colnames = F, cluster_cols = F, scale = 'row')
