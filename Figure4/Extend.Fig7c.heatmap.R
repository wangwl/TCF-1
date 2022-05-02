data=read.table("data/NOD_BL6.Ifi.csv", header=F)
rownames(data) = data$V1
data=data[,-1]
pheatmap(data, show_colnames = F, cluster_cols = F, scale = 'row')
