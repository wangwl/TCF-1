library(ggplot2)
library(ggpubr)
my_comparisons <- list( c("CTCF+TCF-1+", "CTCF+TCF-1-"), c("CTCF+TCF-1+", "CTCF-TCF-1+"), c("CTCF+TCF-1+", "CTCF-TCF-1-"), c("CTCF-TCF-1+", "CTCF-TCF-1-"), c("CTCF-TCF-1+", "CTCF+TCF-1-") )
data=read.csv("data/TFs_insulation.20K.csv")
data[data$Thy_TCF1==1 & data$DP_CTCF==1,]$X = "CTCF+TCF-1+"
data[data$Thy_TCF1==0 & data$DP_CTCF==1,]$X = "CTCF+TCF-1-"
data[data$Thy_TCF1==0 & data$DP_CTCF==0,]$X = "CTCF-TCF-1-"
data[data$Thy_TCF1==1 & data$DP_CTCF==0,]$X = "CTCF-TCF-1+"
data = data[,c("X", "CLP", "ETP", "DN2", "DN3", "DN4", "DP")]
#data_melt = melt(data, id.vars = c("X", "DP"))
data$diff = data$DP - data$CLP
nadata = data[!is.na(data$diff),]
nadata[nadata$diff>2.5,]$diff = 2.5
nadata[nadata$diff <= -2,]$diff = -2
ggplot(nadata, aes(X, diff)) + geom_violin() + geom_boxplot(width=0.1) +  stat_compare_means(comparisons = my_comparisons,  method = 't.test', method.args = list(alternative = "greater", exact = TRUE)) + scale_x_discrete(limits=c("CTCF+TCF-1-", "CTCF+TCF-1+", "CTCF-TCF-1+", "CTCF-TCF-1-")) + xlab("") + ylab("Increase of Insulation Score") + theme_classic() + ggsave("FigureS2B.pdf", width = 4, height = 4)

