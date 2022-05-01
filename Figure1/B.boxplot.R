data=read.table("data/Domain_PETs.binding.tab", header = T)
data[data$Peak>150,]$Peak = 150
data_melt = melt(data, id.vars = c("Chrom", "Start", "End", "DP", "Peak"))
ggplot(data_melt, aes(Peak, log2(DP/value))) + geom_boxplot(aes(group=cut_width(Peak, 2.5)), outlier.size = 0.3) + facet_grid(variable~., scales = 'free') + ylab("log2 Fold Change in DP T cells") + xlab("TCF-1 Peak density") + theme_classic() + ggsave("Figure1B.boxplot.pdf", width = 6, height = 8)
