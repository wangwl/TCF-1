my_comparisons <- list( c("CTCF+TCF-1+", "CTCF+TCF-1-"), c("CTCF+TCF-1+", "CTCF-TCF-1+"), c("CTCF+TCF-1+", "CTCF-TCF-1-"), c("CTCF-TCF-1+", "CTCF-TCF-1-"), c("CTCF-TCF-1+", "CTCF+TCF-1-") )
data=read.table("data/S5G.violin.tab", header = F)
ggplot(data, aes(V6, V4-V5)) + geom_violin() + geom_boxplot(width=0.2) + stat_compare_means(comparisons = my_comparisons, method = 't.test', method.args = list(alternative = "less")) +stat_compare_means() + scale_x_discrete(limits=c("CTCF+TCF-1-", "CTCF+TCF-1+", "CTCF-TCF-1+", "CTCF-TCF-1-")) + xlab("") + ylab("Increase of Loops at Shared Anchors")+ theme_classic() + ggsave("FigureS5G.pdf")

