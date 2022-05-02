data=read.table("data/Domain_score_binding.tab", header = T)

ggplot(data, aes(Change, TCF1)) + geom_boxplot(aes(group=cut_width(Change, 0.005)), outlier.size = 0.3) + stat_cor() + xlab("Domain Score Change") + ylab("TCF-1 Peak Density (per Mb)") + theme_classic()  + ggsave("FigureS5H.pdf")
