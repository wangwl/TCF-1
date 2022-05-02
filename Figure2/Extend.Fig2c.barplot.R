data=read.table("data/barplot.tab", header = F)
ggplot(data, aes(V1, V2)) + geom_bar(stat="identity") + xlab("") + ylab("Number of Boundaries") + scale_x_discrete(limits=c("CLP", "ETP", "DN2", "DN3", "DN4", "DP")) + theme_classic() + ggsave("FigureS2C.pdf", width = 6, height = 6)
