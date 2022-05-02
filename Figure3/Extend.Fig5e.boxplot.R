data=read.table("data/S5C.boxplot.tab", header = T)
data_melt=melt(data, id.vars = c("Group"))
pdf("FigureS5C.pdf", useDingbats = F, width = 8, height = 6)
ggplot(data_melt, aes(Group, value, fill=Group)) + geom_boxplot() + facet_wrap(~variable, nrow = 1, ncol = 4, scales = 'free') + theme_classic()
dev.off()

