data=read.table("data/S7E.volcano.tab", header=F)
ggplot(data, aes(V2, -log10(V3))) + geom_point(alpha=0.5) + geom_text(data = subset(data, V3<0.05 & abs(V2)>1), aes(label=V1), color='red') + xlab("log2 Fold Change (WT/KO)") + ylab("-log10(Pvalue)")  + theme_classic() + ggsave("FigureS7E.pdf", width = 4, height = 4)

