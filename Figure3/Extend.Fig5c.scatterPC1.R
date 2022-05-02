data=read.table("data/S5A.scatterPC1.bedgraph", header=F)
ggplot(data, aes(V4, V5)) + geom_point(alpha=0.1) + xlab("Norm 3T3") + ylab("TCF-1 3T3") + geom_abline(intercept = 0, slope = 1, color='red') + geom_hline(yintercept = 0, color='blue') + geom_vline(xintercept = 0, color='blue') + stat_cor() + scale_x_continuous(limits = c(-3,3)) + scale_y_continuous(limits = c(-3,3)) + theme_classic() + ggsave("FigureS5A.pdf", height = 4, width = 4)

