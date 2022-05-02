data=read.table("data/A.scatterPC1.tab", header=F)
ggplot(data, aes(V4, V5)) + geom_point(alpha=0.1) + xlab("WT PC1") + ylab("KO PC1") + geom_abline(intercept = 0, slope = 1, color='red') + geom_hline(yintercept = 0, color='blue') + geom_vline(xintercept = 0, color='blue') + stat_cor() + scale_x_continuous(limits = c(-3,3)) + scale_y_continuous(limits = c(-3, 3)) + theme_classic() + ggsave("Figure4A.pdf", width = 4, height = 4)

