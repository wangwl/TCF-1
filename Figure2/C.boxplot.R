data=read.table("data/GeneDist.stat")
my_comparisons <- list( c("Tconv1", "Tconv2"), c("Tconv2", "Tconv3"), c("Tconv3", "Tconv4"), c("Tconv4", "Tconv5"))
ggplot(data, aes(V1, V2)) + geom_boxplot() + xlab("T cell development stages") + ylab("Distance to the TCF-1 and CTCF co-bound peaks") + theme_classic() + stat_compare_means(comparisons = my_comparisons) + ggsave("Figure2C.GeneDist.pdf")

