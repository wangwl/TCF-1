data=read.table("data/G.diff_bins.TCF1_counts.bed", header=F)
ggplot(data, aes(as.factor(V8), V6, fill=V7)) + geom_boxplot(outlier.size = 0.3) + xlab("Number of TCF1 peaks within 50 Kb") + ylab("Log2 Fold Change (TCF1 induced 3T3/untreated 3T3)") + ggsave("Figure3G.pdf")
