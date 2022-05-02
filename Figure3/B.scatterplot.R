data=read.table("data/Domain_score_binding.tab", header = T)

ggplot(data, aes(Norm_3T3, TCF1_3T3)) + geom_point(alpha=0.5) + geom_abline(slope = 1, intercept = 0, color='red') + xlab("Norm 3T3") + ylab("TCF-1 3T3") + ggtitle("Domain Score") + theme_classic() + ggsave("Figure3B.pdf")

