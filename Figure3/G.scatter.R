data=read.table("data/G.scatter.tab", header=T)
ggplot(data, aes(DN3EV, DN3KO)) + geom_point(alpha=0.3) + ylab("TCF-1 KO") + xlab("WT") + geom_abline(slope = 1, intercept = 0, color='red') + theme_classic()
ggsave("Figure3G.pdf", width = 7, height = 7)

