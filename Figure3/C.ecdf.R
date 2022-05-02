data=read.table("data/E.CTCF_TCF1.bed", header=F)
ggplot(data, aes(V4, color=paste(V5,V7))) + stat_ecdf() + ggsave("Figure3E.pdf")

