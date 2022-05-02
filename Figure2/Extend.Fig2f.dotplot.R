data = read.table("data/upright_enrich.tab", header=F)
ggplot(data, aes(V2, V3, color=V1)) + geom_point() + ggsave("FigureS2F.pdf")
