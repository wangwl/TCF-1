res=read.table("data/CTCF_KO_vs_EV_results.csv", header=T)
pdf("FigureS6B.pdf")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-7,7), ylim=c(0,40), col="darkgray"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="darkgray"))
with(subset(res, abs(log2FoldChange)>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="darkgray"))
with(subset(res, padj<.05 & (log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick"))
with(subset(res, padj<.05 & (log2FoldChange)< -1), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
dev.off()

