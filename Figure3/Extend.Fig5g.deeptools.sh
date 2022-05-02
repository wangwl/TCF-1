computeMatrix scale-regions -R data/Norm.bed data/common.bed data/TCF1.bed -S data/TCF1.3T3.bw --regionBodyLength 100000 -bs 10000 -a 400000 -b 400000 --skipZeros --missingDataAsZero -o uniq.bound.gz -p 30
plotHeatmap -m uniq.bound.gz -o FigureS5E.pdf --colorMap YlGnBu --perGroup --heatmapHeight 10
