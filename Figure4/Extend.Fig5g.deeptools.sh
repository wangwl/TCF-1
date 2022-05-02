computeMatrix scale-regions -R data/DN3EV.bed data/common.bed data/DN3KO.bed -S data/TCF1.bw --regionBodyLength 100000 -bs 10000 -a 400000 -b 400000 --skipZeros --missingDataAsZero -o uniq.bound.gz -p 30
plotHeatmap -m uniq.bound.gz -o FigureS7C.pdf --colorMap YlGnBu --perGroup --heatmapHeight 10
