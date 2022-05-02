computeMatrix reference-point -R data/Both.decrease0.05.bed data/TCF1.decrease0.05.bed -S data/DN3EV.IS.bw data/DN3KO.IS.bw --samplesLabel WT KO -a 400000 -b 400000 -bs 1000 --skipZeros --missingDataAsZero -o DN3.IS.gz --referencePoint center -p 30
plotHeatmap -m DN3.IS.gz --colorMap PiYG --perGroup --heatmapHeight 10 -o Figure4D.pdf
