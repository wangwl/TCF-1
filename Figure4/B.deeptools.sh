computeMatrix reference-point -R data/Both.Increase0.05.bed data/TCF1.Increase0.05.bed -S data/Norm_3T3.insulation.bw data/TCF1_3T3.insulation.bw --samplesLabel Norm TCF1 -a 400000 -b 400000 -bs 1000 --skipZeros --missingDataAsZero -o all.IS.gz --referencePoint center -p 30
plotHeatmap -m all.IS.gz --colorMap PiYG --perGroup --heatmapHeight 10 -o Figure4B.pdf
