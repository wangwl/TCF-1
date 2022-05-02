computeMatrix reference-point -R data/CTCF.noTCF1.bed data/TCF1.withCTCF.bed data/TCF1.noCTCF.bed data/noCTCF.noTCF1.bed -S mm10.60way.phastCons.bw -a 2500 -b 2500 --skipZeros --missingDataAsZero -o all_phastCons.gz --referencePoint center -p 30

plotHeatmap -m all_phastCons.gz --colorMap YlGnBu --perGroup --heatmapHeight 10

