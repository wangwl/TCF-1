python pileup.py --bed data/NIH3T3.TCF1.withCTCF.bed data/NIH3T3.TCF1.noCTCF.bed --bname Both.3T3 TCF1.3T3 --outdir data/ --cool data/Norm_3T3.cool data/TCF1_3T3.cool

awk '$5-$4>0.05 {print $1,$2,$3}' OFS="\t" data/Both.3T3.enrich.csv > data/Both.3T3.Increase0.05.bed
awk '$5-$4>0.05 {print $1,$2,$3}' OFS="\t" data/TCF1.3T3.enrich.csv > data/TCF1.3T3.Increase0.05.bed

coolpup.py data/TCF1_3T3.cool data/Both.3T3.Increase0.05.bed --outname TCF1_3T3.Both.Increase0.05.txt --pad 250 --local
coolpup.py data/Norm_3T3.cool data/Both.3T3.Increase0.05.bed --outname Norm_3T3.Both.Increase0.05.txt --pad 250 --local
coolpup.py data/TCF1_3T3.cool data/TCF1.3T3.Increase0.05.bed --outname TCF1_3T3.TCF1.Increase0.05.txt --pad 250 --local
coolpup.py data/Norm_3T3.cool data/TCF1.3T3.Increase0.05.bed --outname Norm_3T3.TCF1.Increase0.05.txt --pad 250 --local

plotpup.py TCF1_3T3*.txt Norm_3T3.*.txt --row_names TCF1,Norm --col_names Both,TCF1 --enrichment 0 --output Figure5A.pdf

