python pileup.py --cool data/DN3EV.cool data/DN3KO.cool --bed data/DP.TCF1.withCTCF.bed data/DP.TCF1.noCTCF.bed --bname Both.DN3 TCF1.DN3 --outdir data/

awk '$4-$5>0.05 {print $1,$2,$3}' OFS="\t" data/Both.DN3.enrich.csv > data/Both.decrease0.05.bed
awk '$4-$5>0.05 {print $1,$2,$3}' OFS="\t" data/TCF1.DN3.enrich.csv > data/TCF1.decrease0.05.bed

coolpup.py data/DN3EV.cool data/Both.decrease0.05.bed --outname DN3EV.Both.decrease0.05.txt --pad 250 --local
coolpup.py data/DN3KO.cool data/Both.decrease0.05.bed --outname DN3KO.Both.decrease0.05.txt --pad 250 --local
coolpup.py data/DN3EV.cool data/TCF1.decrease0.05.bed --outname DN3EV.TCF1.decrease0.05.txt --pad 250 --local
coolpup.py data/DN3KO.cool data/TCF1.decrease0.05.bed --outname DN3KO.TCF1.decrease0.05.txt --pad 250 --local

plotpup.py DN3EV*.txt DN3KO.*.txt --row_names WT,KO --col_names Both,TCF1 --enrichment 0 --output Figure4C.pdf
