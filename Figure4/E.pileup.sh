python3 pileup.py --cool data/DP.cool data/DP_TCF1_KO.cool --bed data/DP.TCF1.withCTCF.bed data/DP.TCF1.noCTCF.bed --bname Both.DP TCF1.DP --outdir data/

awk '$4-$5>0.2 {print $1,$2,$3}' OFS="\t" data/Both.DP.enrich.csv > data/Both.decrease0.2.bed
awk '$4-$5>0.2 {print $1,$2,$3}' OFS="\t" data/TCF1.DP.enrich.csv > data/TCF1.decrease0.2.bed

coolpup.py ../data/hic/DP_10k.cool data/Both.decrease0.05.bed --outname DP_WT.Both.decrease0.05.txt --pad 250 --local --ignore_diags 0
coolpup.py ../data/hic/DP_TCF1_KO_10k.cool data/Both.decrease0.05.bed --outname DP_KO.Both.decrease0.05.txt --pad 250 --local --ignore_diags 0
coolpup.py ../data/hic/DP_10k.cool data/TCF1.decrease0.05.bed --outname DP_WT.TCF1.decrease0.05.txt --pad 250 --local --ignore_diags 0
coolpup.py ../data/hic/DP_TCF1_KO_10k.cool data/TCF1.decrease0.05.bed --outname DP_KO.TCF1.decrease0.05.txt --pad 250 --local --ignore_diags 0

plotpup.py DP_WT.*0.05.txt DP_KO.*0.05.txt --row_names WT,KO --col_names Both,TCF1 --enrichment 1  --output Figure5C_0.05_10k.pdf
