coolpup.py data/CLP.cool data/CTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname CLP.CTCF.noTCF1.all_genes.txt
coolpup.py data/DN2.cool data/CTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN2.CTCF.noTCF1.all_genes.txt
coolpup.py data/DN3.cool data/CTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN3.CTCF.noTCF1.all_genes.txt
coolpup.py data/DN4.cool data/CTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN4.CTCF.noTCF1.all_genes.txt
coolpup.py data/DP.cool data/CTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DP.CTCF.noTCF1.all_genes.txt
coolpup.py data/ETP.cool data/CTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname ETP.CTCF.noTCF1.all_genes.txt
coolpup.py data/CLP.cool data/TCF1.CTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname CLP.TCF1.CTCF.all_genes.txt
coolpup.py data/DN2.cool data/TCF1.CTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN2.TCF1.CTCF.all_genes.txt
coolpup.py data/DN3.cool data/TCF1.CTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN3.TCF1.CTCF.all_genes.txt
coolpup.py data/DN4.cool data/TCF1.CTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN4.TCF1.CTCF.all_genes.txt
coolpup.py data/DP.cool data/TCF1.CTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DP.TCF1.CTCF.all_genes.txt
coolpup.py data/ETP.cool data/TCF1.CTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname ETP.TCF1.CTCF.all_genes.txt
coolpup.py data/CLP.cool data/TCF1.noCTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname CLP.TCF1.noCTCF.all_genes.txt
coolpup.py data/DN2.cool data/TCF1.noCTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN2.TCF1.noCTCF.all_genes.txt
coolpup.py data/DN3.cool data/TCF1.noCTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN3.TCF1.noCTCF.all_genes.txt
coolpup.py data/DN4.cool data/TCF1.noCTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN4.TCF1.noCTCF.all_genes.txt
coolpup.py data/DP.cool data/TCF1.noCTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DP.TCF1.noCTCF.all_genes.txt
coolpup.py data/ETP.cool data/TCF1.noCTCF.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname ETP.TCF1.noCTCF.all_genes.txt
coolpup.py data/CLP.cool data/noCTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname CLP.noCTCF.noTCF1.all_genes.txt
coolpup.py data/DN2.cool data/noCTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN2.noCTCF.noTCF1.all_genes.txt
coolpup.py data/DN3.cool data/noCTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN3.noCTCF.noTCF1.all_genes.txt
coolpup.py data/DN4.cool data/noCTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DN4.noCTCF.noTCF1.all_genes.txt
coolpup.py data/DP.cool data/noCTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname DP.noCTCF.noTCF1.all_genes.txt
coolpup.py data/ETP.cool data/noCTCF.noTCF1.bed --bed2 data/all_genes.bed --mindist 200000 --maxdist 2000000 --outname ETP.noCTCF.noTCF1.all_genes.txt

plotpup.py CLP.*.txt ETP.*.txt D*.txt --row_names CLP,ETP,DN2,DN3,DN4,DP --col_names CTCF,TCF1.CTCF,TCF1,None --n_cols 4 --output peak_genes.pileup.pdf 
