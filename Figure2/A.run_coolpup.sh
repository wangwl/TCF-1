coolpup.py data/CLP.cool data/CTCF.noTCF1.bed --outname CLP.CTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DN2.cool data/CTCF.noTCF1.bed --outname DN2.CTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DN3.cool data/CTCF.noTCF1.bed --outname DN3.CTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DN4.cool data/CTCF.noTCF1.bed --outname DN4.CTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DP.cool data/CTCF.noTCF1.bed --outname DP.CTCF.noTCF1.txt --pad 250 --local
coolpup.py data/ETP.cool data/CTCF.noTCF1.bed --outname ETP.CTCF.noTCF1.txt --pad 250 --local
coolpup.py data/CLP.cool data/TCF1.CTCF.bed --outname CLP.TCF1.CTCF.txt --pad 250 --local
coolpup.py data/DN2.cool data/TCF1.CTCF.bed --outname DN2.TCF1.CTCF.txt --pad 250 --local
coolpup.py data/DN3.cool data/TCF1.CTCF.bed --outname DN3.TCF1.CTCF.txt --pad 250 --local
coolpup.py data/DN4.cool data/TCF1.CTCF.bed --outname DN4.TCF1.CTCF.txt --pad 250 --local
coolpup.py data/DP.cool data/TCF1.CTCF.bed --outname DP.TCF1.CTCF.txt --pad 250 --local
coolpup.py data/ETP.cool data/TCF1.CTCF.bed --outname ETP.TCF1.CTCF.txt --pad 250 --local
coolpup.py data/CLP.cool data/TCF1.noCTCF.bed --outname CLP.TCF1.noCTCF.txt --pad 250 --local
coolpup.py data/DN2.cool data/TCF1.noCTCF.bed --outname DN2.TCF1.noCTCF.txt --pad 250 --local
coolpup.py data/DN3.cool data/TCF1.noCTCF.bed --outname DN3.TCF1.noCTCF.txt --pad 250 --local
coolpup.py data/DN4.cool data/TCF1.noCTCF.bed --outname DN4.TCF1.noCTCF.txt --pad 250 --local
coolpup.py data/DP.cool data/TCF1.noCTCF.bed --outname DP.TCF1.noCTCF.txt --pad 250 --local
coolpup.py data/ETP.cool data/TCF1.noCTCF.bed --outname ETP.TCF1.noCTCF.txt --pad 250 --local
coolpup.py data/CLP.cool data/noCTCF.noTCF1.bed --outname CLP.noCTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DN2.cool data/noCTCF.noTCF1.bed --outname DN2.noCTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DN3.cool data/noCTCF.noTCF1.bed --outname DN3.noCTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DN4.cool data/noCTCF.noTCF1.bed --outname DN4.noCTCF.noTCF1.txt --pad 250 --local
coolpup.py data/DP.cool data/noCTCF.noTCF1.bed --outname DP.noCTCF.noTCF1.txt --pad 250 --local
coolpup.py data/ETP.cool data/noCTCF.noTCF1.bed --outname ETP.noCTCF.noTCF1.txt --pad 250 --local

plotpup.py CLP.*.txt ETP.*.txt D*.txt --row_names CLP,ETP,DN2,DN3,DN4,DP --col_names CTCF,TCF1.CTCF,TCF1,None --n_cols 4 --output Figure2A.pileup.pdf --enrichment 0
