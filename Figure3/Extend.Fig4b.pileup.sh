coolpup.py data/ETP_TALL1.cool data/CTCF.noTCF1.mid.bed --outname ETP_TALL1.CTCF.noTCF1.mid.txt --pad 250 --local
coolpup.py data/HPAP027.cool data/CTCF.noTCF1.mid.bed --outname HPAP027.CTCF.noTCF1.mid.txt --pad 250 --local
coolpup.py data/TALL1.cool data/CTCF.noTCF1.mid.bed --outname TALL1.CTCF.noTCF1.mid.txt --pad 250 --local
coolpup.py data/ETP_TALL1.cool data/TCF1.noCTCF.mid.bed --outname ETP_TALL1.TCF1.noCTCF.mid.txt --pad 250 --local
coolpup.py data/HPAP027.cool data/TCF1.noCTCF.mid.bed --outname HPAP027.TCF1.noCTCF.mid.txt --pad 250 --local
coolpup.py data/TALL1.cool data/TCF1.noCTCF.mid.bed --outname TALL1.TCF1.noCTCF.mid.txt --pad 250 --local
coolpup.py data/ETP_TALL1.cool data/TCF1.withCTCF.mid.bed --outname ETP_TALL1.TCF1.withCTCF.mid.txt --pad 250 --local
coolpup.py data/HPAP027.cool data/TCF1.withCTCF.mid.bed --outname HPAP027.TCF1.withCTCF.mid.txt --pad 250 --local
coolpup.py data/TALL1.cool data/TCF1.withCTCF.mid.bed --outname TALL1.TCF1.withCTCF.mid.txt --pad 250 --local

plotpup.py HPAP*.txt TALL1.*.txt ETP_TALL*.txt --row_names HPAP,TALL,ETP_TALL --col_names CTCF,TCF1,cobind --enrichment 0 --output Figure6B.pdf --n_cols 3
