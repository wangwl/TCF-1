coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/Bcell1.mcool::/resolutions/5000 data/CTCF.noTCF1.bed --outname Bcell1.CTCF.noTCF1.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/CD4.mcool::/resolutions/5000 data/CTCF.noTCF1.bed --outname CD4.CTCF.noTCF1.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/Cardio.mcool::/resolutions/5000 data/CTCF.noTCF1.bed --outname Cardio.CTCF.noTCF1.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/DP.mcool::/resolutions/5000 data/CTCF.noTCF1.bed --outname DP.CTCF.noTCF1.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/NIH3T3.mcool::/resolutions/5000 data/CTCF.noTCF1.bed --outname NIH3T3.CTCF.noTCF1.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/mESC.mcool::/resolutions/5000 data/CTCF.noTCF1.bed --outname mESC.CTCF.noTCF1.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/Bcell1.mcool::/resolutions/5000 data/TCF1.noCTCF.bed --outname Bcell1.TCF1.noCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/CD4.mcool::/resolutions/5000 data/TCF1.noCTCF.bed --outname CD4.TCF1.noCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/Cardio.mcool::/resolutions/5000 data/TCF1.noCTCF.bed --outname Cardio.TCF1.noCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/DP.mcool::/resolutions/5000 data/TCF1.noCTCF.bed --outname DP.TCF1.noCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/NIH3T3.mcool::/resolutions/5000 data/TCF1.noCTCF.bed --outname NIH3T3.TCF1.noCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/mESC.mcool::/resolutions/5000 data/TCF1.noCTCF.bed --outname mESC.TCF1.noCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/Bcell1.mcool::/resolutions/5000 data/TCF1.withCTCF.bed --outname Bcell1.TCF1.withCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/CD4.mcool::/resolutions/5000 data/TCF1.withCTCF.bed --outname CD4.TCF1.withCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/Cardio.mcool::/resolutions/5000 data/TCF1.withCTCF.bed --outname Cardio.TCF1.withCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/DP.mcool::/resolutions/5000 data/TCF1.withCTCF.bed --outname DP.TCF1.withCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/NIH3T3.mcool::/resolutions/5000 data/TCF1.withCTCF.bed --outname NIH3T3.TCF1.withCTCF.txt --pad 250 --local
coolpup.py /mnt/data0/wenliang/collaborate/Golnaz/Tcellpileup/01.mouse/00.cool/mESC.mcool::/resolutions/5000 data/TCF1.withCTCF.bed --outname mESC.TCF1.withCTCF.txt --pad 250 --local

coolpup.py mESC.*.txt NIH3T3.*.txt Cardio.*.txt Bcell1.*.txt CD4*.txt DP*.txt --row_names mESC,NIH3T3,Cardio,Bcell,CD4,DP --col_names CTCF,TCF1,cobind --n_cols 3 --enrichment 0 --output Figure6C.pdf
