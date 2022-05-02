# 2021 12 22
# Sushi plot for cd8a region.

library(Sushi)

setwd("/mnt/data0/sora/TCF1_revision/Reproduce/Figure5/")

cd8a = "chr6:71000000-71850000"


WTcool = '/mnt/alvand/sora/TCF1_revision/NovaSeq/Arima-HiC_VavCre_DP_223231022/workdir/HiC_pro/hic_results/data/rawdata/DP.cool'
KOcool = '/mnt/alvand/sora/TCF1_revision/NovaSeq/Arima-HiC_fl_fl_VavCre_DP_223231022/workdir/HiC_pro/hic_results/data/rawdata/TCF1_KO_DP.cool'

geneinfo = read.delim('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed',header=F)
TAD = read.delim('TADs.bed')
# Bcl6
co = cd8a
name="cd8a"
CHR = unlist(strsplit(co, split=":"))[1]
START = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[1])
END = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[2])
OUT = 'DP_WT_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",WTcool," --coord ",co," --norm VC_SQRT --out ", OUT," --res 5000"))
OUT = 'DP_KO_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",KOcool," --coord ",co," --norm VC_SQRT --out ", OUT," --res 5000"))

WT = read.delim(paste0('DP_WT_0.txt'),row.names=1)
colnames(WT) = rownames(WT)
KO = read.delim(paste0('DP_KO_0.txt'),row.names=1)
colnames(KO) = rownames(KO)
#CHR = paste0('chr',CHR)

pdf('Cd8a_R_heatmap4.pdf', width = 10, height = 8)
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.9,0.6,1))
phic_WT = plotHic(WT, CHR ,START,END, max_y = 100,
                   zrange=c(0,4), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.9,0.56,0.6), new=T)
plotBed(beddata = TAD,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)


par(fig=c(0.1,0.9,0.2,0.6), new=T)
phic_KO = plotHic(KO, CHR, START, END, max_y = 100,
                   zrange=c(0,4), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.9,0.16,0.2), new=T)
plotBed(beddata = TAD,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.1,0.9,0.05,0.16), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==CHR &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.0, arrowlength = 0.001, bheight = 0.1, packrow=F)


labelgenome(CHR,START,END,n=5,scale="Mb")

dev.off()


#####################################################################################
# DN3

WTcool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/DN3EV.cool'
KOcool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/DN3KO.cool'
setwd("/mnt/data0/sora/TCF1_revision/Reproduce/Figure5/")

cd8a = "6:71000000-71850000"

geneinfo = read.delim('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed',header=F)
TAD = read.delim('TADs.bed')
# Bcl6
co = cd8a
name="cd8a"
CHR = unlist(strsplit(co, split=":"))[1]
START = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[1])
END = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[2])
OUT = 'DN3_WT_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",WTcool," --coord ",co," --norm VC_SQRT --out ", OUT," --res 5000"))
OUT = 'DN3_KO_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",KOcool," --coord ",co," --norm VC_SQRT --out ", OUT," --res 5000"))

WT = read.delim(paste0('DN3_WT_0.txt'),row.names=1)
colnames(WT) = rownames(WT)
KO = read.delim(paste0('DN3_KO_0.txt'),row.names=1)
colnames(KO) = rownames(KO)
#CHR = paste0('chr',CHR)

pdf('Cd8a_DN3_R_heatmap4.pdf', width = 10, height = 8)
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.9,0.6,1))
phic_WT = plotHic(WT, CHR ,START,END, max_y = 100,
                  zrange=c(0,4.3), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.9,0.56,0.6), new=T)
plotBed(beddata = TAD,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)


par(fig=c(0.1,0.9,0.2,0.6), new=T)
phic_KO = plotHic(KO, CHR, START, END, max_y = 100,
                  zrange=c(0,4.3), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.9,0.16,0.2), new=T)
plotBed(beddata = TAD,chrom = paste0("chr",CHR),chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.1,0.9,0.05,0.16), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==paste0('chr',CHR) &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.0, arrowlength = 0.001, bheight = 0.1, packrow=F)

labelgenome(CHR,START,END,n=5,scale="Mb")

dev.off()

######################################################################################################3
#3T3

WTcool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/TCF1_3T3_5k.cool'
KOcool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/Norm_3T3_5k.cool'
setwd("/mnt/data0/sora/TCF1_revision/Reproduce/Figure5/")

cd8a = "6:71000000-71850000"

geneinfo = read.delim('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed',header=F)
TAD = read.delim('TADs.bed')
# Bcl6
co = cd8a
name="cd8a"
CHR = unlist(strsplit(co, split=":"))[1]
START = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[1])
END = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[2])
OUT = 'TCF1_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",WTcool," --coord ",co," --norm VC_SQRT --out ", OUT," --res 5000"))
OUT = 'NORM_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",KOcool," --coord ",co," --norm VC_SQRT --out ", OUT," --res 5000"))

WT = read.delim(paste0('TCF1_0.txt'),row.names=1)
colnames(WT) = rownames(WT)
KO = read.delim(paste0('NORM_0.txt'),row.names=1)
colnames(KO) = rownames(KO)

pdf('Cd8a_3T3_R_heatmap4.pdf', width = 10, height = 8)
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.9,0.6,1))
phic_WT = plotHic(WT, CHR ,START,END, max_y = 100,
                  zrange=c(0,1.7), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.9,0.56,0.6), new=T)
plotBed(beddata = TAD,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)


par(fig=c(0.1,0.9,0.2,0.6), new=T)
phic_KO = plotHic(KO, CHR, START, END, max_y = 100,
                  zrange=c(0,1.7), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.9,0.16,0.2), new=T)
plotBed(beddata = TAD,chrom = paste0('chr',CHR),chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.1,0.9,0.05,0.16), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==paste0('chr',CHR) &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.0, arrowlength = 0.001, bheight = 0.1, packrow=F)

labelgenome(CHR,START,END,n=5,scale="Mb")

dev.off()



###### Comparison of contact frequency (in 3T3, DN3 and DP)
setwd('/mnt/data0/sora/TCF1_revision/Reproduce/Figure4/Sushi/selected/')
options(scipen=-100)

DP_WTcool = '/mnt/alvand/sora/TCF1_revision/NovaSeq/Arima-HiC_VavCre_DP_223231022/workdir/HiC_pro/hic_results/data/rawdata/DP.cool'
DP_KOcool = '/mnt/alvand/sora/TCF1_revision/NovaSeq/Arima-HiC_fl_fl_VavCre_DP_223231022/workdir/HiC_pro/hic_results/data/rawdata/TCF1_KO_DP.cool'

DN_WTcool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/DN3EV.cool'
DN_KOcool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/DN3KO.cool'

TCF_3T3_cool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/TCF1_3T3_5k.cool'
Norm_3T3_cool = '/mnt/data0/sora/TCF1_revision/Reproduce/data/hic/Norm_3T3_5k.cool'

# DN3
co1 = "6:71435000-71850000"
co2 = "6:71000000-71435000"

OUT = 'Diff_DN_WT_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",DN_WTcool," --coord1 ",co1," --coord2 ",co2," --norm VC_SQRT --out ", OUT," --res 5000"))
OUT = 'Diff_DN_KO_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",DN_KOcool," --coord1 ",co1," --coord2 ",co2," --norm VC_SQRT --out ", OUT," --res 5000"))
A = read.delim('Diff_DN_WT_0.txt',row.names = 1)
B = read.delim('Diff_DN_KO_0.txt',row.names = 1)

A = data.matrix(A)
B = data.matrix(B)
par(fig=c(0.1,0.9,0.1,0.9))
name  ='Cd8a'
pdf(file = paste0(name,'_DN3_yellow_box_difference.pdf'), width = 3,height = 5)
boxplot(list(DN3_WT = A, DN3_KO = B), las=1, ylab="Normalized contact frequency (VC_SQRT)",pch=19,cex=0.5, ylim=c(0,10))
p = wilcox.test(as.numeric(A),as.numeric(B), paired = T)
p$p.value
text(x=1.5, y=0.01, labels = paste0('P = ',signif(p$p.value, 3)))
dev.off()

# DP
co1 = "chr6:71435000-71850000"
co2 = "chr6:71000000-71435000"
OUT = 'Diff_DP_WT_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",DP_WTcool," --coord1 ",co1," --coord2 ",co2," --norm VC_SQRT --out ", OUT," --res 5000"))
OUT = 'Diff_DP_KO_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",DP_KOcool," --coord1 ",co1," --coord2 ",co2," --norm VC_SQRT --out ", OUT," --res 5000"))
A = read.delim('Diff_DP_WT_0.txt',row.names = 1)
B = read.delim('Diff_DP_KO_0.txt',row.names = 1)

A = data.matrix(A)
B = data.matrix(B)
par(fig=c(0.1,0.9,0.1,0.9))
name  ='Cd8a'
pdf(file = paste0(name,'_DP_yellow_box_difference.pdf'), width = 3,height = 5)
boxplot(list(DP_WT = A, DP_KO = B), las=1, ylab="Normalized contact frequency (VC_SQRT)",pch=19,cex=0.5, ylim=c(0,10))
p = wilcox.test(as.numeric(A),as.numeric(B), paired = T)
p$p.value
text(x=1.5, y=0.01, labels = paste0('P = ',signif(p$p.value, 3)))
dev.off()

# 3T3

co1 = "6:71435000-71850000"
co2 = "6:71000000-71435000"
OUT = 'Diff_TCF1_3T3_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",TCF_3T3_cool," --coord1 ",co1," --coord2 ",co2," --norm VC_SQRT --out ", OUT," --res 5000"))
OUT = 'Diff_Norm_3T3_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",Norm_3T3_cool," --coord1 ",co1," --coord2 ",co2," --norm VC_SQRT --out ", OUT," --res 5000"))
A = read.delim('Diff_TCF1_3T3_0.txt',row.names = 1)
B = read.delim('Diff_Norm_3T3_0.txt',row.names = 1)

A = data.matrix(A)
B = data.matrix(B)
par(fig=c(0.1,0.9,0.1,0.9))
name  ='Cd8a'
pdf(file = paste0(name,'_3T3_yellow_box_difference.pdf'), width = 3,height = 5)
boxplot(list(TCF1_3T3 = A, Norm_3T3 = B), las=1, ylab="Normalized contact frequency (VC_SQRT)",pch=19,cex=0.5, ylim=c(0,10))
p = wilcox.test(as.numeric(A),as.numeric(B), paired = T)
p$p.value
text(x=1.5, y=0.01, labels = paste0('P = ',signif(p$p.value, 3)))
dev.off()