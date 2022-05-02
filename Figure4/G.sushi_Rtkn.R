# 2021 12 22
# Sushi plot for cd8a region.

library(Sushi)

setwd("/mnt/data0/sora/TCF1_revision/Reproduce/Figure4/Sushi/selected/")

DP_WT_cool = '/mnt/alvand/sora/TCF1_revision/NovaSeq/Arima-HiC_VavCre_DP_223231022/workdir/HiC_pro/hic_results/data/rawdata/DP.cool'
DP_KO_cool = '/mnt/alvand/sora/TCF1_revision/NovaSeq/Arima-HiC_fl_fl_VavCre_DP_223231022/workdir/HiC_pro/hic_results/data/rawdata/TCF1_KO_DP.cool'

Fib_ind_cool = '/mnt/data0/wenliang/project/04.TCF1/16.GEO_submission/HiChIP/processed_data/TCF1_3T3.mcool::resolutions/5000'
Fib_wt_cool = '/mnt/data0/wenliang/project/04.TCF1/16.GEO_submission/HiChIP/processed_data/Norm_3T3.mcool::resolutions/5000'

DN_WT_cool = '/mnt/data0/wenliang/project/04.TCF1/09.DN3_KO/02.hic_cool/DN3EV.mcool::resolutions/5000'
DN_KO_cool = '/mnt/data0/wenliang/project/04.TCF1/09.DN3_KO/02.hic_cool/DN3KO.mcool::resolutions/5000'
geneinfo = read.delim('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed',header=F)

# 
co = "chr6:82950000-83350000"
name="Wdr54"
CHR = unlist(strsplit(co, split=":"))[1]
START = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[1])
END = as.numeric(unlist(strsplit(unlist(strsplit(co, split=":"))[2], split="-"))[2])
OUT = 'DP_WT_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",DP_WT_cool," --coord ",co," --norm weight --out ", OUT," --res 5000"))
OUT = 'DP_KO_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",DP_KO_cool," --coord ",co," --norm weight --out ", OUT," --res 5000"))

OUT = 'DN3_WT_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",DN_WT_cool," --coord ",co," --norm weight --out ", OUT," --res 5000"))
OUT = 'DN3_KO_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",DN_KO_cool," --coord ",co," --norm weight --out ", OUT," --res 5000"))

OUT = 'TCF1_3T3_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",Fib_ind_cool," --coord ",co," --norm weight --out ", OUT," --res 5000"))
OUT = 'Norm_3T3_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",Fib_wt_cool," --coord ",co," --norm weight --out ", OUT," --res 5000"))

DP_WT = read.delim(paste0('DP_WT_0.txt'),row.names=1)
colnames(DP_WT) = rownames(DP_WT)
DP_KO = read.delim(paste0('DP_KO_0.txt'),row.names=1)
colnames(DP_KO) = rownames(DP_KO)

DN_WT = read.delim(paste0('DN3_WT_0.txt'),row.names=1)
colnames(DN_WT) = rownames(DN_WT)
DN_KO = read.delim(paste0('DN3_KO_0.txt'),row.names=1)
colnames(DN_KO) = rownames(DN_KO)

TCF1_3T3 = read.delim(paste0('TCF1_3T3_0.txt'),row.names=1)
colnames(TCF1_3T3) = rownames(TCF1_3T3)
NORM_3T3 = read.delim(paste0('Norm_3T3_0.txt'),row.names=1)
colnames(NORM_3T3) = rownames(NORM_3T3)

## BigWIg to Bedgraph
system(
paste0("
BB=/mnt/data0/wenliang/software/UCSC/linux.x86_64/bigWigToBedGraph;
name=",name,";
gene=/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed;
folder=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ChIPseq/Thymus/03.bigwig/;
TCF1_DP=ChIP_seq_DP_TCF1_195748578_S7;
folder_dn=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ChIPseq/DN3_Scid.adh/03.bigwig/;
TCF1_DN=ChIPseq_TCF1_DN3-EV_Rep1_204407220_S1;
folder_3t3=/mnt/alvand/VahediLab/GenomicsData/Mouse/NIH3T3/ChIPseq/03.bigwig/;
TCF1_3T3=ChIP_seq_3T3_dox_72hrs_TCF1_rep1_195748578_S5;

CTCF_DP=ChIP-seq_CTCF_DP_BL6_Rep1_87983902_S3;
CTCF_DN=ChIPseq_CTCF_DN3_EV_10M_Rep1_203364293_S4;
CTCF_3T3=ChIP_seq_CTCF_3T3_pIND20_TCF_Rep1T_191047891_S3;


chr=",gsub("chr","",CHR),";
start=",START,";
end=",END,";
cd /mnt/data0/sora/TCF1_revision/Reproduce/Figure5/
$BB ${folder}${TCF1_DP}.bw ${TCF1_DP}.${name}.bedgraph -chrom=$chr -start=$start -end=$end; 
$BB ${folder_dn}${TCF1_DN}.bw ${TCF1_DN}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder_3t3}${TCF1_3T3}.bw ${TCF1_3T3}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder}${CTCF_DP}.bw ${CTCF_DP}.${name}.bedgraph -chrom=$chr -start=$start -end=$end; 
$BB ${folder_dn}${CTCF_DN}.bw ${CTCF_DN}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder_3t3}${CTCF_3T3}.bw ${CTCF_3T3}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;") 
)

TCF1_ChIP_DP= read.delim(paste0('../../../Figure5/ChIP_seq_DP_TCF1_195748578_S7.',name,'.bedgraph'), header=F)
TCF1_ChIP_DN= read.delim(paste0('../../../Figure5/ChIPseq_TCF1_DN3-EV_Rep1_204407220_S1.',name,'.bedgraph'), header=F)
TCF1_ChIP_3T3= read.delim(paste0('../../../Figure5/ChIP_seq_3T3_dox_72hrs_TCF1_rep1_195748578_S5.',name,'.bedgraph'), header=F)
CTCF_ChIP_DP= read.delim(paste0('../../../Figure5/ChIP-seq_CTCF_DP_BL6_Rep1_87983902_S3.',name,'.bedgraph'), header=F)
CTCF_ChIP_DN= read.delim(paste0('../../../Figure5/ChIPseq_CTCF_DN3_EV_10M_Rep1_203364293_S4.',name,'.bedgraph'), header=F)
CTCF_ChIP_3T3= read.delim(paste0('../../../Figure5/ChIP_seq_CTCF_3T3_pIND20_TCF_Rep1T_191047891_S3.',name,'.bedgraph'), header=F)



TCF1_ChIP_DP$V1 = paste0('chr',TCF1_ChIP_DP$V1)
TCF1_ChIP_DN$V1 = paste0('chr',TCF1_ChIP_DN$V1)
TCF1_ChIP_3T3$V1 = paste0('chr',TCF1_ChIP_3T3$V1)
CTCF_ChIP_DP$V1 = paste0('chr',CTCF_ChIP_DP$V1)
CTCF_ChIP_DN$V1 = paste0('chr',CTCF_ChIP_DN$V1)
CTCF_ChIP_3T3$V1 = paste0('chr',CTCF_ChIP_3T3$V1)

## TAD INFO 
TAD.TCF1 = read.delim('/mnt/data0/wenliang/project/04.TCF1/03.HiChIP/03.TADs/05.cooltools/02.TADs/TCF1_3T3.TADs', header=F)
TAD.Norm = read.delim('/mnt/data0/wenliang/project/04.TCF1/03.HiChIP/03.TADs/05.cooltools/02.TADs/Norm_3T3.TADs', header=F)
TAD.DN3EV = read.delim('/mnt/data0/wenliang/project/04.TCF1/09.DN3_KO/05.TAD/01.DN3EV_DN3KO/DN3EV.TADs', header=F)
TAD.DN3KO = read.delim('/mnt/data0/wenliang/project/04.TCF1/09.DN3_KO/05.TAD/01.DN3EV_DN3KO/DN3KO.TADs', header=F)
TAD.DPWT = read.delim('/mnt/data0/sora/TCF1_revision/Reproduce/Figure4/DP_WT.10k.TAD.bed', header=F)
TAD.DPKO = read.delim('/mnt/data0/sora/TCF1_revision/Reproduce/Figure4/DP_TCF1_KO.10k.TAD.bed', header=F)
TAD.DPWT$V1 = paste0('chr',TAD.DPWT$V1)
TAD.DPKO$V1 = paste0('chr',TAD.DPKO$V1)

zmax_3t3=0.003
zmax_dn=0.0023
zmax_dp=0.0008
max_y = 40
pdf(paste0(name,'_R_heatmap.pdf'), width = 12, height = 8)
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.4,0.7,1))
phic_WT = plotHic(TCF1_3T3, CHR ,START,END, max_y = max_y,
                   zrange=c(0,zmax_3t3), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.4,0.7,0.75), new=T)
plotBed(beddata = TAD.TCF1,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)



par(fig=c(0.1,0.4,0.4,0.7), new=T)
phic_KO = plotHic(NORM_3T3, CHR, START, END, max_y = max_y,
                   zrange=c(0,zmax_3t3), palette=colorRampPalette(c("white",'red'))
)

par(fig=c(0.1,0.4,0.4,0.45), new=T)
plotBed(beddata = TAD.Norm,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.1,0.4,0.3,0.4), new=T)
plotBedgraph(TCF1_ChIP_3T3, chrom = CHR, chromstart = START, chromend = END, addscale = T, range = c(0,2), color = "#3A53A4")

par(fig=c(0.1,0.4,0.2,0.3), new=T)
plotBedgraph(CTCF_ChIP_3T3, chrom = CHR, chromstart = START, chromend = END, addscale = T, range = c(0,8.5), color = "#ED2224")


par(fig=c(0.1,0.4,0.05,0.2), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==CHR &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.0, arrowlength = 0.001, bheight = 0.1, packrow=F)




labelgenome(CHR,START,END,n=5,scale="Mb")

par(fig=c(0.4,0.7,0.7,1), new=T)
phic_WT = plotHic(DN_WT, CHR ,START,END, max_y = max_y,
                  zrange=c(0,zmax_dn), palette=colorRampPalette(c("white",'red'))
)

par(fig=c(0.4,0.7,0.7,0.75), new=T)
plotBed(beddata = TAD.DN3EV,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.4,0.7,0.4,0.7), new=T)
phic_KO = plotHic(DN_KO, CHR, START, END, max_y = max_y,
                  zrange=c(0,zmax_dn), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.4,0.7,0.4,0.45), new=T)
plotBed(beddata = TAD.DN3KO,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.4,0.7,0.3,0.4), new=T)
plotBedgraph(TCF1_ChIP_DN, chrom = CHR, chromstart = START, chromend = END, addscale = T, range = c(0,1.5))

par(fig=c(0.4,0.7,0.2,0.3), new=T)
plotBedgraph(CTCF_ChIP_DN, chrom = CHR, chromstart = START, chromend = END, addscale = T, range = c(0,1.5))


par(fig=c(0.4,0.7,0.05,0.2), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==CHR &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.0, arrowlength = 0.001, bheight = 0.1, packrow=F)
labelgenome(CHR,START,END,n=5,scale="Mb")

par(fig=c(0.7,0.99,0.7,1), new=T)
phic_WT = plotHic(DP_WT, CHR ,START,END, max_y = max_y,
                  zrange=c(0,zmax_dp), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.7,0.99,0.7,0.75), new=T)
plotBed(beddata = TAD.DPWT,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.7,0.99,0.4,0.7), new=T)
phic_KO = plotHic(DP_KO, CHR, START, END, max_y = max_y,
                  zrange=c(0,zmax_dp), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.7,0.99,0.4,0.45), new=T)
plotBed(beddata = TAD.DPKO,chrom = CHR,chromstart = START,
        chromend =END,
        colorbycol = SushiColors(2),row = "auto",wiggle=0.001)

par(fig=c(0.7,0.99,0.3,0.4), new=T)
plotBedgraph(TCF1_ChIP_DP, chrom = CHR, chromstart = START, chromend = END, addscale = T, range = c(0,0.8))

par(fig=c(0.7,0.99,0.2,0.3), new=T)
plotBedgraph(CTCF_ChIP_DP, chrom = CHR, chromstart = START, chromend = END, addscale = T, range = c(0,1.5))

par(fig=c(0.7,0.99,0.05,0.2), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==CHR &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.0, arrowlength = 0.001, bheight = 0.1, packrow=F)

# par(fig=c(0.4,0.7,0.4,0.5), new=T)
# plotBedgraph(nod_ctcf ,CHR,START,END,color='#22854a',range=c(0,3),addscale = T)
# 
# par(fig=c(0.4,0.7,0.3,0.4), new=T)
# plotBedgraph(nod_smc1, CHR,START,END,color='#ff2185',range=c(0,4),addscale = T)
# 
# par(fig=c(0.4,0.7,0.2,0.3), new=T)
# plotBedgraph(nod_k27, CHR,START,END,color='#0600b5',range=c(0,3),addscale = T)
# 
# par(fig=c(0.4,0.7,0.1,0.2), new=T)
# plotBedgraph(nod_atac, CHR,START,END,color='black',range=c(0,2),addscale = T)

labelgenome(CHR,START,END,n=5,scale="Mb")

dev.off()


###### Comparison of contact frequency (in 3T3)
setwd('/mnt/data0/sora/TCF1_revision/Reproduce/Figure4/Sushi/selected/')
options(scipen=-100)

Fib_ind_cool = '/mnt/data0/wenliang/project/04.TCF1/16.GEO_submission/HiChIP/processed_data/TCF1_3T3.mcool::resolutions/5000'
Fib_wt_cool = '/mnt/data0/wenliang/project/04.TCF1/16.GEO_submission/HiChIP/processed_data/Norm_3T3.mcool::resolutions/5000'

co1 = "chr6:82950000-83130000"
co2 = "chr6:83200000-83350000"

name="Wdr54"
CHR1 = unlist(strsplit(co1, split=":"))[1]
START1 = as.numeric(unlist(strsplit(unlist(strsplit(co1, split=":"))[2], split="-"))[1])
END1 = as.numeric(unlist(strsplit(unlist(strsplit(co1, split=":"))[2], split="-"))[2])
CHR2 = unlist(strsplit(co2, split=":"))[1]
START2 = as.numeric(unlist(strsplit(unlist(strsplit(co2, split=":"))[2], split="-"))[1])
END2 = as.numeric(unlist(strsplit(unlist(strsplit(co2, split=":"))[2], split="-"))[2])


OUT = 'Diff_TCF1_3T3_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",Fib_ind_cool," --coord1 ",co1," --coord2 ",co2," --norm weight --out ", OUT," --res 5000"))
OUT = 'Diff_Norm_3T3_'
system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq_noSymmetric.py --cool ",Fib_wt_cool," --coord1 ",co1," --coord2 ",co2," --norm weight --out ", OUT," --res 5000"))

A = read.delim('Diff_TCF1_3T3_0.txt',row.names = 1)
B = read.delim('Diff_Norm_3T3_0.txt',row.names = 1)

A = data.matrix(A)
B = data.matrix(B)
par(fig=c(0.1,0.9,0.1,0.9))
pdf(file = 'Wdr54_yellow_box_difference.pdf', width = 3,height = 5)
boxplot(list(TCF1_3T3 = A, Norm_3T3 = B), las=1, ylab="Normalized contact frequency (balancing)",pch=19,cex=0.5)
p = wilcox.test(A,B, paired=T)
p$p.value
text(x=1, y=0.008, labels = paste0('P = ',signif(p$p.value, 3)))
dev.off()
