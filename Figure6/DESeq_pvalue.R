####################################################################
# Date: December 13, 2021
# Golnaz Vahedi, PhD
# Revision of TCF1 3D paper for Nature Immunology
# The UNION is the list of all peaks in Nipbl and Smc1 in WT and KO DN3 T cells
####################################################################
# WARNING WARNING WARNING WARNING about DESEq
# The order in sampleSheet.txt has to be exactly the order of countdata which pretty much depends on how dir orders file names....
####################################################################
####################################################################
rm(list=ls())
####################################################################
library('DESeq2')
library('ggplot2')
library("gplots")
library("corrplot")
library('pheatmap')
####################################################################
foldChange <- 0.5
pvalue_cutoff <- 5e-2
species='mm10'
min_count = 5 #remove rows with low entries
####################################################################
ReadFiles <- function(count_file)
{
  tag_counts = read.table(paste(filedir,count_file,sep=""),
                           header=F,stringsAsFactors=F,sep='\t')
  return(tag_counts[,4])
}
####################################################################
mydir='/mnt/alvand/golnaz/projects/TCF1_3D/Smc1_Nipbl_union/'
setwd(mydir)
foldername_annotation='Union_Nipbl_Smc1'
####################################################################
filedir=paste(mydir,"counts/",sep='')
####################################################################
coord_file=dir(filedir)[1] # This gives the name of one of the count files
coordinates = read.table(paste(filedir,coord_file,sep=""),
                         header=F,stringsAsFactors=F,sep='\t')
colnames(coordinates)=c('chr','start','end')
####################################################################
# WARNING WARNING WARNING WARNING
#Important note: the order of these should be the same as the order of columns in countdata
sampleTable = read.table(paste(mydir,"/sampleSheet.txt",sep=''), 
                         header=T, sep="\t", stringsAsFactors=T)
####################################################################
# Preparing data for DESeq
myfiles = dir(filedir)
list_files_all = myfiles[grep("_All_Nipbl_Smc1.count",myfiles)]
countdata = (sapply(list_files_all,ReadFiles,simplify="array"))
TPM=countdata
rownames(countdata) = paste('peakid_',1:dim(countdata)[1],sep='')
rownames(TPM) = paste('peakid_',1:dim(countdata)[1],sep='')
####################################################################
####################### DESeq calculations #########################
####################################################################
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleTable,
                                 design = ~ Genotype)
# estimate the size factors
cds = estimateSizeFactors(ddsMat, locfunc =  median)

ddsMat$Vector <- factor(ddsMat$Genotype, levels=c("KO","EV"))
ddsMat$Vector <- relevel(ddsMat$Genotype, ref="EV")
# Major part of DESeq calculation
dds <- DESeq(ddsMat)
coordinates = coordinates[rowSums(counts(dds)) > min_count,]
dds <- dds[ rowSums(counts(dds)) > min_count, ]

DESeq_normalization_factors <- sizeFactors(dds) # this is the size factor sj calculated by DESeq
#rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
coordinates_oredered <- coordinates[order(res$pvalue),1:3]

# Take NA out
resOrdered_filterNA <- resOrdered[is.na(resOrdered$pvalue)!=T,]
coordinates_oredered_filterNA <- coordinates_oredered[is.na(resOrdered$pvalue)!=T,1:3]

resOrdered = resOrdered_filterNA
coordinates_oredered = coordinates_oredered_filterNA

####################################################################
####################### Coordinate calculations #########################
####################################################################
#Coordinates of all enahncers
coordinates_FC_pvalue_ordered = cbind(coordinates_oredered,resOrdered)

write.csv(as.data.frame(coordinates_FC_pvalue_ordered),
          file=paste("DESeq_",foldername_annotation,"_results",".csv",sep=''))
write.table(as.data.frame(coordinates_FC_pvalue_ordered[order(coordinates_FC_pvalue_ordered[,1],coordinates_FC_pvalue_ordered[,2],decreasing=F),]),
          file=paste("DESeq_",foldername_annotation,"_results",".bed",sep=''),row.names=F,quote=F,sep='\t',col.names=F)
#
Down_significant = coordinates_FC_pvalue_ordered[(resOrdered$log2FoldChange > foldChange & 
                                         resOrdered$pvalue < pvalue_cutoff),]
#
Up_significant = coordinates_FC_pvalue_ordered[(resOrdered$log2FoldChange < -foldChange & 
                                          resOrdered$pvalue < pvalue_cutoff),]

No_change = coordinates_FC_pvalue_ordered[!(resOrdered$log2FoldChange > foldChange & 
                                                  resOrdered$pvalue < pvalue_cutoff) &
                                                 !(resOrdered$log2FoldChange < -foldChange & 
                                                    resOrdered$pvalue < pvalue_cutoff),]
# For checking pvalue and FC
write.csv(as.data.frame(Down_significant),file=paste("DESeq_",foldername_annotation,"_Down_significant_",foldChange,"_",pvalue_cutoff,".csv",sep=''),row.names=F)
write.csv(as.data.frame(Up_significant),file=paste("DESeq_",foldername_annotation,"_Up_significant_",foldChange,"_",pvalue_cutoff,".csv",sep=''),row.names=F)
# For IGV

Down_significant = Down_significant[order(Down_significant[,1],Down_significant[,2],decreasing=F),]
Up_significant = Up_significant[order(Up_significant[,1],Up_significant[,2],decreasing=F),]
No_change = No_change[order(No_change[,1],No_change[,2],decreasing=F),]

write.table(as.data.frame(Down_significant[,1:3]),file=paste("DESeq_",foldername_annotation,"_Down_significant"
                                                               ,"_",foldChange,"_",pvalue_cutoff,".bed",sep=''),row.names=F,quote=F,sep='\t',col.names=F)
write.table(as.data.frame(Up_significant[,1:3]),file=paste("DESeq_",foldername_annotation,"_Up_significant"
                                                              ,"_",foldChange,"_",pvalue_cutoff,".bed",sep=''),row.names=F,quote=F,sep='\t',col.names=F)

write.table(as.data.frame(No_change[,1:3]),file=paste("DESeq_",foldername_annotation,"_No_change"
                                                            ,"_",foldChange,"_",pvalue_cutoff,".bed",sep=''),row.names=F,quote=F,sep='\t',col.names=F)

####################################################################
####################################################################
####################################################################

####################################################################
pdf("All_ChIPSeq_PCA_DESeq.pdf", useDingbats=FALSE,height=8, width=8)
p = plotPCA(vsd, intgroup=c("Genotype"))
plot(p)
dev.off()
####################################################################

####################################################################
####################################################################
#res <- results(dds)
#table(res$pvalue<pvalue_cutoff)
## Order by adjusted p-value
resOrdered <- res[order(res$pvalue),]
coordinates_oredered <- coordinates[order(res$pvalue),1:3]
## Merge with normalized count data
rownames(coordinates_oredered)=rownames(resOrdered)
resdata1 <- merge(as.data.frame(coordinates_oredered), as.data.frame(resOrdered),  by="row.names", sort=FALSE)
rownames(resdata1)=rownames(resOrdered)
resdata <- merge(resdata1,as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file=paste("DESeq_",foldername_annotation,"_diff-results_",format(Sys.time(), "%b_%d_%Y_%H_%M%p"),".csv",sep = ''))

####################################################################
####################################################################
####################################################################
####################################################################
## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=foldChange, sigthresh=pvalue_cutoff, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, pvalue<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  # with(subset(res, pvalue<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  #with(subset(res, pvalue<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  with(subset(res, pvalue<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  
  if (labelsig) {
    require(calibrate)
    with(subset(res, pvalue<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend( x="topleft", y="topleft", legend=c(paste("Number of Down regions in WT is ",dim(Down_significant)[1],sep=""), 
                                                   paste("Number of Up regions in WT is ",dim(Up_significant)[1],sep="")), pch=20, col=c("green","green"))
}
pdf(paste("DESeq_",foldername_annotation,"-volcanoplot.pdf",sep=''))
volcanoplot(resdata, lfcthresh=foldChange, sigthresh=pvalue_cutoff, textcx=.8, ylim=c(0, 
                                                                                      30),
            xlim=c(-5, 5))
dev.off()
####################################################################
####################################################################
# Save the data
list_files_all = myfiles[grep("*.RData",myfiles)]
file.remove(list_files_all) #remove previous clusters specially if K changes.
save.image(file=paste(mydir,'RData_',foldername_annotation,"_",format(Sys.time(), "%b_%d_%Y_%H_%M%p"),'.RData',sep=''))



