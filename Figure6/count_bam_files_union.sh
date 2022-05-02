#!/bin/bash
#-----------------------------------------------------------------------------
# December 13, 2021
# Golnaz Vahedi, PhD
# Counting and DESeq analysis for TCF-1 3D revision for Nature Immunology
#-----------------------------------------------------------------------------
genome=mm10
macs_genome=mm
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
mainDir=/mnt/alvand/golnaz/projects/TCF1_3D/Smc1_Nipbl_union
Nipbl_bamDir=/mnt/data0/aditi/TCF1_3Dgenome/TCF1-KO_DN3/ChIP_CUTnRUN/ChIPseq_NIPBL/bam
Smc1_bamDir=/mnt/data0/aditi/TCF1_3Dgenome/TCF1-KO_DN3/ChIP_CUTnRUN/ChIPseq_SMC1/bam
K27ac_bamDir=/mnt/data0/aditi/TCF1_3Dgenome/TCF1-KO_DN3/ChIP_CUTnRUN/ChIPseq_H3K27ac/bam
peakDir=/mnt/data0/aditi/TCF1_3Dgenome/TCF1-KO_DN3/ChIP_CUTnRUN/
#
countDir=${mainDir}/counts
if [ ! -d ${countDir} ]; then
 mkdir ${countDir}
fi
#-----------------------------------------------------------------------------
#echo "Create Cataloge Map combining all peaks in NIPBL and SMC1 in WT and KO DN3"
cat `ls ${peakDir}/ChIPseq_NIPBL/macs2/*NIPBL_DN3*_macs2_peaks_unblacklisted.bed` > union_all_ChIP-seq_peaks_NIPBL.bed
cat `ls ${peakDir}/ChIPseq_SMC1/macs2/*SMC1_DN3*_macs2_peaks_unblacklisted.bed` > union_all_ChIP-seq_peaks_SMC1.bed
cat union_all_ChIP-seq_peaks_NIPBL.bed union_all_ChIP-seq_peaks_SMC1.bed > union_all_ChIP-seq_peaks_cat.bed
sort -k1,1 -k2,2n union_all_ChIP-seq_peaks_cat.bed > union_all_ChIP-seq_peaks_cat_sorted.bed
bedtools merge -i union_all_ChIP-seq_peaks_cat_sorted.bed > union_SMC1_NIPBL_DN3_all.bed
rm union_all_ChIP-seq_peaks_cat.bed
rm union_all_ChIP-seq_peaks_cat_sorted.bed
rm union_all_ChIP-seq_peaks_NIPBL.bed
rm union_all_ChIP-seq_peaks_SMC1.bed
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
REFFILE=${mainDir}/union_SMC1_NIPBL_DN3_all.bed
# This reference peak list has all SMC1 and NIPBL peaks in WT and Tcf-1 KO DN3 T cells
cd $Nipbl_bamDir
echo "Counting Nipbl Library"
for i in `ls *NIPBL*.bam`;
do
 filename="${i%.*}"
 echo $filename
 bedtools coverage -counts -a $REFFILE -b ${filename}.bam  > ${countDir}/${filename}_All_Nipbl_Smc1.count
done
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
cd $Smc1_bamDir
echo "Counting Smc1 Library"
for i in `ls *SMC1*.bam`;
do
 filename="${i%.*}"
 echo $filename
 bedtools coverage -counts -a $REFFILE -b ${filename}.bam  > ${countDir}/${filename}_All_Nipbl_Smc1.count
done
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
cd $K27ac_bamDir
echo "Counting K27ac Library"
for i in `ls *K27ac*.bam`;
do
 filename="${i%.*}"
 echo $filename
 bedtools coverage -counts -a $REFFILE -b ${filename}.bam  > ${countDir}/${filename}_All_Nipbl_Smc1.count
done
