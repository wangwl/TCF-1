# Golnaz Vahedi PhD
# December 13, 2021
# Bed files are generated by count_bam_files_union.sh and DESeq_pvalue.R
cd /mnt/alvand/golnaz/projects/TCF1_3D/Smc1_Nipbl_union
background=DESeq_Union_Nipbl_Smc1_No_change_0.5_0.05.bed
findMotifsGenome.pl DESeq_Union_Nipbl_Smc1_Up_significant_0.5_0.05.bed mm10 motif_lost_in_KO/ -bg /mnt/alvand/golnaz/projects/TCF1_3D/Smc1_Nipbl_union/${background} -size given
findMotifsGenome.pl DESeq_Union_Nipbl_Smc1_Down_significant_0.5_0.05.bed mm10 motif_gained_in_KO/ -bg /mnt/alvand/golnaz/projects/TCF1_3D/Smc1_Nipbl_union/${background} -size given
