# Minimum distance
library(ggplot2)
library(ggsignif)

get_dist = function(df)
{
  res = tapply(X = df$dist, INDEX = paste0(df$Label,"_",df$nuc.idx), FUN = function(x){a=sort(x,decreasing=F);return(a[1])})
  res = unlist(res)
  res = res[which(!is.infinite(res))]
  return(res)
}

WT1 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_WT1_colocalization/t2.csv')
WT2 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_WT2_colocalization/t2.csv')
KO1 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_KO1_colocalization/t2.csv')
KO2 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_KO2_colocalization/t2.csv')
KO3 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_KO3_colocalization/t2.csv')

dist_WT1 = get_dist(WT1)
dist_WT2 = get_dist(WT2)
dist_KO1 = get_dist(KO1)
dist_KO2 = get_dist(KO2)
dist_KO3 = get_dist(KO3)

dist_df = data.frame(Distance = c(dist_WT1, dist_WT2, dist_KO1, dist_KO2, dist_KO3), 
                     Label = c(rep('WT1', length(dist_WT1)),
                               rep("WT2", length(dist_WT2)),
                               rep('KO1', length(dist_KO1)),
                               rep('KO2', length(dist_KO2)),
                               rep('KO3', length(dist_KO3))),
                     type = c(rep("WT", length(dist_WT1)+length(dist_WT2)),
                              rep("KO", length(dist_KO1)+length(dist_KO2)+length(dist_KO3))) 
)

pdf('/mnt/data0/sora/TCF1_revision/Figure/boxplot_DP_TAD_distance.pdf', width = 5, height =5)
ggplot(dist_df, aes(x=Label, y=Distance, fill = type))+geom_boxplot()+theme(axis.line = element_line(colour = "black"),
                                                                           panel.grid.major = element_blank(),
                                                                           panel.grid.minor = element_blank(),
                                                                           panel.border = element_blank(),
                                                                           panel.background = element_blank())#+coord_cartesian(ylim = c(0, 3.5))
dev.off()
pdf('/mnt/data0/sora/TCF1_revision/Figure/boxplot_DP_TAD_distance_pooled.pdf', width = 5, height = 5)
ggplot(dist_df, aes(x=type, y=Distance, fill = type))+geom_boxplot() +   theme(axis.line = element_line(colour = "black"),
                                                                               panel.grid.major = element_blank(),
                                                                               panel.grid.minor = element_blank(),
                                                                               panel.border = element_blank(),
                                                                               panel.background = element_blank()) + 
  geom_signif(
  comparisons = list(c("KO", "WT")),
  map_signif_level = FALSE
)
dev.off()

nCells = data.frame(table(dist_df$Label),yloc=1)

ggplot(dist_df, aes(x=Label, y=Distance, col=type))+geom_boxplot()+ 
  geom_text(data=nCells, aes(x=Var1, y=yloc, label=Freq), col='black', size=10)

pdf('/mnt/data0/sora/TCF1_revision/Figure/cdf_DP_TAD_distance.pdf', width = 5, height = 5)
library(ggplot2)
ggplot(dist_df, aes(Distance, color=type)) +stat_ecdf(geom='line',size=3)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic()+
  xlab('TAD2-TAD3 distance (um)')+
  ylab("Percent") + theme(axis.text = element_text(size=15), 
                          axis.title=element_text(size=17,face="bold"), 
                          legend.text = element_text(size=15),
                          legend.position=c(.85,.65),
                          legend.title = element_blank())+
  annotate("text", x=3.5, y=0.5, label= "K-S test P = 2.6e-14",cex=6) 
dev.off()


########################################################################################################################33
get_overlaps = function(df)
{
  res = tapply(X = df$pcOverlap1, INDEX = paste0(df$Label,"_",df$nuc.idx), FUN = function(x){a=sort(x,decreasing=T);return(a[1])})
  res = unlist(res)
  res = res[which(!is.infinite(res))]
  return(res)
}

WT1 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_WT1_colocalization/t2_t3.csv')
WT2 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_WT2_colocalization/t2_t3.csv')
KO1 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_KO1_colocalization/t2_t3.csv')
KO2 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_KO2_colocalization/t2_t3.csv')
KO3 = read.delim('/mnt/data0/sora/TCF1_revision/11.19.2021 - Slide B1/12.18.2021-Slide B1_KO3_colocalization/t2_t3.csv')

ov_WT1 = get_overlaps(WT1)
ov_WT2 = get_overlaps(WT2)
ov_KO1 = get_overlaps(KO1)
ov_KO2 = get_overlaps(KO2)
ov_KO3 = get_overlaps(KO3)

ov_df = data.frame(Overlap = c(ov_WT1, ov_WT2, ov_KO1, ov_KO2, ov_KO3), 
                     Label = c(rep('WT1', length(ov_WT1)),
                               rep("WT2", length(ov_WT2)),
                               rep('KO1', length(ov_KO1)),
                               rep('KO2', length(ov_KO2)),
                               rep('KO3', length(ov_KO3))),
                     type = c(rep("WT", length(ov_WT1)+length(ov_WT2)),
                              rep("KO", length(ov_KO1)+length(ov_KO2)+length(ov_KO3)))
)


pdf('/mnt/data0/sora/TCF1_revision/Figure/boxplot_DP_TAD_overlap.pdf', width = 5, height =5)
ggplot(ov_df, aes(x=Label, y=Overlap, fill = type))+geom_boxplot()+
  ylab("Overlap volume fraction (normalized with TAD2)")+theme(axis.line = element_line(colour = "black"),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               panel.border = element_blank(),
                                                               panel.background = element_blank())#+coord_cartesian(ylim = c(0, 3.5))
dev.off()

pdf('/mnt/data0/sora/TCF1_revision/Figure/boxplot_DP_TAD_overlap_pooled.pdf', width = 5, height =5)
ggplot(ov_df, aes(x=type, y=Overlap, fill = type))+geom_boxplot()+
  ylab("Overlap volume fraction (normalized with TAD2)") + theme(axis.line = element_line(colour = "black"),
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank(),
                                                                  panel.border = element_blank(),
                                                                  panel.background = element_blank())+
  geom_signif(comparisons = list(c("KO","WT")))
#+coord_cartesian(ylim = c(0, 3.5))
dev.off()

pdf('/mnt/data0/sora/TCF1_revision/Figure/cdf_DP_TAD_overlap.pdf', width = 5, height = 5)
ggplot(ov_df, aes(Overlap, color=type)) +stat_ecdf(geom='line',size=3)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic()+
  xlab('Overlap volume fraction (normalized with upstream TAD)')+
  ylab("Percent") + theme(axis.text = element_text(size=15), 
                          axis.title=element_text(size=17,face="bold"), 
                          legend.text = element_text(size=15),
                          legend.position=c(.85,.25),
                          legend.title = element_blank()) +
  annotate("text", x=0.5, y=0.5, label= "K-S test P = 0.00057",cex=6)
dev.off()
