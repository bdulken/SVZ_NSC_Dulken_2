#============================================================================================================================================#
#============================================================================================================================================#
#Group composition
#FIGURE S3A
rm(list=ls())
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])



group1_percents<-c((length(group1[grepl("qNSC",group1)])/length(group1)),(length(group1[grepl("aNSC",group1)])/length(group1)),(length(group1[grepl("NPC",group1)])/length(group1)))

group2_percents<-c((length(group2[grepl("qNSC",group2)])/length(group2)),(length(group2[grepl("aNSC",group2)])/length(group2)),(length(group2[grepl("NPC",group2)])/length(group2)))

group3_percents<-c((length(group3[grepl("qNSC",group3)])/length(group3)),(length(group3[grepl("aNSC",group3)])/length(group3)),(length(group3[grepl("NPC",group3)])/length(group3)))

group4_percents<-c((length(group4[grepl("qNSC",group4)])/length(group4)),(length(group4[grepl("aNSC",group4)])/length(group4)),(length(group4[grepl("NPC",group4)])/length(group4)))

group5_percents<-c((length(group5[grepl("qNSC",group5)])/length(group5)),(length(group5[grepl("aNSC",group5)])/length(group5)),(length(group5[grepl("NPC",group5)])/length(group5)))

groups<-rbind(group1_percents,
              group2_percents,
              group3_percents,
              group4_percents,
              group5_percents)

names<-c("qNSC-like","aNSC-early","aNSC-mid","aNSC-late","NPC-like")

groups_data<-data.frame(names,groups)

library(reshape)
groups_melt<-melt(groups_data,id.var="names")


p<-ggplot(groups_melt)

p<-p+geom_bar(aes(x=factor(groups_melt$names,levels=c("qNSC-like","aNSC-early","aNSC-mid","aNSC-late","NPC-like"),ordered=T),
                  y=as.numeric(groups_melt$value),fill=groups_melt$variable),stat="identity",width=0.6)
p<-p+labs(x="Cell grouping by monocle ordering",y="Cumulative fraction of FACS-sorted cell type")
p<-p+theme_classic()
p<-p+theme(legend.position="none")
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+scale_fill_manual(values=c("#9966FF","#ff3333","#00CCFF"))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
pdf("GroupComposition_5groups_50minModels_11052015.pdf",height=9,width=11)
print(p)
dev.off()

