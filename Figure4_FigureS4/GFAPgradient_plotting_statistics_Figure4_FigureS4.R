#Ben Dulken
#4/26/16
#This script performs plotting and statistical testing for the experiments relating to validating
#the finding that the aNSC gate contains multiple populations.

#Multiple groups of aNSCs (GFAP-high, GFAP-mid, GFAP-low) were sorted along with NPCs
#qPCRs were performed on genes of interest.

#============================================================================================================================================#
#Plots box/dot plots for specific genes
#Figure 4B, S4B
rm(list=ls())

library(ggplot2)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
expression<-read.table("allSorts_average.txt")

expression_compressed<-aggregate(t(as.matrix(expression)),list(c(rep(c("aNSC_early_1","aNSC_mid_1","aNSC-late_1","NPC_1"),4))),mean)

factor<-factor(c(rep(c("GFAP-high\naNSC","GFAP-mid\naNSC","GFAP-low\naNSC","NPC"),4)))
color<-factor(c(rep(c("#FFFFFF","#FFFFFF","#ff3333","#00CCFF"),4)))
exp_color<-factor(c(rep("#000000",16)))

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure4_Output")
for(i in 1:length(expression[,1])){
  expression_expanded<-c(rep(expression_compressed[,i+1],4))
  data<-data.frame(t(expression[i,]),expression_expanded,factor=factor,color=color)
  data$factor<-as.character(data$factor)
  data$factor <- factor(data$factor, levels=unique(data$factor), ordered = T) 
  data$color<-as.character(data$color)
  data$color <- factor(data$color, levels=unique(data$color), ordered = T) 
  p<-ggplot(data)
  p<-p+geom_boxplot(aes(x=factor,y=data[,1], fill=color))
  p<-p+geom_point(aes(x=factor,y=data[,1]),color=exp_color,size=5,alpha=0.7)
  p<-p+geom_line(aes(x=factor,y=data[,2],group=exp_color,color=exp_color),size=1.5)
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p+labs(y="Gene expression\n(normalized to B-Actin)")
  p<-p+theme(legend.position = "none") 
  p<-p+theme(axis.text.x=element_text(size=18))
  p<-p+theme(axis.title.x=element_blank())
  p<-p+theme(axis.text.y=element_text(size=22))
  p<-p+theme(axis.title.y=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=22))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+scale_fill_manual(values=c("#FFFFFF","#ff3333","#00CCFF"))
  p<-p+scale_color_manual(values=c("#000000"))
  pdf(paste(colnames(data)[1],"_allsorts_average.pdf",sep=""),width=8,height=4)
  print(p)
  dev.off()
}

#============================================================================================================================================#
#Statistical testing
ttest_results<-c()
for(i in 1:length(expression[,1])){
  curr<-expression[i,]
  ttest1<-wilcox.test(as.numeric(curr[c(1,5,9,13)]),as.numeric(curr[c(2,6,10,14)]),alternative="less")
  ttest2<-wilcox.test(as.numeric(curr[c(1,5,9,13,2,6,10,14)]),as.numeric(curr[c(3,7,11,15)]),alternative="less")
  ttest3<-wilcox.test(as.numeric(curr[c(3,7,11,15)]),as.numeric(curr[c(4,8,12,16)]),alternative="less")
  ttest_results<-rbind(ttest_results,c(ttest1$p.value,ttest2$p.value,ttest3$p.value))
}
ttest_results

ttest_results<-c()
for(i in 1:length(expression[,1])){
  curr<-expression[i,]
  ttest1<-wilcox.test(as.numeric(curr[c(1,5,9,13)]),as.numeric(curr[c(2,6,10,14)]),alternative="greater")
  ttest2<-wilcox.test(as.numeric(curr[c(1,5,9,13,2,6,10,14)]),as.numeric(curr[c(3,7,11,15)]),alternative="greater")
  ttest3<-wilcox.test(as.numeric(curr[c(3,7,11,15)]),as.numeric(curr[c(4,8,12,16)]),alternative="greater")
  ttest_results<-rbind(ttest_results,c(ttest1$p.value,ttest2$p.value,ttest3$p.value))
}
ttest_results



#============================================================================================================================================#
#Correlation plot for only aNSC-mid and aNSC-late cells
#Figure 4C

rm(list=ls())

library(ggplot2)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
expression<-read.table("allSorts_average_mod_noNPC_noEarly.txt")


#Correlation

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(expression),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for qPCR
library(gplots)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure4_Output")
pdf("Correlation_qPCR_GFAPgradient_allsorts_average_noNPC_noEarly.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()



#============================================================================================================================================#
#Correlation plot for only all replicates (aNSC-early, aNSC-mid, aNSC-late, NPC)
#Figure S4C

rm(list=ls())

library(ggplot2)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
expression<-read.table("allSorts_average_mod.txt")

expression_compressed<-aggregate(t(as.matrix(expression)),list(c(rep(c("aNSC_early_1","aNSC_mid_1","aNSC-late_1","NPC_1"),4))),mean)


#Correlation

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(expression),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for qPCR
library(gplots)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure4_Output")
pdf("Correlation_qPCR_GFAPgradient_allsorts_average_mod.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none")

dev.off()


