#Ki67 levels
rm(list=ls())
library(ggplot2)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Ki67_analysis/")
ki67<-read.table("Ki67Levels.txt",header=T)

Ast_mean<-mean(as.numeric(ki67[2,]))
qNSC_mean<-mean(as.numeric(ki67[4,]))
aNSC_mean<-mean(as.numeric(ki67[6,]))
NPC_mean<-mean(as.numeric(ki67[8,]))

Ast_sd<-sd(as.numeric(ki67[2,]))
qNSC_sd<-sd(as.numeric(ki67[4,]))
aNSC_sd<-sd(as.numeric(ki67[6,]))
NPC_sd<-sd(as.numeric(ki67[8,]))

ki67mean<-c(Ast_mean,
        qNSC_mean,
        aNSC_mean,
        NPC_mean)

sd<-c(Ast_sd,
        qNSC_sd,
        aNSC_sd,
        NPC_sd)

upper<-ki67mean+sd
lower<-ki67mean-sd
top<-100-ki67mean

wilcox.test(as.numeric(ki67[6,]),as.numeric(ki67[8,]))

data<-data.frame(names=c("Astrocytes","qNSCs","aNSCs","NPCs"),ki67mean,top)
groups_melt<-melt(data)

groups_melt<-data.frame(groups_melt,upper=c(upper,upper),lower=c(lower,lower))
p<-ggplot(groups_melt)

p<-p+geom_bar(aes(x=factor(groups_melt$names,levels=c("Astrocytes","qNSCs","aNSCs","NPCs"),ordered=T),
                  y=as.numeric(groups_melt$value),fill=groups_melt$variable),stat="identity",width=0.6)
# p<p+geom_errorbar(aes(x=factor(groups_melt$names,levels=c("Astrocytes","qNSCs","aNSCs","NPCs"),ordered=T),
#                       ymax=as.numeric(groups_melt$upper),ymin=as.numeric(groups_melt$lower)),stat="identity",width=0.6)
p<-p+labs(x="Cell type",y="Percent")
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
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
pdf("Ki67Values.pdf")
print(p)
dev.off()




setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Ki67_analysis/")
ki67_neg<-c(as.numeric(ki67[6,]),as.numeric(ki67[8,]))
ki67_neg_fac<-c(rep("aNSC",4),rep("NPC",4))

d<-data.frame(neg=ki67_neg,fac=ki67_neg_fac)

p<-ggplot(d)
p<-p+geom_boxplot(aes(x=d$fac,y=d$neg),fill=c("#FF3300","#00CCFF"))
# p<p+geom_errorbar(aes(x=factor(groups_melt$names,levels=c("Astrocytes","qNSCs","aNSCs","NPCs"),ordered=T),
#                       ymax=as.numeric(groups_melt$upper),ymin=as.numeric(groups_melt$lower)),stat="identity",width=0.6)
p<-p+labs(x="Cell type",y="Percent Ki67 negative")
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+theme(legend.position="none")
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+scale_fill_manual(values=c("#FF3300","#00CCFF"))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
pdf("Ki67Values_aNSC_NPC_only.pdf")
print(p)
dev.off()

















#Ki67 levels
rm(list=ls())
library(ggplot2)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Ki67_analysis/")
ki67<-read.table("Ki67Levels_ANSC_NPC.txt",header=T)


aNSC_mean<-mean(as.numeric(ki67[2,]))
NPC_mean<-mean(as.numeric(ki67[4,]))

aNSC_sd<-sd(as.numeric(ki67[2,]))
NPC_sd<-sd(as.numeric(ki67[4,]))

ki67mean<-c(aNSC_mean,
            NPC_mean)

sd<-c(ANSC_sd,
      NPC_sd)

upper<-ki67mean+sd
lower<-ki67mean-sd
top<-100-ki67mean

wilcox.test(as.numeric(ki67[2,]),as.numeric(ki67[4,]))



setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Ki67_analysis/")
ki67_neg<-c(as.numeric(ki67[2,]),as.numeric(ki67[4,]))
ki67_neg_fac<-c(rep("aNSC",6),rep("NPC",6))

d<-data.frame(neg=ki67_neg,fac=ki67_neg_fac)

p<-ggplot(d)
p<-p+geom_boxplot(aes(x=d$fac,y=d$neg),fill=c("#FF3300","#00CCFF"))
# p<p+geom_errorbar(aes(x=factor(groups_melt$names,levels=c("Astrocytes","qNSCs","aNSCs","NPCs"),ordered=T),
#                       ymax=as.numeric(groups_melt$upper),ymin=as.numeric(groups_melt$lower)),stat="identity",width=0.6)
p<-p+labs(x="Cell type",y="Percent Ki67 negative")
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+theme(legend.position="none")
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+scale_fill_manual(values=c("#FF3300","#00CCFF"))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
pdf("Ki67Values_aNSC_NPC_only_extrareps.pdf")
print(p)
dev.off()