rm(list=ls())
library(ggplot2)

#This code generates the plot for the enrichments of Go Terms in Figure S1G

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
go_terms<-read.table("Oligo_go_terms.txt",header=T,sep="\t")

sig_terms<-as.vector(go_terms[5:1,1])
sig_vals<-as.vector(go_terms[5:1,4])

dat<-data.frame(terms=sig_terms,vals=-log10(sig_vals))

p<-ggplot(dat)
p<-p+geom_bar(aes(x=dat$terms,y=dat$vals),stat="identity",fill="darkgreen",color="darkgreen")
dat$terms<-as.character(dat$terms)
dat$terms <- factor(dat$terms, levels=unique(dat$terms), ordered = T)
p<-p+coord_flip()
p<-p+labs(y="-log10(Adjusted p-value)",x="Go Term",title="Go Term Enrichment for PC2 loadings < -0.05")
p<-p+theme_classic()
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=18,face="bold"))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+theme(plot.title=element_blank())
p
#FIGURE S1G
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Oligo_loadings<-0.05.pdf",height=4,width=10)
print(p)
dev.off()
