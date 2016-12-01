
rm(list=ls())
library(ggplot2)
library(edgeR)

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
allcounts_allcells<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
neurosphere<-read.table("AllCounts_Neurosphere_read_gene_filt_FINAL.txt")

allcells<-cbind(allcounts_allcells,neurosphere)

allcells_col<-c(rep("red",length(allcounts_allcells[1,])),rep("blue",length(neurosphere[1,])))

allcells_geneNum<-colSums(allcounts_allcells!= 0)   #Number of nonzero FPKM values for each cell
neurosphere_geneNum<-colSums(neurosphere!= 0)   #Number of nonzero FPKM values for each cell

d=data.frame(Genes_Detected=allcells_geneNum[2:length(allcells_geneNum)])
e=data.frame(Genes_Detected=neurosphere_geneNum[2:length(neurosphere_geneNum)])
p<-ggplot(d)
p<-p+geom_histogram(aes(x=d$Genes_Detected),smooth=T,binwidth=100,fill="red",alpha=0.6)
p<-p+geom_histogram(data=e,aes(x=e$Genes_Detected),smooth=T,binwidth=100,fill="blue",alpha=0.6)
p<-p+geom_density(aes(x=d$Genes_Detected))
p<-p+labs(x="# of Genes Detected")
p<-p+theme_classic()
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=18))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1B
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("NumGenesDetect_allcells_withneurosphere_miseq.pdf", width = 8, height=4)
print(p)
dev.off()
