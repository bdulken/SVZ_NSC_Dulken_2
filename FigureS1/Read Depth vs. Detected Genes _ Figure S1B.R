
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#The following script generates the plot comparing the gene counts to read depth as found in FIGURE S1B.
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#

#Reading in the counts for all cells that were sequenced on miseq and on hiseq and removing one aNSC cell which did not have
#the threshold level of 500 genes detected on MiSeq.

#The hiseq counts which are being read in have already been combined for the indicated runs.
rm(list=ls())
library(ggplot2)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
miseq<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")
miseq<-miseq[,!grepl("NPC",colnames(miseq))]

hi1<-read.table("AllCounts_HiSeq_Run1_highqual.txt")
hi1<-hi1[,!grepl("NPC",colnames(hi1))]

hi2<-read.table("AllCounts_HiSeq_Run2_highqual.txt")

#Pulling out cells that are in both miseq and hiseq1
miseq<-miseq[,match(colnames(hi1),colnames(miseq))]

miseq_1hi<-miseq+hi1
miseq_2hi<-miseq+hi1+hi2

#Separating miseq based on cell type
miseq_Ast<-miseq[,grep("Ast",colnames(miseq))]
miseq_aNSC<-miseq[,grep("aNSC",colnames(miseq))]
miseq_mod<-miseq[,-grep("aNSC",colnames(miseq))]
miseq_qNSC<-miseq_mod[,grep("qNSC",colnames(miseq_mod))]

#Separating hiseq1+miseq based on cell type
miseq_1hi_Ast<-miseq_1hi[,grep("Ast",colnames(miseq_1hi))]
miseq_1hi_aNSC<-miseq_1hi[,grep("aNSC",colnames(miseq_1hi))]
miseq_1hi_mod<-miseq_1hi[,-grep("aNSC",colnames(miseq_1hi))]
miseq_1hi_qNSC<-miseq_1hi_mod[,grep("qNSC",colnames(miseq_1hi_mod))]

#Separating hiseq1+hiseq2+miseq based on cell type
miseq_2hi_Ast<-miseq_2hi[,grep("Ast",colnames(miseq_2hi))]
miseq_2hi_aNSC<-miseq_2hi[,grep("aNSC",colnames(miseq_2hi))]
miseq_2hi_mod<-miseq_2hi[,-grep("aNSC",colnames(miseq_2hi))]
miseq_2hi_qNSC<-miseq_2hi_mod[,grep("qNSC",colnames(miseq_2hi_mod))]


#Determining Read Counts
miseq_Ast_readcounts<-apply(miseq_Ast[,1:length(miseq_Ast[1,])],2,sum)
miseq_qNSC_readcounts<-apply(miseq_qNSC[,1:length(miseq_qNSC[1,])],2,sum)
miseq_aNSC_readcounts<-apply(miseq_aNSC[,1:length(miseq_aNSC[1,])],2,sum)

miseq_1hi_Ast_readcounts<-apply(miseq_1hi_Ast[,1:length(miseq_1hi_Ast[1,])],2,sum)
miseq_1hi_qNSC_readcounts<-apply(miseq_1hi_qNSC[,1:length(miseq_1hi_qNSC[1,])],2,sum)
miseq_1hi_aNSC_readcounts<-apply(miseq_1hi_aNSC[,1:length(miseq_1hi_aNSC[1,])],2,sum)

miseq_2hi_Ast_readcounts<-apply(miseq_2hi_Ast[,1:length(miseq_2hi_Ast[1,])],2,sum)
miseq_2hi_qNSC_readcounts<-apply(miseq_2hi_qNSC[,1:length(miseq_2hi_qNSC[1,])],2,sum)
miseq_2hi_aNSC_readcounts<-apply(miseq_2hi_aNSC[,1:length(miseq_2hi_aNSC[1,])],2,sum)

#Determining gene counts
miseq_Ast_genecounts<-colSums(miseq_Ast != 0)
miseq_qNSC_genecounts<-colSums(miseq_qNSC != 0)
miseq_aNSC_genecounts<-colSums(miseq_aNSC != 0)

miseq_1hi_Ast_genecounts<-colSums(miseq_1hi_Ast != 0)
miseq_1hi_qNSC_genecounts<-colSums(miseq_1hi_qNSC != 0)
miseq_1hi_aNSC_genecounts<-colSums(miseq_1hi_aNSC != 0)

miseq_2hi_Ast_genecounts<-colSums(miseq_2hi_Ast != 0)
miseq_2hi_qNSC_genecounts<-colSums(miseq_2hi_qNSC != 0)
miseq_2hi_aNSC_genecounts<-colSums(miseq_2hi_aNSC != 0)


#Plotting correlation - aNSC - log10
readcounts_aNSC_total<-log10(c(miseq_aNSC_readcounts,miseq_1hi_aNSC_readcounts,miseq_2hi_aNSC_readcounts)+1)
genecounts_aNSC_total<-log10(c(miseq_aNSC_genecounts,miseq_1hi_aNSC_genecounts,miseq_2hi_aNSC_genecounts)+1)
cell_ident<-c(rep(seq(1,length(miseq_aNSC_readcounts)),3))
data=data.frame(readcounts=readcounts_aNSC_total,genecounts=genecounts_aNSC_total,cellident=cell_ident)
p<-ggplot(data)
p<-p+geom_line(aes(x=data$readcounts,y=data$genecounts,group=data$cellident,col=data$cellident))
p<-p+geom_point(aes(x=data$readcounts,y=data$genecounts,group=data$cellident))
p<-p+theme_classic()
p<- p + theme(legend.position = "none")
p<-p+labs(x="log10(Total Reads + 1)",y="log10(#Genes Detected + 1)")
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p<-p+scale_y_continuous(limits=c(2,4))
p
#FIGURE S1B
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Depth of Sequencing Correlation - aNSC - All three runs_log_nozeros.pdf",height=5,width=10)
print(p)
dev.off()


