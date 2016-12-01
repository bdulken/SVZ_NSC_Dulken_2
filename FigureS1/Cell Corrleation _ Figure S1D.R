
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#The following code generates the FIGURE S1C, the correlation between individual cells and 50 cell pseudopopulations
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
rm(list=ls())
library(ggplot2)
library(edgeR)

# Loading Cell Matrix ------------------------------------------


#This reading in of the cells includes all aNSC, qNSC, NPC, Ast, Young, and Old cells which were expressing at least 500 genes (and for the cells 
#for which it is relevant have less than 25% ERCC spike ins.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
allcounts_allcells<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

#Filtering genes for those that are dected at 10 counts in 5 cells.
greaterthan0<-allcounts_allcells>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_genefilt<-allcounts_allcells[greaterthan0sum>=5,]

#Select qNSC for correlation
miseq_qNSC_2<-allcounts_allcells_genefilt[,grepl("qNSC_B_80_1g",colnames(allcounts_allcells_genefilt))]
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(allcounts_allcells_genefilt)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(allcounts_allcells_genefilt),
                    exonic.size=exonic.gene.sizes.ord)
miseq_qNSC_2_fpkm <- rpkm(miseq_qNSC_2, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

#Select two aNSCs for correlation
miseq_aNSC_1<-allcounts_allcells_genefilt[,grepl("aNSC_B_93",colnames(allcounts_allcells_genefilt))]
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(allcounts_allcells_genefilt)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(allcounts_allcells_genefilt),
                    exonic.size=exonic.gene.sizes.ord)
miseq_aNSC_1_fpkm <- rpkm(miseq_aNSC_1, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

miseq_aNSC_2<-allcounts_allcells_genefilt[,grepl("aNSC_B_92",colnames(allcounts_allcells_genefilt))]
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(allcounts_allcells_genefilt)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(allcounts_allcells_genefilt),
                    exonic.size=exonic.gene.sizes.ord)
miseq_aNSC_2_fpkm <- rpkm(miseq_aNSC_2, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

#Correlation and p values for correlation between the two selected cells.
aNSC.cor<-cor(miseq_aNSC_1_fpkm,miseq_aNSC_2_fpkm,method="pearson")
aNSC.cor_pval<-cor.test(miseq_aNSC_1_fpkm,miseq_aNSC_2_fpkm,method="pearson")

#Plotting scatter of FPKMs for two cells
miseq.data<-data.frame(x=log2(miseq_aNSC_2_fpkm+1),y=log2(miseq_aNSC_1_fpkm+1))
p<-ggplot(miseq.data)
p<-p+geom_point(aes(x=miseq.data$x,y=miseq.data$y))
p<-p+labs(x="log2(FPKM+1)\nSingle aNSC #1 ",y="log2(FPKM+1)\nSingle aNSC #2")
p<-p+theme_classic()
p<-p+theme(axis.text=element_text(size=28))
p<-p+theme(axis.title=element_text(size=28))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1D (left panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Correlation Between aNSC cells - Miseq.pdf",height=7,width=7)
print(p)
dev.off()

#aNSC.cor=0.602

#----------------------------------------------------------------------------------------------------------------------------#
#Outersect function to exclude cells when selecting pseudopopualtions
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}
#----------------------------------------------------------------------------------------------------------------------------#

#Correlation of individual cells in miSeq and HiSeq

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/")
run1<-read.table("AllCounts_HiSeq_run1.txt")
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2/")
run2<-read.table("AllCounts_HiSeq_run2.txt")

run1_NPC<-run1[,grepl("NPC",colnames(run1))]
run1_noNPC<-run1[,!grepl("NPC",colnames(run1))]

comb<-run1_noNPC+run2

comb_all<-cbind(comb,run1_NPC)

comb_all_genefilt<-comb_all[match(rownames(allcounts_allcells_genefilt),rownames(comb_all)),]

hiseq_aNSC_1<-comb_all_genefilt[,grepl("aNSC_B_93",colnames(comb_all_genefilt))]

load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(allcounts_allcells_genefilt)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(allcounts_allcells_genefilt),
                    exonic.size=exonic.gene.sizes.ord)
hiseq_aNSC_1_fpkm <- rpkm(hiseq_aNSC_1, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

#Plotting scatter of FPKMs for two cells
miseq.data<-data.frame(x=log2(hiseq_aNSC_1_fpkm+1),y=log2(miseq_aNSC_1_fpkm+1))
p<-ggplot(miseq.data)
p<-p+geom_point(aes(x=miseq.data$x,y=miseq.data$y))
p<-p+labs(x="log2(FPKM+1)\nSingle aNSC - MiSeq ",y="log2(FPKM+1)\nSingle aNSC - HiSeq")
p<-p+theme_classic()
p<-p+theme(axis.text=element_text(size=28))
p<-p+theme(axis.title=element_text(size=28))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1D (left panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Correlation Between aNSC cell 1 - Miseq vs Hiseq.pdf",height=7,width=7)
print(p)
dev.off()


hiseq_aNSC_2<-comb_all_genefilt[,grepl("aNSC_B_92",colnames(comb_all_genefilt))]

load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(allcounts_allcells_genefilt)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(allcounts_allcells_genefilt),
                    exonic.size=exonic.gene.sizes.ord)
hiseq_aNSC_2_fpkm <- rpkm(hiseq_aNSC_2, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

#Plotting scatter of FPKMs for two cells
miseq.data<-data.frame(x=log2(hiseq_aNSC_2_fpkm+1),y=log2(miseq_aNSC_2_fpkm+1))
p<-ggplot(miseq.data)
p<-p+geom_point(aes(x=miseq.data$x,y=miseq.data$y))
p<-p+labs(x="log2(FPKM+1)\nSingle aNSC #2 - MiSeq ",y="log2(FPKM+1)\nSingle aNSC #2 - HiSeq")
p<-p+theme_classic()
p<-p+theme(axis.text=element_text(size=28))
p<-p+theme(axis.title=element_text(size=28))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1D (left panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Correlation Between aNSC cell 2 - Miseq vs Hiseq.pdf",height=7,width=7)
print(p)
dev.off()

hiseq_qNSC_2<-comb_all_genefilt[,grepl("qNSC_B_80_1g",colnames(comb_all_genefilt))]

load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(allcounts_allcells_genefilt)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(allcounts_allcells_genefilt),
                    exonic.size=exonic.gene.sizes.ord)
hiseq_qNSC_2_fpkm <- rpkm(hiseq_qNSC_2, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

#Plotting scatter of FPKMs for two cells
miseq.data<-data.frame(x=log2(hiseq_qNSC_2_fpkm+1),y=log2(miseq_qNSC_2_fpkm+1))
p<-ggplot(miseq.data)
p<-p+geom_point(aes(x=miseq.data$x,y=miseq.data$y))
p<-p+labs(x="log2(FPKM+1)\nSingle qNSC - MiSeq ",y="log2(FPKM+1)\nSingle qNSC - HiSeq")
p<-p+theme_classic()
p<-p+theme(axis.text=element_text(size=28))
p<-p+theme(axis.title=element_text(size=28))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1D (left panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Correlation Between qNSC cell 1 - Miseq vs Hiseq.pdf",height=7,width=7)
print(p)
dev.off()


#Plotting scatter of FPKMs for two cells
miseq.data<-data.frame(x=log2(hiseq_aNSC_2_fpkm+1),y=log2(hiseq_aNSC_1_fpkm+1))
p<-ggplot(miseq.data)
p<-p+geom_point(aes(x=miseq.data$x,y=miseq.data$y))
p<-p+labs(x="log2(FPKM+1)\nSingle aNSC #1 - HiSeq ",y="log2(FPKM+1)\nSingle aNSC #2 - HiSeq")
p<-p+theme_classic()
p<-p+theme(axis.text=element_text(size=28))
p<-p+theme(axis.title=element_text(size=28))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1D (left panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Correlation Between aNSC cells Hiseq.pdf",height=7,width=7)
print(p)
dev.off()

aNSC.cor.hiseq<-cor(hiseq_aNSC_1_fpkm,hiseq_aNSC_2_fpkm,method="pearson")




#----------------------------------------------------------------------------------------------------------------------------#
#Generating two pseudopopulations of 50 aNSC miseq cells
miseq_PGE<-allcounts_allcells_genefilt[,grep("aNSC",colnames(allcounts_allcells_genefilt))]
miseq_noPGE<-allcounts_allcells_genefilt[,-grep("aNSC",colnames(allcounts_allcells_genefilt))]
miseq_PG<-allcounts_allcells_genefilt[,grep("qNSC",colnames(allcounts_allcells_genefilt))]

#Sample two populations of 50 aNSCs
sample1<-sample(colnames(miseq_PGE),50,replace=F)
miseq_PGE_nosample1<-outersect(sample1,colnames(miseq_PGE))
sample2<-sample(miseq_PGE_nosample1,50,replace=F)

#Get counts for two pseudopopulations
sample1_counts<-allcounts_allcells_genefilt[,match(sample1,colnames(allcounts_allcells_genefilt))]
sample2_counts<-allcounts_allcells_genefilt[,match(sample2,colnames(allcounts_allcells_genefilt))]

#Summing counts for each pseudopopulation and generating FPKM for two pseudopopulations
sample1_comb<-apply(sample1_counts,1,sum)
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(names(sample1_comb)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=names(sample1_comb),
                    exonic.size=exonic.gene.sizes.ord)
sample1_fpkm <- rpkm(sample1_comb, gene.length=genes$exonic.size,
                     normalized.lib.sizes=T, log=FALSE)

sample2_comb<-apply(sample2_counts,1,sum)
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(names(sample2_comb)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=names(sample2_comb),
                    exonic.size=exonic.gene.sizes.ord)
sample2_fpkm <- rpkm(sample2_comb, gene.length=genes$exonic.size,
                     normalized.lib.sizes=T, log=FALSE)

aNSC_pool.cor<-cor(sample1_fpkm,sample2_fpkm,method="pearson")
aNSC_pool.cor_pval<-cor.test(sample1_fpkm,sample2_fpkm,method="pearson")

#plotting
miseq.data<-data.frame(x=log2(sample1_fpkm+1),y=log2(sample2_fpkm+1))
p<-ggplot(miseq.data)
p<-p+geom_point(aes(x=miseq.data$x,y=miseq.data$y))
p<-p+theme_classic()
p<-p+labs(x="log2(FPKM+1)\n50 aggregated aNSCs - #1",y="log2(FPKM+1)\n50 aggregated aNSCs - #2")
p<-p+theme(axis.text=element_text(size=28))
p<-p+theme(axis.title=element_text(size=28))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S1D (right panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("Correlation Between 50 aNSC cells - miseq.pdf",height=7,width=7)
print(p)
dev.off()

