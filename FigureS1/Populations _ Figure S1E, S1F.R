
#Generation of Figure S1E - PCA with single cell aggregated pseudopopulations and population datasets

rm(list=ls())
library(ggplot2)
library(edgeR)
#Sources special pvclust code capable of running a spearman correlation
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/pvclust_spearman.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/pvclust_internal.R")

#spearman dist function, this is the function that pvclust will call as its distance method
spearman <- function(x, ...) {
  x <- as.matrix(x)
  res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
  res <- as.dist(res)
  attr(res, "method") <- "spearman"
  return(res)
}

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

#Get all cells in particular FACS group
miseq_aNSC<-allcounts_allcells[,grep("aNSC",colnames(allcounts_allcells))]
miseq_qNSC<-allcounts_allcells[,grep("qNSC",colnames(allcounts_allcells))]
miseq_Ast<-allcounts_allcells[,grep("Ast",colnames(allcounts_allcells))]
miseq_NPC<-allcounts_allcells[,grep("NPC",colnames(allcounts_allcells))]

#Load population datasets
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
allcounts_pop<-read.table("AllDenaPopulation_allgenecounts.txt")
allcounts_pop<-allcounts_pop[match(rownames(allcounts_allcells),rownames(allcounts_pop)),]

#Sum all counts for individual cell types
aNSC_sum<-apply(miseq_aNSC,1,sum)
qNSC_sum<-apply(miseq_qNSC,1,sum)
Ast_sum<-apply(miseq_Ast,1,sum)
NPC_sum<-apply(miseq_NPC,1,sum)
allpop_counts<-cbind(aNSC_sum,qNSC_sum,Ast_sum,NPC_sum,allcounts_pop)
colnames(allpop_counts)<-c("aNSC_sing","qNSC_sing","Ast_sing","NPC_sing","Endothelial_1","Endothelial_2","Endothelial_3","Endothelial_4","NPC_1","NPC_2","NPC_3","NPC_4","Ast_1","Ast_2","Ast_3","Ast_4","qNSC_1","qNSC_2","aNSC_1","aNSC_2","aNSC_3")

#Filter genes gased on genes that are detected in at least one of the population datasets OR one of the single cell pseudopopulations
greaterthan0<-allpop_counts>0
greaterthan0sum<-rowSums(greaterthan0)
allpop_counts_filt<-allpop_counts[greaterthan0sum>=1,]

#Normalize with glm normalization
glm <- DGEList(counts=allpop_counts_filt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm <- rpkm(glm.norm, gene.length=genes$exonic.size,
                 normalized.lib.sizes=T, log=FALSE)


#PCA
PCA_int<-prcomp(t(log2(fpkm_glm+1)), scale = T, center = T, retx=T)
summa<-summary(PCA_int)
d<-data.frame(PCA_int$x,factors=colnames(fpkm_glm))
d$factors<-as.character(d$factors)
d$factors <- factor(d$factors, levels=unique(d$factors), ordered = T)
p<-ggplot(d, aes(x=d$PC1, y=d$PC2)) 
p<-p+geom_point(label=d$factors,color=c("#FF3300","#9966FF","#009900","#00CCFF",rep("#FF9900",4),rep("#00CCFF",4),rep("#009900",4),rep("#9966FF",2),rep("#FF3300",3)),shape=c(rep(17,4),rep(16,17)),size=8,alpha=0.8)#color=c(1,2,3,3,3,3,4,4,4,4,5,5,6,6,6))#color=c(1,2,3,3,4,4,4))
p<-p+theme_classic()
p<-p+scale_color_brewer(palette="Dark2") #COLOR BREWER CALL
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 3)-round(summa$importance[3,1],digits = 3))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#Not used for paper
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("qNSC_aNSC_comparingCombinedSingleCellstoPopulation.pdf",width=7,height=7)
print(p)
dev.off()

#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#


#Generation of Figure S1F - Clustering of single cells with pseudopopulations and populations

rm(list=ls())
library(ggplot2)
library(edgeR)

#Sources special pvclust code capable of running a spearman correlation
source("/Volumes/guacamole/Software/R_Files_Packages_Functions/pvclust_spearman.R")
source("/Volumes/guacamole/Software/R_Files_Packages_Functions/pvclust_internal.R")

#spearman dist function, this is the function that pvclust will call as its distance method
spearman <- function(x, ...) {
  x <- as.matrix(x)
  res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
  res <- as.dist(res)
  attr(res, "method") <- "spearman"
  return(res)
}

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

#Pulling out all cells in FACS group
miseq_aNSC<-allcounts_allcells[,grep("aNSC",colnames(allcounts_allcells))]
miseq_qNSC<-allcounts_allcells[,grep("qNSC",colnames(allcounts_allcells))]
miseq_Ast<-allcounts_allcells[,grep("Ast",colnames(allcounts_allcells))]
miseq_NPC<-allcounts_allcells[,grep("NPC",colnames(allcounts_allcells))]

#Reading in Populations
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
allcounts_pop<-read.table("AllDenaPopulation_allgenecounts.txt")
allcounts_pop<-allcounts_pop[match(rownames(allcounts_allcells),rownames(allcounts_pop)),]

#Summing all counts
aNSC_sum<-apply(miseq_aNSC,1,sum)
qNSC_sum<-apply(miseq_qNSC,1,sum)
Ast_sum<-apply(miseq_Ast,1,sum)
NPC_sum<-apply(miseq_NPC,1,sum)
allpop_counts<-cbind(aNSC_sum,qNSC_sum,Ast_sum,NPC_sum,allcounts_pop)
colnames(allpop_counts)<-c("aNSC_sing","qNSC_sing","Ast_sing","NPC_sing","Endothelial_1","Endothelial_2","Endothelial_3","Endothelial_4","NPC_1","NPC_2","NPC_3","NPC_4","Ast_1","Ast_2","Ast_3","Ast_4","qNSC_1","qNSC_2","aNSC_1","aNSC_2","aNSC_3")

#Filter genes gased on genes that are detected in at least one of the population datasets OR one of the single cell pseudopopulations
greaterthan0<-allpop_counts>0
greaterthan0sum<-rowSums(greaterthan0)
allpop_counts_filt<-allpop_counts[greaterthan0sum>=1,]

#Select random cell from each population
rand_qNSC<-miseq_qNSC[,sample(c(1:length(miseq_qNSC[1,])),1)]
rand_Ast<-miseq_Ast[,sample(c(1:length(miseq_Ast[1,])),1)]
rand_aNSC<-miseq_aNSC[,sample(c(1:length(miseq_aNSC[1,])),1)]
rand_NPC<-miseq_NPC[,sample(c(1:length(miseq_NPC[1,])),1)]

rand<-cbind(rand_aNSC,rand_qNSC,rand_Ast,rand_NPC)
rand_filt<-rand[match(rownames(allpop_counts_filt),rownames(allcounts_allcells)),]

allpop_counts_filt_mod<-cbind(rand_filt,allpop_counts_filt)
colnames(allpop_counts_filt_mod)<-c("aNSC_sing","qNSC_sing","Ast_sing","NPC_sing","aNSC_pseudopopulation","qNSC_pseudopopulation","Ast_pseudopopulation",
                                    "NPC_pseudopopulation","Endothelial_1","Endothelial_2","Endothelial_3","Endothelial_4",
                                    "NPC_1","NPC_2","NPC_3","NPC_4","Ast_1","Ast_2","Ast_3","Ast_4","qNSC_1","qNSC_2","aNSC_1","aNSC_2","aNSC_3")

#Perform normalization
glm <- DGEList(counts=allpop_counts_filt_mod)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm <- rpkm(glm.norm, gene.length=genes$exonic.size,
                 normalized.lib.sizes=T, log=FALSE)


#FIGURE S1E
library(Hmisc)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/FigureS1_Output")
pdf("spearman_varclust_withsinglecells.pdf", height=5, width=6)
plot( varclus(log2(fpkm_glm+1), similarity="spearman") , hang = -1)
dev.off()
