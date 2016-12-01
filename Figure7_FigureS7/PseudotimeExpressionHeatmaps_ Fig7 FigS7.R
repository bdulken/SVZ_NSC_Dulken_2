#This script contains the code for generating pseudotime spectra plots for all conditions.
#Figure 7B

rm(list=ls())
library(edgeR)
library(reshape)

# ====================================================================================================


#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]

#remove astrocytes for regression learning
allcounts_allcells_nooligo_noast<-allcounts_allcells_nooligo[,!grepl("Ast",colnames(allcounts_allcells_nooligo))]

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_allcells_nooligo_noast>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_nooligo_noast_genefilt<-allcounts_allcells_nooligo_noast[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_nooligo_noast_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)



# ====================================================================================================
# ====================================================================================================
# ====================================================================================================
# ====================================================================================================
#My data divided by 10 groups

#Read in pseudotime as defined by ordering by consensus ordering genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
pseudotime<-read.table("In50models_pseudotime_mydat.txt")

#Calculate "Average pseudotime of Expression" for each gene
average_gene_pseudo<-vector()
for(i in 1:length(fpkm_glm_genefilt_nooligo[,1])){
  av_pseudo<-mean(pseudotime[fpkm_glm_genefilt_nooligo[i,]!=0,1])
  average_gene_pseudo<-c(average_gene_pseudo,av_pseudo)
}

#Sort genes based on Average Time of Expression
names(average_gene_pseudo)<-rownames(fpkm_glm_genefilt_nooligo)
gene_pseudo_sort<-sort(average_gene_pseudo,decreasing=T)

#Divide pseudotime spectrum into 10 equally sized groups
fpkm_glm_corrsort<-fpkm_glm_genefilt_nooligo[na.omit(match(names(gene_pseudo_sort),rownames(fpkm_glm_genefilt_nooligo))),]
pseudotime_factor_5<-cut_number(pseudotime[,1],10,labels=F)
fpkm_glm_pseudo<-data.frame(t(fpkm_glm_genefilt_nooligo),pseudo=pseudotime_factor_5)

#Find mean expression of each gene for each group of cells
pseudo_agg<-aggregate(t(fpkm_glm_corrsort), by=list(pseudotime_factor_5) ,FUN=mean)
pseudo_agg_mat<-as.matrix(pseudo_agg[,2:length(pseudo_agg[1,])],type="numeric")
pseudo_agg_mat_log<-log2(pseudo_agg_mat+1)

#Establish subsets of the genes ordered by Average pseudotime of Expression for plotting (Last 20)
pseudo_agg_mat_log_extextbrev<-cbind(pseudo_agg_mat_log[,1:20])

#Figure 7B (left panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
pseudo_agg_mat_log_re<-apply(pseudo_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("MyData_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()

#Figure 7C (left panel)
library(scales)
pseudo_agg_mat_log_re<-apply(pseudo_agg_mat_log_extextbrev,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("MyData_pseudotimegrouped_averagepospseudotime_mean_extextbrev.pdf",height=3.5,width=3.5)
print(p)
dev.off()

#========================================================================================================================#
#========================================================================================================================#

#Figure 7B, Figure 7C
#SVZ, Llorens-Bobadilla, divided by 10 groups
Llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

#Remove neuroblasts from analysis
allcounts_allcells<-Llorens_allcounts
#allcounts_allcells_notaps<-allcounts_allcells[!grepl("tap",colnames(allcounts_allcells))]
allcounts_allcells_noblasts<-allcounts_allcells[!grepl("PSA",colnames(allcounts_allcells))]

#Remove oligodendrocytes from analysis
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
oligos<-as.vector(read.table("comp_oligos.txt")[,1])
allcounts_allcells_noblasts_nooligo<-allcounts_allcells_noblasts[,-na.omit(match(oligos,colnames(allcounts_allcells_noblasts)))]
allcounts_allcells_noblasts_nooligo_noERCC<-allcounts_allcells_noblasts_nooligo[!grepl("ERCC-",rownames(allcounts_allcells_noblasts_nooligo)),]

#Filtering for genes detected as expressed in our dataset
allcounts_allcells_noblasts_nooligo_noERCC_genefilt<-allcounts_allcells_noblasts_nooligo_noERCC[match(rownames(fpkm_glm_genefilt_nooligo),rownames(allcounts_allcells_noblasts_nooligo_noERCC)),]

glm <- DGEList(counts=allcounts_allcells_noblasts_nooligo_noERCC_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_comp_cells <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                     normalized.lib.sizes=T, log=F)

#Read in pseudotime of Llorens-Bobadilla data ordered by consensus ordering genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
pseudotime_comp<-read.table("pseudotime_comp_table50genes.txt")

#Order genes by Average Pseudotime of expression from our dataset
fpkm_glm_corrsort_comp<-fpkm_glm_genefilt_comp_cells[na.omit(match(names(gene_pseudo_sort),rownames(fpkm_glm_genefilt_comp_cells))),]
#Divide into 10 groups
pseudotime_factor_comp<-cut_number(pseudotime_comp[,1],10,labels=F)

#Aggregate by group
pseudo_comp_agg<-aggregate(t(fpkm_glm_corrsort_comp), by=list(pseudotime_factor_comp) ,FUN=mean)
pseudo_comp_agg_mat<-as.matrix(pseudo_comp_agg[,2:length(pseudo_comp_agg[1,])],type="numeric")

pseudo_comp_agg_mat_log<-log2(pseudo_comp_agg_mat+1)

pseudo_comp_agg_mat_log_extextbrev<-cbind(pseudo_comp_agg_mat_log[,1:20])

#Figure 7B (right panel)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
pseudo_comp_agg_mat_log_re<-apply(pseudo_comp_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_comp_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("Llorens_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()

#Figure 7C (right panel)
library(scales)
pseudo_agg_mat_log_re<-apply(pseudo_comp_agg_mat_log_extextbrev,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("Llorens_pseudotimegrouped_averagepospseudotime_mean_extextbrev.pdf",height=3.5,width=3.5)
print(p)
dev.off()


#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#Hippocampus, divided by 10 groups
#Using pseudotime values provided by authors
#Figure S7D

#Read in hippocampal data as provided by the authors, pseudotime provided by authors
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
hippo_data<-read.table("Hippocampus_fpkm_withPsuedo.txt")
pseudotime_hippo<-as.numeric(hippo_data[1,])
hippo_fpkm<-hippo_data[2:length(hippo_data[,1]),]

#Sort genes by average pseudotime of expression from our data
hippo_corrsort<-hippo_fpkm[na.omit(match(names(gene_pseudo_sort),rownames(hippo_fpkm))),]
pseudotime_factor_hippo<-cut_number(pseudotime_hippo,10,labels=F)

#Aggregate into 10 groups
pseudo_hippo_agg<-aggregate(t(hippo_corrsort), by=list(pseudotime_factor_hippo) ,FUN=mean)
pseudo_hippo_agg_mat<-as.matrix(pseudo_hippo_agg[,2:length(pseudo_hippo_agg[1,])],type="numeric")

pseudo_hippo_agg_mat_log<-log2(pseudo_hippo_agg_mat+1)


setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
pseudo_hippo_agg_mat_log_re<-apply(pseudo_hippo_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_hippo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("hippocampus_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()




#========================================================================================================================#
#========================================================================================================================#
#Myoblast, divided by 10 groups
#Using pseudotime values provided by authors
#Figure S7A
#Read in myoblast data as provided by the authors
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/datasetHSMM.R")
#myo_fpkm<-exprs(HSMM)
# rownames(myo_fpkm)<-fData(HSMM)[,1]
# myo_fpkm_agg<-aggregate(myo_fpkm,by=list(fData(HSMM)[,1]),FUN=sum)
# myo_fpkm_agg_mod<-myo_fpkm_agg[,2:length(myo_fpkm_agg[1,])]
# rownames(myo_fpkm_agg_mod)<-myo_fpkm_agg[,1]
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")

myo_fpkm_agg_mod<-read.table("fpkm_myoblasts_aggregated_final.txt")
myo_fpkm_agg_mod<-myo_fpkm_agg_mod[,pData(HSMM)$State!=3]

pseudotime_myo<-pData(HSMM)[,6] #Provided by authors
pseudotime_myo<-pseudotime_myo[pData(HSMM)$State!=3]

#Filter for genes detected in our dataset and order by average pseudotime of expression in our dataset
myo_corrsort<-myo_fpkm_agg_mod[na.omit(match(tolower(names(gene_pseudo_sort)),tolower(rownames(myo_fpkm_agg_mod)))),]
pseudotime_factor_myo<-cut_number(pseudotime_myo,10,labels=F)

pseudo_myo_agg<-aggregate(t(myo_corrsort), by=list(pseudotime_factor_myo) ,FUN=mean)
pseudo_myo_agg_mat<-as.matrix(pseudo_myo_agg[,2:length(pseudo_myo_agg[1,])],type="numeric")


pseudo_myo_agg_mat_log<-log2(pseudo_myo_agg_mat+1)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
library(reshape)
pseudo_myo_agg_mat_log_re<-apply(pseudo_myo_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_myo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
#This reversed the pseudotime order, so that the proliferative cells were ordered on the right side of the plot.
p<-p+scale_x_reverse()
pdf("myoblast_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()







#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#

#Ordering based on consensus ordering genes for control datasets


#========================================================================================================================#
#========================================================================================================================#
#Llorens, divided by 10 groups, FPKM provided by llorens-bobadilla et al
# Figure S7B


#The following chunk of code is simply to get the cells that would be used for the subsequent analyses (eliminating oligodendrocytes
#And neuroblasts)
#========================================================================================================================#
Llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

allcounts_allcells<-Llorens_allcounts
#allcounts_allcells_notaps<-allcounts_allcells[!grepl("tap",colnames(allcounts_allcells))]
allcounts_allcells_noblasts<-allcounts_allcells[!grepl("PSA",colnames(allcounts_allcells))]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
oligos<-as.vector(read.table("comp_oligos.txt")[,1])
allcounts_allcells_noblasts_nooligo<-allcounts_allcells_noblasts[,-na.omit(match(oligos,colnames(allcounts_allcells_noblasts)))]
allcounts_allcells_noblasts_nooligo_noERCC<-allcounts_allcells_noblasts_nooligo[!grepl("ERCC-",rownames(allcounts_allcells_noblasts_nooligo)),]
#========================================================================================================================#
split_name<-c()
for(i in 1:length(colnames(allcounts_allcells_noblasts_nooligo_noERCC))){
  split<-strsplit(colnames(allcounts_allcells_noblasts_nooligo_noERCC)[i],split="_STAR")[[1]][1]
  split_name<-c(split_name,split)
}
#Extract only those cells that are not oligodendrocytes or neuroblasts for analysis from the FPKM matrix which they provided
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
comp_data<-read.table("Llorens_aggregated_fpkm.txt")
comp_data_filt<-comp_data[,match(split_name,colnames(comp_data))]


fpkm_glm_genefilt_comp_cells<-comp_data_filt

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
pseudotime_comp<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/pseudotime_comp_theirFPKM_table50genes.txt")

#Filter for genes detected in our analysis and order by pseudotime of expression in our dataset
fpkm_glm_corrsort_comp<-fpkm_glm_genefilt_comp_cells[na.omit(match(names(gene_pseudo_sort),rownames(fpkm_glm_genefilt_comp_cells))),]
pseudotime_factor_comp<-cut_number(pseudotime_comp[,1],10,labels=F)

#Aggregate expression by group
pseudo_comp_agg<-aggregate(t(fpkm_glm_corrsort_comp), by=list(pseudotime_factor_comp) ,FUN=mean)
pseudo_comp_agg_mat<-as.matrix(pseudo_comp_agg[,2:length(pseudo_comp_agg[1,])],type="numeric")

pseudo_comp_agg_mat_log<-log2(pseudo_comp_agg_mat+1)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
pseudo_comp_agg_mat_log_re<-apply(pseudo_comp_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_comp_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("Llorens_consensusorderinggenes_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()


#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#Hippocampus, divided by 10 groups, ordered using consensus ordering genes
# Figure S7E
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
hippo_data<-read.table("Hippocampus_fpkm_withPsuedo.txt")

#Use pseudotime for hippocampus generated by ordering cells by consensus ordering genes.
pseudotime_hippo<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/pseudotime_hippo_table50genes.txt")
hippo_fpkm<-hippo_data[2:length(hippo_data[,1]),]

#Filter for genes detected in our analysis and order by pseudotime of expressin in our dataset
hippo_corrsort<-hippo_fpkm[na.omit(match(names(gene_pseudo_sort),rownames(hippo_fpkm))),]
pseudotime_factor_hippo<-cut_number(pseudotime_hippo[,1],10,labels=F)

pseudo_hippo_agg<-aggregate(t(hippo_corrsort), by=list(pseudotime_factor_hippo) ,FUN=mean)
pseudo_hippo_agg_mat<-as.matrix(pseudo_hippo_agg[,2:length(pseudo_hippo_agg[1,])],type="numeric")


pseudo_hippo_agg_mat_log<-log2(pseudo_hippo_agg_mat+1)


setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
pseudo_hippo_agg_mat_log_re<-apply(pseudo_hippo_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_hippo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
pdf("hippocampus_consensusorderinggenes_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()


#========================================================================================================================#
#========================================================================================================================#
#Myoblast, divided by 10 groups, ordered using consensus ordering genes
# Figure S7C
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/datasetHSMM.R")
#myo_fpkm<-exprs(HSMM)
# rownames(myo_fpkm)<-fData(HSMM)[,1]
# myo_fpkm_agg<-aggregate(myo_fpkm,by=list(fData(HSMM)[,1]),FUN=sum)
# myo_fpkm_agg_mod<-myo_fpkm_agg[,2:length(myo_fpkm_agg[1,])]
# rownames(myo_fpkm_agg_mod)<-myo_fpkm_agg[,1]
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")

myo_fpkm_agg_mod<-read.table("fpkm_myoblasts_aggregated_final.txt")
myo_fpkm_agg_mod<-myo_fpkm_agg_mod[,pData(HSMM)$State!=3]

#Use pseudotime values obtained by ordering myoblast cells by consensus ordering genes
pseudotime_myo<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/pseudotime_myo_table50genes.txt")

#Filter for genes detected in our analysis and order by pseudotime of expressin in our dataset
myo_corrsort<-myo_fpkm_agg_mod[na.omit(match(tolower(names(gene_pseudo_sort)),tolower(rownames(myo_fpkm_agg_mod)))),]
pseudotime_factor_myo<-cut_number(pseudotime_myo[,1],10,labels=F)

pseudo_myo_agg<-aggregate(t(myo_corrsort), by=list(pseudotime_factor_myo) ,FUN=mean)
pseudo_myo_agg_mat<-as.matrix(pseudo_myo_agg[,2:length(pseudo_myo_agg[1,])],type="numeric")


pseudo_myo_agg_mat_log<-log2(pseudo_myo_agg_mat+1)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
library(scales)
library(reshape)
pseudo_myo_agg_mat_log_re<-apply(pseudo_myo_agg_mat_log,2,rescale)
#ggplot2 plotted version
data<-melt(pseudo_myo_agg_mat_log_re)
data$X2<-as.character(data$X2)
data$X2 <- factor(data$X2, levels=unique(data$X2), ordered = T)
p <- ggplot(data, aes(X1, X2))
p<-p+geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "tomato")
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank())
p<-p + theme(  axis.text.x = element_blank(),  axis.text.y = element_blank(),  axis.ticks = element_blank())
p<-p + labs(x="Relative Pseudotime", y="Genes")
p<-p+scale_x_reverse()
pdf("myoblast_consensusorderinggenes_pseudotimegrouped_averagepospseudotime_mean.pdf",height=60,width=5)
print(p)
dev.off()




