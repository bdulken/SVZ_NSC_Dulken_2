#This script generates the correlation plots for Average Time of Expression rankings for supplementary Figure 7
#For this analysis, the FPKM values provided by Llorens-Bobadilla, and the ordering of the hippocampal and
#myoblast datasets by the consensus ordering genes are used.

#Figure S7F,S7G,S7H
# ====================================================================================================

rm(list=ls())
library(edgeR)

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

setwd("/Volumes/guacamole/Analyses/09012015_Figure1")
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
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
pseudotime<-read.table("In50models_pseudotime_mydat.txt")

average_gene_pseudo<-vector()
for(i in 1:length(fpkm_glm_genefilt_nooligo[,1])){
  av_pseudo<-mean(pseudotime[fpkm_glm_genefilt_nooligo[i,]!=0,1])
  average_gene_pseudo<-c(average_gene_pseudo,av_pseudo)
}

names(average_gene_pseudo)<-rownames(fpkm_glm_genefilt_nooligo)
gene_pseudo_sort<-sort(average_gene_pseudo,decreasing=T)








#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#Hippocampus, divided by 10 groups
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
hippo_data<-read.table("Hippocampus_fpkm_withPsuedo.txt")

#Pseudotime for cells ordered by consensus ordered genes.
pseudotime_hippo<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/pseudotime_hippo_table50genes.txt")
hippo_fpkm<-hippo_data[2:length(hippo_data[,1]),]


hippo_fpkm_filt<-hippo_fpkm[match(rownames(fpkm_glm_genefilt_nooligo),rownames(hippo_fpkm)),]

#Ordering genes by average pseudotime of expression in hippocampal dataset
average_gene_pseudo_hippo<-vector()
for(i in 1:length(hippo_fpkm_filt[,1])){
  av_pseudo<-mean(pseudotime_hippo[hippo_fpkm_filt[i,]!=0])
  average_gene_pseudo_hippo<-c(average_gene_pseudo_hippo,av_pseudo)
}

names(average_gene_pseudo_hippo)<-rownames(hippo_fpkm_filt)
gene_pseudo_sort_hippo<-sort(average_gene_pseudo_hippo,decreasing=T)


#========================================================================================================================#
#========================================================================================================================#
#Myoblast, divided by 10 groups
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/datasetHSMM.R")
#myo_fpkm<-exprs(HSMM)
# rownames(myo_fpkm)<-fData(HSMM)[,1]
# myo_fpkm_agg<-aggregate(myo_fpkm,by=list(fData(HSMM)[,1]),FUN=sum)
# myo_fpkm_agg_mod<-myo_fpkm_agg[,2:length(myo_fpkm_agg[1,])]
# rownames(myo_fpkm_agg_mod)<-myo_fpkm_agg[,1]
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")

myo_fpkm_agg_mod<-read.table("fpkm_myoblasts_aggregated_final.txt")
myo_fpkm_agg_mod<-myo_fpkm_agg_mod[,pData(HSMM)$State!=3]

#Pseudotime for cells ordered by consensus ordered genes.
pseudotime_myo<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/pseudotime_myo_table50genes.txt")

myo_fpkm_filt<-myo_fpkm_agg_mod[na.omit(match(tolower(rownames(fpkm_glm_genefilt_nooligo)),tolower(rownames(myo_fpkm_agg_mod)))),]

#Ordering genes by average pseudotime of expression in myoblast dataset
average_gene_pseudo_myo<-vector()
for(i in 1:length(myo_fpkm_filt[,1])){
  av_pseudo<-mean(pseudotime_myo[myo_fpkm_filt[i,]!=0])
  average_gene_pseudo_myo<-c(average_gene_pseudo_myo,av_pseudo)
}

names(average_gene_pseudo_myo)<-rownames(myo_fpkm_filt)
gene_pseudo_sort_myo<-sort(average_gene_pseudo_myo,decreasing=T)

#========================================================================================================================#
#========================================================================================================================#
#Competitor, using their own FPKM

#Competitor, divided by 10 groups
Llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

allcounts_allcells<-Llorens_allcounts
#allcounts_allcells_notaps<-allcounts_allcells[!grepl("tap",colnames(allcounts_allcells))]
allcounts_allcells_noblasts<-allcounts_allcells[!grepl("PSA",colnames(allcounts_allcells))]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("comp_oligos.txt")[,1])
allcounts_allcells_noblasts_nooligo<-allcounts_allcells_noblasts[,-na.omit(match(oligos,colnames(allcounts_allcells_noblasts)))]
allcounts_allcells_noblasts_nooligo_noERCC<-allcounts_allcells_noblasts_nooligo[!grepl("ERCC-",rownames(allcounts_allcells_noblasts_nooligo)),]

split_name<-c()
for(i in 1:length(colnames(allcounts_allcells_noblasts_nooligo_noERCC))){
  split<-strsplit(colnames(allcounts_allcells_noblasts_nooligo_noERCC)[i],split="_STAR")[[1]][1]
  split_name<-c(split_name,split)
}
#Extract only those cells that are not oligodendrocytes or neuroblasts for analysis from the FPKM matrix which they provided
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
comp_data<-read.table("Llorens_aggregated_fpkm.txt")
comp_data_filt<-comp_data[,match(split_name,colnames(comp_data))]

pseudotime_comp<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/pseudotime_comp_theirFPKM_table50genes.txt")

comp_fpkm_filt<-comp_data_filt[na.omit(match(tolower(rownames(fpkm_glm_genefilt_nooligo)),tolower(rownames(comp_data_filt)))),]

#Ordering genes by average pseudotime of expression in Llorens-Bobadilla dataset
average_gene_pseudo_comp<-vector()
for(i in 1:length(comp_fpkm_filt[,1])){
  av_pseudo<-mean(pseudotime_comp[comp_fpkm_filt[i,]!=0,1])
  average_gene_pseudo_comp<-c(average_gene_pseudo_comp,av_pseudo)
}

names(average_gene_pseudo_comp)<-rownames(comp_fpkm_filt)
gene_pseudo_sort_comp<-sort(average_gene_pseudo_comp,decreasing=T)




#Plotting

#========================================================================================================================#
#========================================================================================================================#

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/")

#Filter for genes that are also detected in hippocampal dataset (for comparison to my dataset)
gene_pseudo_sort_hippofilt<-gene_pseudo_sort[!is.na(match(names(gene_pseudo_sort),names(gene_pseudo_sort_hippo)))]

#Generating rank vectors
gene_pseudo_sort_hippofilt_nums<-c(1:length(gene_pseudo_sort_hippofilt))
names(gene_pseudo_sort_hippofilt_nums)<-names(gene_pseudo_sort_hippofilt)

gene_pseudo_sort_hippo_nums<-c(1:length(gene_pseudo_sort_hippo))
names(gene_pseudo_sort_hippo_nums)<-names(gene_pseudo_sort_hippo)

#Ordering the hippocampal rankings in the same order as my dataset for plotting
gene_pseudo_sort_hippo_order<-gene_pseudo_sort_hippo[match(tolower(names(gene_pseudo_sort_hippofilt)),tolower(names(gene_pseudo_sort_hippo)))]
gene_pseudo_sort_hippo_nums_order<-gene_pseudo_sort_hippo_nums[match(tolower(names(gene_pseudo_sort_hippofilt_nums)),tolower(names(gene_pseudo_sort_hippo_nums)))]


png("Pseudotime_Biased_Expression_MyData_consnsusordered_hippo_rank_smoothScatter.png")
smoothScatter(gene_pseudo_sort_hippofilt_nums,gene_pseudo_sort_hippo_nums_order,transformation = function(x) x^.5, nrpoints = 0,xlab="Psuedotime-Biased Gene Ranking - SVZ - This Study",ylab="Psuedotime-Biased Gene Ranking - Hippocampus")
dev.off()

#========================================================================================================================#
#========================================================================================================================#

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/")

#Filter for genes that are also detected in myoblast dataset (for comparison to my dataset)
gene_pseudo_sort_myofilt<-gene_pseudo_sort[!is.na(match(tolower(names(gene_pseudo_sort)),tolower(names(gene_pseudo_sort_myo))))]

#Establishing gene rankings
gene_pseudo_sort_myofilt_nums<-c(1:length(gene_pseudo_sort_myofilt))
names(gene_pseudo_sort_myofilt_nums)<-names(gene_pseudo_sort_myofilt)

gene_pseudo_sort_myo_nums<-c(1:length(gene_pseudo_sort_myo))
names(gene_pseudo_sort_myo_nums)<-names(gene_pseudo_sort_myo)

#Ordering the myoblast rankings in the same order as my dataset for plotting
gene_pseudo_sort_myo_order<-gene_pseudo_sort_myo[match(tolower(names(gene_pseudo_sort_myofilt)),tolower(names(gene_pseudo_sort_myo)))]
gene_pseudo_sort_myo_nums_order<-gene_pseudo_sort_myo_nums[match(tolower(names(gene_pseudo_sort_myofilt_nums)),tolower(names(gene_pseudo_sort_myo_nums)))]

png("Pseudotime_Biased_Expression_MyData_consnsusordered_myo_rank_smoothScatter.png")
smoothScatter(gene_pseudo_sort_myofilt_nums,gene_pseudo_sort_myo_nums_order,transformation = function(x) x^.5, nrpoints = 0,xlab="Psuedotime-Biased Gene Ranking - SVZ - This Study",ylab="Psuedotime-Biased Gene Ranking - Differentiating Myoblast")
dev.off()

#========================================================================================================================#
#========================================================================================================================#
#========================================================================================================================#
#Competitor with their own FPKM
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output/")

#Filter for genes that are also detected in Llorens SVZ dataset (for comparison to my dataset)
gene_pseudo_sort_compfilt<-gene_pseudo_sort[!is.na(match(tolower(names(gene_pseudo_sort)),tolower(names(gene_pseudo_sort_comp))))]

#Establishing gene ranking vector
gene_pseudo_sort_compfilt_nums<-c(1:length(gene_pseudo_sort_compfilt))
names(gene_pseudo_sort_compfilt_nums)<-names(gene_pseudo_sort_compfilt)

gene_pseudo_sort_comp_nums<-c(1:length(gene_pseudo_sort_comp))
names(gene_pseudo_sort_comp_nums)<-names(gene_pseudo_sort_comp)

#Ordering the Llorens-Bobadilla rankings in the same order as my dataset for plotting
gene_pseudo_sort_comp_order<-gene_pseudo_sort_comp[match(tolower(names(gene_pseudo_sort_compfilt)),tolower(names(gene_pseudo_sort_comp)))]
gene_pseudo_sort_comp_nums_order<-gene_pseudo_sort_comp_nums[match(tolower(names(gene_pseudo_sort_compfilt_nums)),tolower(names(gene_pseudo_sort_comp_nums)))]


png("Pseudotime_Biased_Expression_MyData_LlorensFPKM_comp_rank_smoothScatter.png")
smoothScatter(gene_pseudo_sort_compfilt_nums,gene_pseudo_sort_comp_nums_order,transformation = function(x) x^.5, nrpoints = 0,xlab="Psuedotime-Biased Gene Ranking - SVZ - This Study",ylab="Psuedotime-Biased Gene Ranking - SVZ - Llorens-Bobadilla, et al.")
dev.off()






#Calculation of correlations for these comparisons
corrs<-c(cor(gene_pseudo_sort_compfilt_nums,gene_pseudo_sort_comp_nums_order,method="spearman"),
         cor(gene_pseudo_sort_hippofilt_nums,gene_pseudo_sort_hippo_nums_order,method="spearman"),
         cor(gene_pseudo_sort_myofilt_nums,gene_pseudo_sort_myo_nums_order,method="spearman"))

corr.test.llorens<-cor.test(gene_pseudo_sort_compfilt_nums,gene_pseudo_sort_comp_nums_order,method="spearman")
corr.test.hippo<-cor.test(gene_pseudo_sort_hippofilt_nums,gene_pseudo_sort_hippo_nums_order,method="spearman")
corr.test.myo<-cor.test(gene_pseudo_sort_myofilt_nums,gene_pseudo_sort_myo_nums_order,method="spearman")

write.table(corrs,"correlation_all_compOTHERFPKM_hippoCONSENSUS_myoCONSENSUS.txt")


