#The following code performs Monocle ordering for Hippocampal and myoblast datasets using the consensus ordering genes
#This is for generating the supplemental figures in which ordering is done using these genes.

rm(list=ls())
library(monocle)
library(ggplot2)
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")


#Hippocampus
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
hippo_data<-read.table("Hippocampus_fpkm_withPsuedo.txt")
pseudotime_hippo<-as.numeric(hippo_data[1,])
hippo_fpkm<-hippo_data[2:length(hippo_data[,1]),]

comb_fpkm_comp<-hippo_fpkm
#Monocle ordering with llorenss cells
comb_fpkm_comp_col<-c(rep("Hippo",length(comb_fpkm_comp[1,])))

#This pheno vector can be referenced by calling "color_by='type'" in any of the plotting 
#functions in monocle to color the cells plotted by their type.
pheno<-comb_fpkm_comp_col

pheno.data.df <- data.frame(type=pheno)
rownames(pheno.data.df) <- colnames(comb_fpkm_comp)

feature<-data.frame(rownames(comb_fpkm_comp))


pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm_comp, phenoData = pd)


#Ordering using the consensus ordering genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
curated_ordering_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])
match(curated_ordering_genes,rownames(comb_fpkm_comp))
HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = F) # Order cells

pdf(paste("spanningtree_looped_hippo_table50.pdf",sep=""))
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#FF9900","#6600CC"))
print(p)
dev.off()
library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

write.table(pseudotime_mod,"pseudotime_hippo_table50genes.txt")


# ====================================================================================================
# ====================================================================================================
# ====================================================================================================
# ====================================================================================================

#Myoblasts
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/datasetHSMM.R")
#myo_fpkm<-exprs(HSMM)
# rownames(myo_fpkm)<-fData(HSMM)[,1]
# myo_fpkm_agg<-aggregate(myo_fpkm,by=list(fData(HSMM)[,1]),FUN=sum)
# myo_fpkm_agg_mod<-myo_fpkm_agg[,2:length(myo_fpkm_agg[1,])]
# rownames(myo_fpkm_agg_mod)<-myo_fpkm_agg[,1]
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")

myo_fpkm_agg_mod<-read.table("fpkm_myoblasts_aggregated_final.txt")
myo_fpkm_agg_mod_noint<-myo_fpkm_agg_mod[,pData(HSMM)$State!=3]

myo_fpkm<-myo_fpkm_agg_mod_noint

comb_fpkm_comp<-myo_fpkm
#Monocle ordering with llorenss cells
comb_fpkm_comp_col<-c(rep("myo",length(comb_fpkm_comp[1,])))

#This pheno vector can be referenced by calling "color_by='type'" in any of the plotting 
#functions in monocle to color the cells plotted by their type.
pheno<-comb_fpkm_comp_col

pheno.data.df <- data.frame(type=pheno)
rownames(pheno.data.df) <- colnames(comb_fpkm_comp)

feature<-data.frame(rownames(comb_fpkm_comp))


pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm_comp, phenoData = pd)


#This is the list of ordering genes which is used to order the qNSCs, all have known roles in the differentiation of neurons,
#and because we have the qNSCs and Astrocytes in this as well I will also use some genes associated with quiescence
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
#I use a different file for the consensus ordering genes because the names of the genes were different for the two annotations
curated_ordering_genes<-as.vector(read.table("10282015_names_in50models_formyo_3.txt")[,1])
match(toupper(curated_ordering_genes),rownames(comb_fpkm_comp))
HSMM_q <- setOrderingFilter(HSMM_q, toupper(curated_ordering_genes)) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = F) # Order cells

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
pdf(paste("spanningtree_looped_myo_table50.pdf",sep=""))
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#FF9900","#6600CC"))
print(p)
dev.off()
library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

write.table(pseudotime_mod,"pseudotime_myo_table50genes.txt")

my_genes<-toupper(c("Id3","Rpl32","Cdk4","Cdk1"))
#my_genes<-c("Dlx1","Dlx2","Dcx","Notch2","Notch1","Malat1","Jag1","Fgfr3","Apoe","Tubb3","Ascl1","Dlx5")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("aggregated_ordering_top100_mean_myo.pdf",height=6,width=5)
plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=1.5,color_custom=c("#FF9900","#6600CC"),trend_formula=NULL)
dev.off()







#Llorens, using their own FPKM

#Llorens, divided by 10 groups
llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

allcounts_allcells<-llorens_allcounts
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

comb_fpkm_comp<-comp_data_filt
comb_fpkm_comp_col<-c(rep("NSC",length(comb_fpkm_comp[1,])))
comb_fpkm_comp_col[grepl("tap",colnames(comb_fpkm_comp))]<-"TAP"

#This pheno vector can be referenced by calling "color_by='type'" in any of the plotting 
#functions in monocle to color the cells plotted by their type.
pheno<-comb_fpkm_comp_col

pheno.data.df <- data.frame(type=pheno)
rownames(pheno.data.df) <- colnames(comb_fpkm_comp)

feature<-data.frame(rownames(comb_fpkm_comp))


pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm_comp, phenoData = pd)


#Ordering using the consensus ordering genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
#I use a different file for the consensus ordering genes because the names of some genes were different for the two annotations
curated_ordering_genes<-as.vector(read.table("10282015_names_in50models_forcomp_3.txt")[,1])
match(curated_ordering_genes,rownames(comb_fpkm_comp))
HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = F) # Order cells

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure7_Output")
pdf(paste("spanningtree_looped_comp_theirFPKM_table50.pdf",sep=""))
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#FF9900","#6600CC"))
print(p)
dev.off()
library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

write.table(pseudotime_mod,"pseudotime_comp_theirFPKM_table50genes.txt")

