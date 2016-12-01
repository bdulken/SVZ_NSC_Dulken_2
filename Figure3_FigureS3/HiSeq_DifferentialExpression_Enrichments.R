
rm(list=ls())
library(edgeR)
library(monocle)
library(scde)
library("Hmisc")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")

# ====================================================================================================



#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]


setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/")
run1<-read.table("AllCounts_HiSeq_run1.txt")
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2/")
run2<-read.table("AllCounts_HiSeq_run2.txt")

run1_NPC<-run1[,grepl("NPC",colnames(run1))]
run1_noNPC<-run1[,!grepl("NPC",colnames(run1))]

comb<-run1_noNPC+run2

comb_all<-cbind(comb,run1_NPC)

comb_all_cellfilt<-comb_all[,!is.na(match(colnames(comb_all),colnames(allcounts_allcells_nooligo)))]

comb_all_cellfilt_noast<-comb_all_cellfilt[,!grepl("Ast",colnames(comb_all_cellfilt))]
#Filtering for expressed by 5 cells at 10 counts

greaterthan0<-comb_all_cellfilt_noast>10
greaterthan0sum<-rowSums(greaterthan0)
comb_all_cellfilt_noast_genefilt<-comb_all_cellfilt_noast[greaterthan0sum>=5,]


glm <- DGEList(counts=comb_all_cellfilt_noast_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)


#============================================================================================================================================#
#Establishing groups
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])



group1_counts<-comb_all_cellfilt_noast_genefilt[,na.omit(match(group1,colnames(comb_all_cellfilt_noast_genefilt)))]
group2_counts<-comb_all_cellfilt_noast_genefilt[,na.omit(match(group2,colnames(comb_all_cellfilt_noast_genefilt)))]
group3_counts<-comb_all_cellfilt_noast_genefilt[,na.omit(match(group3,colnames(comb_all_cellfilt_noast_genefilt)))]
group4_counts<-comb_all_cellfilt_noast_genefilt[,na.omit(match(group4,colnames(comb_all_cellfilt_noast_genefilt)))]
group5_counts<-comb_all_cellfilt_noast_genefilt[,na.omit(match(group5,colnames(comb_all_cellfilt_noast_genefilt)))]

comb_counts<-cbind(group1_counts,group2_counts,group3_counts,group4_counts,group5_counts)
comb_fac<-factor(c(rep("group1",length(group1_counts[1,])),
                   rep("group2",length(group2_counts[1,])),
                   rep("group3",length(group3_counts[1,])),
                   rep("group4",length(group4_counts[1,])),
                   rep("group5",length(group5_counts[1,]))))

#============================================================================================================================================#
#Setting up SCDE
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output/HiSeqDiff")
n.cores <- 10;
# calculate models
o.ifm <- scde.error.models(counts=comb_counts,groups=comb_fac,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1);


valid.cells <- o.ifm$corr.a >0;
table(valid.cells)
o.ifm <- o.ifm[valid.cells,];

o.prior <- scde.expression.prior(models=o.ifm,counts=comb_counts,length.out=400,show.plot=F)



#Group1 v Group2

comb_fac<-factor(c(rep("group1",length(group1_counts[1,])),rep("group2",length(group2_counts[1,])),rep(NA,length(group3_counts[1,])),rep(NA,length(group4_counts[1,])),rep(NA,length(group5_counts[1,]))))

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group1_group2_ediff_full_10282015_round2_4groups_allgenesINIT_table50_hiseq.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group1_group2_10282015_round2_4groups_allgenesINIT_table50_hiseq_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")



#Group2 vs. Group3

comb_fac<-factor(c(rep(NA,length(group1_counts[1,])),rep("group2",length(group2_counts[1,])),rep("group3",length(group3_counts[1,])),rep(NA,length(group4_counts[1,])),rep(NA,length(group5_counts[1,]))))

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group2_group3_ediff_full_10282015_round2_4groups_allgenesINIT_table50_hiseq.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group2_group3_10282015_round2_4groups_allgenesINIT_table50_hiseq_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")




#Group3 vs. Group4

comb_fac<-factor(c(rep(NA,length(group1_counts[1,])),rep(NA,length(group2_counts[1,])),rep("group3",length(group3_counts[1,])),rep("group4",length(group4_counts[1,])),rep(NA,length(group5_counts[1,]))))

# define two groups of cells

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group3_group4_ediff_full_10282015_round2_4groups_allgenesINIT_table50_hiseq.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group3_group4_10282015_round2_4groups_allgenesINIT_table50_hiseq_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")



#Group4 vs. Group5

comb_fac<-factor(c(rep(NA,length(group1_counts[1,])),rep(NA,length(group2_counts[1,])),rep(NA,length(group3_counts[1,])),rep("group4",length(group4_counts[1,])),rep("group5",length(group5_counts[1,]))))

# define two groups of cells

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group4_group5_ediff_full_10282015_round2_4groups_allgenesINIT_table50_hiseq.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group4_group5_10282015_round2_4groups_allgenesINIT_table50_hiseq_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")









