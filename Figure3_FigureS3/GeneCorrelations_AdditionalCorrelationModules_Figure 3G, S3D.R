
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

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus.txt")[,1])

#Extract fpkm for each group of cells individually
group1_counts<-fpkm_glm_genefilt_nooligo[,match(group1,colnames(fpkm_glm_genefilt_nooligo))]
group2_counts<-fpkm_glm_genefilt_nooligo[,match(group2,colnames(fpkm_glm_genefilt_nooligo))]
group3_counts<-fpkm_glm_genefilt_nooligo[,match(group3,colnames(fpkm_glm_genefilt_nooligo))]
group4_counts<-fpkm_glm_genefilt_nooligo[,match(group4,colnames(fpkm_glm_genefilt_nooligo))]
group5_counts<-fpkm_glm_genefilt_nooligo[,match(group5,colnames(fpkm_glm_genefilt_nooligo))]

#Use only aNSC-mid and aNSC-late cells
comb_fpkm<-cbind(group3_counts,group4_counts)

int_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models.txt")[,1])
#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_consensus50genes_group3_group4.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()

int_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/02012016_names_in25models.txt")[,1])
#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[na.omit(match(int_genes,rownames(comb_fpkm))),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_consensus25genes_group3_group4.pdf",height=15,width=15)
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()


int_genes<-c("Cdk1","Cdk4","Tpx2","Prc1","Ccna2","Ccnb1","Ccnd2","Ube2c","Nono","Atp1a2","Ntsr2","Gja1","Jag1","Fgfr3","Dlx1","Dlx2")
#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_cellcyclemodule_group3_group4.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()


#All Genes
comb_fpkm_curr<-comb_fpkm

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_allgenes_group3_group4.pdf",height=100,width=100)
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()









comb_fpkm<-cbind(group1_counts,
                 group2_counts,
                 group3_counts,
                 group4_counts,
                 group5_counts)

int_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models.txt")[,1])
#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_consensus50genes_allcells.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()

int_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/02012016_names_in25models.txt")[,1])
#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[na.omit(match(int_genes,rownames(comb_fpkm))),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_consensus25genes_allcells.pdf",height=15,width=15)
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()


int_genes<-c("Cdk1","Cdk4","Tpx2","Prc1","Ccna2","Ccnb1","Ccnd2","Ube2c","Nono","Atp1a2","Ntsr2","Gja1","Jag1","Fgfr3","Dlx1","Dlx2")
#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[na.omit(match(int_genes,rownames(comb_fpkm))),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_cellcyclemodule_allcells.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()


#All Genes
comb_fpkm_curr<-fpkm_glm_genefilt_nooligo

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_allgenes_allcells.pdf",height=100,width=100)
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()















