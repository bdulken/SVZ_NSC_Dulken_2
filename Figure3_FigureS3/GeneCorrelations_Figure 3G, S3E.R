#This script calculates correlations between individual genes to demonstrate the mutual exclusivity
#of expression of mediators of neurogenesis and astrocytic markers/self-renewal genes

#Figures 3G, S3D

rm(list=ls())
library(edgeR)
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
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])


#Extract fpkm for each group of cells individually
group1_counts<-fpkm_glm_genefilt_nooligo[,match(group1,colnames(fpkm_glm_genefilt_nooligo))]
group2_counts<-fpkm_glm_genefilt_nooligo[,match(group2,colnames(fpkm_glm_genefilt_nooligo))]
group3_counts<-fpkm_glm_genefilt_nooligo[,match(group3,colnames(fpkm_glm_genefilt_nooligo))]
group4_counts<-fpkm_glm_genefilt_nooligo[,match(group4,colnames(fpkm_glm_genefilt_nooligo))]
group5_counts<-fpkm_glm_genefilt_nooligo[,match(group5,colnames(fpkm_glm_genefilt_nooligo))]

#Use only aNSC-mid and aNSC-late cells
comb_fpkm<-cbind(group3_counts,group4_counts)

#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[match(c("Atp1a2","Ntsr2","Gja1","Jag1","Fgfr3","Dlx1","Dlx2"),rownames(comb_fpkm)),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()








comb_fpkm<-cbind(group3_counts,group4_counts)


correlation_all<-rcorr(log2(t(comb_fpkm+1)),type='spearman')
correlation_all_ast<-correlation_all$r[match(c("Atp1a2","Ntsr2","Gja1","Jag1","Fgfr3"),rownames(correlation_all$r)),]

correlation_ast_sums<-sort(colSums(correlation_all_ast),decreasing=T)

correlation_ast_top<-names(correlation_ast_sums[1:15])

correlation_all_dlx<-correlation_all$r[match(c("Dlx1","Dlx2"),rownames(correlation_all$r)),]

correlation_dlx_sums<-sort(colSums(correlation_all_dlx),decreasing=T)

correlation_dlx_top<-names(correlation_dlx_sums[1:12])







int_genes<-c(correlation_dlx_top,correlation_ast_top)

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
pdf("Correlation_select_genes_heatmap_2_cluster_expandedcluster.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()





# ====================================================================================================
#Repeat analysis with HiSeq cells only.
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



#Use only aNSC-mid and aNSC-late cells
comb_fpkm<-cbind(group3_counts,group4_counts)

#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[match(c("Atp1a2","Ntsr2","Gja1","Jag1","Fgfr3","Dlx1","Dlx2"),rownames(comb_fpkm)),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_hiseq.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()

