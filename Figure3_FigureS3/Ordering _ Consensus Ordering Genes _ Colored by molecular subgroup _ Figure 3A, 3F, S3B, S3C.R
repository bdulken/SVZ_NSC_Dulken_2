#This program generates pseudotime dependent gene expression plots for various genes colored by their putative molecular subgroup
# For Figures 3A, 3F, S3B, S3C


rm(list=ls())
library(edgeR)
library(monocle)
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


#=========================================================================================================================#
#Ordering and plotting genes in pseudotime with putative subpopulations

#For Figure 3A, Figure 3F, Figure S3B, Figure S3C

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])



comb_fpkm<-fpkm_glm_genefilt_nooligo
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"Group1"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"Group2"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"Group3"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"Group4"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"Group5"


#This pheno vector can be referenced by calling "color_by='type'" in any of the plotting 
#functions in monocle to color the cells plotted by their type.
pheno<-comb_fpkm_col

pheno.data.df <- data.frame(type=pheno)
rownames(pheno.data.df) <- colnames(comb_fpkm)

feature<-data.frame(rownames(comb_fpkm))


pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm, phenoData = pd)

#This is the list of ordering genes which is used to order the qNSCs, all have known roles in the differentiation of neurons,
#and because we have the qNSCs and Astrocytes in this as well I will also use some genes associated with quiescence

curated_ordering_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])

HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = T) # Order cells
if(as.character(pData(HSMM_q)[,1])[match(max(pData(HSMM_q)[,2]),pData(HSMM_q)[,2])]=="Group1"){
  HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = F) # Order cells
}

library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

#FIGURE S2C - Spanning Tree - Ordered by consensus ordering genes
pdf("spanningtree_group_separated.pdf",onefile=T)
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"))
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+theme(legend.position = "none")
p<-p+theme(axis.text.x=element_text(size=24))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(size=24))
print(p)
dev.off()

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
#Plotting genes with respect to pseudotime, colored by putative population (qNSC-like, aNSC-early, aNSC-mid, aNSC-late, NPC-like)
my_genes<-c("Atp1a2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Atp1a2_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Dlx6as1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx6as1_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Nrxn3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Nrxn3_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Sp8")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp8_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Sp9")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp9_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Gja1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Gja1_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Ntsr2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ntsr2_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Jag1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Jag1_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Dlx1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx1_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Dlx2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx2ONLY_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Lgr4")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Lgr4_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Fgfr3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Fgfr3_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Fgfr2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Fgfr2_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Slc1a3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Slc1a3_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Nrxn3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Nrxn3_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Dlx6as1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx6as1_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Dcx")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dcx_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Sp8")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp8_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
my_genes<-c("Sp9")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp9_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()





#Longer versions of plots used in Figure 3A
my_genes<-c("Id3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Id3_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Clu")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Clu_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Rpl32")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Rpl32_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Egfr")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Egfr_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Cdk1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Cdk1_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Ccna2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ccna2_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Dlx2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx2_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Dlx1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx1_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Dlx6as1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx6as1_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Dcx")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dcx_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()

my_genes<-c("Mki67")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Mki67_aggregated_ordering_top100_Figure2_groups_long.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(legend.position="none")
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
print(p)
dev.off()
