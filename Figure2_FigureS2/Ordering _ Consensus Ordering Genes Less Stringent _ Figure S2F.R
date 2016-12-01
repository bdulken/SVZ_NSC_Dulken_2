#Identification of consensus ordering genes, ordering using consensus ordering genes


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




#Identification of consensus ordering genes
#Genes that appear in the top 100 genes in at least 50 of the models
library(reshape)
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Machine Learning/NewOutput_5/")
importance<-read.table("importance_matrix_original_all_based.txt")
importance_100<-as.matrix(importance[1:100,])
table_1<-table(importance_100)
table_50<-table_1[table_1>25]
#Writing consensus ordering genes
write.table(names(table_50),"/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/02012016_names_in25models_2.txt")









#Preparing objects for Monocle
comb_fpkm<-fpkm_glm_genefilt_nooligo
#comb_fpkm<-log(fpkm_glm_genefilt_nooligo+1)

comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"Ast"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"qNSC"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"aNSC"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"NPC"

#This pheno vector can be referenced by calling "color_by='type'" in any of the plotting 
#functions in monocle to color the cells plotted by their type.
pheno<-comb_fpkm_col

pheno.data.df <- data.frame(type=pheno)
rownames(pheno.data.df) <- colnames(comb_fpkm)

feature<-data.frame(rownames(comb_fpkm))

#Establishing Monocle object
pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm, phenoData = pd)



#This is the list of ordering genes which is used to order the qNSCs, all have known roles in the differentiation of neurons,
#and because we have the qNSCs and Astrocytes in this as well I will also use some genes associated with quiescence
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/LessStringentConsensus")

#Ordering genes = Consensus Ordering Genes
curated_ordering_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/02012016_names_in25models_2.txt")[,1])

match(curated_ordering_genes,rownames(comb_fpkm))
HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = T) # Order cells

#FIGURE S2C - Spanning Tree - Ordered by consensus ordering genes
pdf("spanningtree_looped_original_4groups_all_based_table50_lessstring.pdf",onefile=T)
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#ff3333","#00CCFF","#9966FF"))
p<-p+theme(axis.text.x=element_text(size=24))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(size=24))
print(p)
dev.off()

library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod


#Plotting genes for Figure 2D and S2D - A subset of the genes below are used - Cells are colored by FACS identity
#larger versions
my_genes<-c("Id3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Id3_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Rpl32")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Rpl32_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Cdk4")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Cdk4_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Cdk1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Cdk1_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Dlx2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Dcx")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dcx_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Egfr")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Egfr_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Ascl1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ascl1_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Pax6")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Pax6_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Olig2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Olig2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Notch2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Notch2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Nrxn3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Nrxn3_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Dlx6as1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx6as1_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Sp8")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp8_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Sp9")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp9_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Ccna2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ccna2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Ccnd2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ccnd2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Atp1a2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Atp1a2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Gja1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Gja1_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Tpx2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Tpx2_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Fgfr3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Fgfr3_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Lgr4")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Lgr4_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Jag1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Jag1_aggregated_ordering_top100_Figure2_lessstring.pdf",height=1.5,width=5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank())
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
#Establishing groups to be used in subsequent analyses
# group1<-rownames(pData(HSMM_q)[pData(HSMM_q)[,2]<=24,])
# group2<-rownames(pData(HSMM_q)[pData(HSMM_q)[,2]>=24&pData(HSMM_q)[,2]<47,])
# group3<-rownames(pData(HSMM_q)[pData(HSMM_q)[,2]>=47&pData(HSMM_q)[,2]<68,])
# group4<-rownames(pData(HSMM_q)[pData(HSMM_q)[,2]>=68&pData(HSMM_q)[,2]<85,])
# group5<-rownames(pData(HSMM_q)[pData(HSMM_q)[,2]>=85,])
# write.table(group1,"group1_10282015_50minModels_quiescent1_consensus.txt")
# write.table(group2,"group2_10282015_50minModels_active_nondividing_consensus.txt")
# write.table(group3,"group3_10282015_50minModels_active_dividing_early_consensus.txt")
# write.table(group4,"group4_10282015_50minModels_active_dividing_late_consensus.txt")
# write.table(group5,"group5_10282015_50minModels_active_dividing_consensus.txt")
# 