#This script contains the code for generating Figure 6C and 6D of the paper
#In which our cells and the Llorens Bobadilla cells are ordered using the consensus ordering genes.


#Ordering my own cells using the consensus ordering genes
rm(list=ls())
library(edgeR)
library(monocle)
source("/Volumes/guacamole/Software/R_Files_Packages_Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/R_Files_Packages_Functions/09282015_plot_genes_in_pseudotime_mod.txt")


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
#=========================================================================================================================#
#Ordering and plotting genes in pseudotime WITHOUT separating based on cell cycle or non cell cycle aNSC

#For Figure 6


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


pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm, phenoData = pd)



#This is the list of ordering genes which is used to order the qNSCs, all have known roles in the differentiation of neurons,
#and because we have the qNSCs and Astrocytes in this as well I will also use some genes associated with quiescence
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
curated_ordering_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])

match(curated_ordering_genes,rownames(comb_fpkm))
HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = T) # Order cells

pdf(paste("spanningtree_looped_ourcells_table50.pdf",sep=""))
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#e6e600","#0066ff"))
print(p)
dev.off()
library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

#Writing pseudotime for llorens-bobadilla data
write.table(pseudotime_mod,"In50models_pseudotime_mydat.txt")


#smaller versions of plots for Figure 6C and S6C - a subset of of the following plots are used
my_genes<-c("Id3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Id3_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Rpl32_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Cdk4_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Cdk1_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Dlx2_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Dcx_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Egfr_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Nrxn3_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Dlx6as1_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Ascl1_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Pax6_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Olig2_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Sp9_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
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
pdf("Sp8_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Clu")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Clu_aggregated_ordering_top100_Figure2_small.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#ff3333","#00CCFF","#9966FF"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


#------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------------#
rm(list=ls())
library(edgeR)
library(monocle)
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")

Llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

allcounts_allcells<-Llorens_allcounts

#Remove neuroblasts from Llorenss data
#allcounts_allcells_notaps<-allcounts_allcells[!grepl("tap",colnames(allcounts_allcells))]
allcounts_allcells_noblasts<-allcounts_allcells[!grepl("PSA",colnames(allcounts_allcells))]

#Remove oligodendrocytes from Llorenss data
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("comp_oligos.txt")[,1])
allcounts_allcells_noblasts_nooligo<-allcounts_allcells_noblasts[,-na.omit(match(oligos,colnames(allcounts_allcells_noblasts)))]
allcounts_allcells_noblasts_nooligo_noERCC<-allcounts_allcells_noblasts_nooligo[!grepl("ERCC-",rownames(allcounts_allcells_noblasts_nooligo)),]

#Genes detected in our dataset
detected_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/detected_genes_mydat.txt")[,1])

#Use only genes which were detected in our dataset
allcounts_allcells_noblasts_nooligo_noERCC_genefilt<-allcounts_allcells_noblasts_nooligo_noERCC[match(detected_genes,rownames(allcounts_allcells_noblasts_nooligo_noERCC)),]

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


#Ordering based on genes which appear in at least 50 models
comb_fpkm_comp<-fpkm_glm_genefilt_comp_cells
#Monocle ordering with Llorenss cells
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


#This is the list of ordering genes which is used to order the qNSCs, all have known roles in the differentiation of neurons,
#and because we have the qNSCs and Astrocytes in this as well I will also use some genes associated with quiescence
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
curated_ordering_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])

match(curated_ordering_genes,rownames(comb_fpkm_comp))
HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = T) # Order cells
if(as.character(pData(HSMM_q)[,1])[match(max(pData(HSMM_q)[,2]),pData(HSMM_q)[,2])]=="NSC"){
  HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = F) # Order cells
}

pdf(paste("spanningtree_looped_comp_table50.pdf",sep=""))
p<-plot_spanning_tree_mod(HSMM_q,color_by="type",color_custom=c("#e6e600","#0066ff"))
print(p)
dev.off()
library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

#Writing pseudotime for llorens-bobadilla data
write.table(pseudotime_mod,"pseudotime_comp_table50genes.txt")



#Plots for Figure 6D, Figure S6D
my_genes<-c("Id3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Id3_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Ascl1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ascl1_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Pax6")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Pax6_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Rpl32")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Rpl32_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Cdk1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Cdk1_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Egfr")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Egfr_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Dcx")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dcx_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Dlx2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx2_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Nrxn3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Nrxn3_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Sp8")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp8_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Dlx6as1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx6as1_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Sp8")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp8_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Sp9")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Sp9_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Bcl11b")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Bcl11b_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Atp1a2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Atp1a2_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Smo")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Smo_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Hes1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Hes1_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Slc1a3")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Slc1a3_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()



my_genes<-c("Gja1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Gja1_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Ntsr2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ntsr2_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Ascl1")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Ascl1_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()


my_genes<-c("Pax6")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Pax6_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Clu")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Clu_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()

my_genes<-c("Olig2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Olig2_aggregated_ordering_top100_mean_comp.pdf",height=1.3,width=4.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#e6e600","#0066ff"),trend_formula=NULL)
p<-p+theme(axis.line=element_blank())
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
