
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
# greaterthan0<-allcounts_allcells_nooligo>10
# greaterthan0sum<-rowSums(greaterthan0)
# allcounts_allcells_nooligo_genefilt<-allcounts_allcells_nooligo[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_nooligo_noast_genefilt)
#glm <- DGEList(counts=allcounts_allcells_nooligo_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)


int_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])
#Genes to be used for correlation
comb_fpkm<-fpkm_glm_genefilt_nooligo
comb_fpkm_curr<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]

dlx2<-comb_fpkm_curr[match("Dlx2",rownames(comb_fpkm_curr)),]
dlx2_fac<-c(rep("#000000",length(dlx2)))
dlx2_fac[dlx2>0]<-"#9966FF"

Atp1a2<-comb_fpkm[match("Atp1a2",rownames(comb_fpkm)),]
Atp1a2_fac<-c(rep("#000000",length(Atp1a2)))
Atp1a2_fac[Atp1a2>0]<-"#9966FF"

Dlx6as1<-comb_fpkm[match("Dlx6as1",rownames(comb_fpkm)),]
Dlx6as1_fac<-c(rep("#000000",length(Dlx6as1)))
Dlx6as1_fac[Dlx6as1>0]<-"#9966FF"

Cdk1<-comb_fpkm[match("Cdk1",rownames(comb_fpkm)),]
Cdk1_fac<-c(rep("#000000",length(Cdk1)))
Cdk1_fac[Cdk1>0]<-"#9966FF"

Ccna2<-comb_fpkm[match("Ccna2",rownames(comb_fpkm)),]
Ccna2_fac<-c(rep("#000000",length(Ccna2)))
Ccna2_fac[Ccna2>0]<-"#9966FF"

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])

comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"#ff9933"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"#ff3333"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"#00CCFF"

comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#ff9933"

library(Hmisc)

#Calculate correlation
correlation<-rcorr(log(comb_fpkm_curr+1),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_consensusOrdering_cellbased.pdf",height=20,width=20)
heatmap.2(correlation$r,col=my_palette,density.info="none",trace="none", ColSideColors=dlx2_fac,keysize = 1)

dev.off()



comb_fpkm<-fpkm_glm_genefilt_nooligo

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])

comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"#ff9933"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"#ff3333"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"#00CCFF"

library(Hmisc)

#Calculate correlation
correlation<-rcorr(log(comb_fpkm+1),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
library(gplots)
pdf("Correlation_allgenes_cellbased.pdf",height=20,width=20)
heatmap.2(correlation$r,col=my_palette,density.info="none",trace="none", ColSideColors=dlx2_fac,keysize = 1)
dev.off()

pdf("Correlation_allgenes_cellbased_monoclecolors.pdf",height=20,width=20)
heatmap.2(correlation$r,col=my_palette,density.info="none",trace="none", ColSideColors=comb_fpkm_col,keysize = 1)
dev.off()

pdf("Correlation_allgenes_cellbased_celltypecolors.pdf",height=20,width=20)
heatmap.2(correlation$r,col=my_palette,density.info="none",trace="none", ColSideColors=comb_fpkm_col,keysize = 1)
dev.off()