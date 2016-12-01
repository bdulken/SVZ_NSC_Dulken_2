#This script is used to generate the figures for Figure 6A and 6B

#------------------------------------------------------------------------------------------------------------------------------#
#In this section I filter genes based on my dataset and normalized the cell sets separately
#------------------------------------------------------------------------------------------------------------------------------#

rm(list=ls())
library(ggplot2)
library(edgeR)

Llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

#Remove neuroblasts from Llorenss data
allcounts_allcells<-Llorens_allcounts
#allcounts_allcells_notaps<-allcounts_allcells[!grepl("tap",colnames(allcounts_allcells))]
allcounts_allcells_noblasts<-allcounts_allcells[!grepl("PSA",colnames(allcounts_allcells))]

#Remove oligos from Llorenss data
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("comp_oligos.txt")[,1])
allcounts_allcells_noblasts_nooligo<-allcounts_allcells_noblasts[,-na.omit(match(oligos,colnames(allcounts_allcells_noblasts)))]
allcounts_allcells_noblasts_nooligo_noERCC<-allcounts_allcells_noblasts_nooligo[!grepl("ERCC-",rownames(allcounts_allcells_noblasts_nooligo)),]

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_allcells_nooligo>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_nooligo_genefilt<-allcounts_allcells_nooligo[greaterthan0sum>=5,]


#Consider genes which are detected in our dataset
filter_genes<-rownames(allcounts_allcells_nooligo_genefilt)

allcounts_allcells_noblasts_nooligo_noERCC_genefilt<-allcounts_allcells_noblasts_nooligo_noERCC[match(filter_genes,rownames(allcounts_allcells_noblasts_nooligo_noERCC)),]

glm <- DGEList(counts=allcounts_allcells_nooligo_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo_mycells <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                          normalized.lib.sizes=T, log=F)

glm <- DGEList(counts=allcounts_allcells_noblasts_nooligo_noERCC_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo_comp_cells <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                             normalized.lib.sizes=T, log=F)




#------------------------------------------------------------------------------------------------------------------------------#
#PCA using cur_genes (either variable genes or curated genes (machine learning)) - cells normalized by experiment, genes filtered by my dataset
#Pca other data, projection with my data
#------------------------------------------------------------------------------------------------------------------------------#

#Figure 6B - Llorens-Bobadilla Projected - 2500 variable genes
cur_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/var.genes_cutoff_0.5_2500genes.txt")[,1])
comb_fpkm<-fpkm_glm_genefilt_nooligo_mycells[na.omit(match(cur_genes,rownames(fpkm_glm_genefilt_nooligo_mycells))),]


comb_fpkm_comp<-fpkm_glm_genefilt_nooligo_comp_cells[na.omit(match(cur_genes,rownames(fpkm_glm_genefilt_nooligo_comp_cells))),]
#comb_fpkm_col<-c(rep("#66FF33",length(comb_fpkm[1,])))
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"


comb_fpkm_comp_col<-c(rep("#e6e600",length(comb_fpkm_comp[1,])))
comb_fpkm_comp_col[grepl("tap",colnames(comb_fpkm_comp))]<-"#0066ff"

PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = T, center = T, retx=T)
summa<-summary(PCA_int)
PCA_results<-PCA_int$x
predict_comp<-predict(PCA_int,t(log2(comb_fpkm_comp+1)))

col<-c(comb_fpkm_col,comb_fpkm_comp_col)
data_2<-data.frame(x=c(PCA_results[,1],predict_comp[,1]),y=c(PCA_results[,2],predict_comp[,2]), col=col)

#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-ggplot()
p<-p+geom_point(data=data_2,aes(x=data_2$x,y=data_2$y),color=data_2$col,size=7,alpha=0.65)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<- p+ labs(y = paste("PC2    (",(round((summa$importance[3,2]-summa$importance[3,1]),digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
pdf("All_miseq_Llorens_PCA_Llorens_projectionMyCells_2500vargenes.pdf",height=7,width=7)
print(p)
dev.off()

pred<-cbind(c(PCA_results[,1],predict_comp[,1]),c(PCA_results[,2],predict_comp[,2]),c(PCA_results[,3],predict_comp[,3]))
library(scatterplot3d)
pdf("All_miseq_Llorens_PCA_Llorens_projectionMyCells_2500vargenes_3d.pdf",height=8,width=8)
scatterplot3d(pred,color=alpha(col,0.7),pch=16,cex.symbols=2.0,cex.axis=1.5,cex.lab=2.0,xlab=paste("PC1    (",round(summa$importance[2,1],digits = 3)*100,"% of variance)", sep = ""),
              ylab=paste("PC2    (",round(summa$importance[2,2],digits = 3)*100,"% of variance)", sep = ""),
              zlab=paste("PC3    (",round(summa$importance[2,3],digits = 3)*100,"% of variance)", sep = ""),
              angle=150)
dev.off()

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#Figures S6A - No projection - 2500 variable genes
comb_fpkm<-cbind(comb_fpkm,comb_fpkm_comp)

PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = T, center = T, retx=T)
summa<-summary(PCA_int)
PCA_results<-PCA_int$x

col<-c(comb_fpkm_col,comb_fpkm_comp_col)
data_2<-data.frame(x=c(PCA_results[,1]),y=c(PCA_results[,2]), col=col)

#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-ggplot()
p<-p+geom_point(data=data_2,aes(x=data_2$x,y=data_2$y),color=data_2$col,size=7,alpha=0.65)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<- p+ labs(y = paste("PC2    (",(round((summa$importance[3,2]-summa$importance[3,1]),digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
pdf("All_miseq_Llorens_PCA_Llorens_noProjection_2500vargenes.pdf",height=7,width=7)
print(p)
dev.off()

library(scatterplot3d)
pdf("All_miseq_Llorens_PCA_Llorens_noProjection_2500vargenes_3d.pdf",height=8,width=8)
scatterplot3d(PCA_results[,1:3],color=alpha(col,0.7),pch=16,cex.symbols=2.0,cex.axis=1.5,cex.lab=2.0,xlab=paste("PC1    (",round(summa$importance[2,1],digits = 3)*100,"% of variance)", sep = ""),
              ylab=paste("PC2    (",round(summa$importance[2,2],digits = 3)*100,"% of variance)", sep = ""),
              zlab=paste("PC3    (",round(summa$importance[2,3],digits = 3)*100,"% of variance)", sep = ""))
dev.off()


#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#Figure S6B - Llorens-Bobadillas Projected - Consensus Ordering Genes
cur_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])
comb_fpkm<-fpkm_glm_genefilt_nooligo_mycells[na.omit(match(cur_genes,rownames(fpkm_glm_genefilt_nooligo_mycells))),]


comb_fpkm_comp<-fpkm_glm_genefilt_nooligo_comp_cells[na.omit(match(cur_genes,rownames(fpkm_glm_genefilt_nooligo_comp_cells))),]
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"


comb_fpkm_comp_col<-c(rep("#e6e600",length(comb_fpkm_comp[1,])))
comb_fpkm_comp_col[grepl("tap",colnames(comb_fpkm_comp))]<-"#0066ff"

PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = T, center = T, retx=T)
summa<-summary(PCA_int)
PCA_results<-PCA_int$x
predict_comp<-predict(PCA_int,t(log2(comb_fpkm_comp+1)))

col<-c(comb_fpkm_col,comb_fpkm_comp_col)
data_2<-data.frame(x=c(PCA_results[,1],predict_comp[,1]),y=c(PCA_results[,2],predict_comp[,2]), col=col)

#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-ggplot()
p<-p+geom_point(data=data_2,aes(x=data_2$x,y=data_2$y),color=data_2$col,size=7,alpha=0.65)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<- p+ labs(y = paste("PC2    (",(round((summa$importance[3,2]-summa$importance[3,1]),digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
pdf("All_miseq_Llorens_PCA_Llorens_projectionMyCells_conesnsusgenes.pdf",height=7,width=7)
print(p)
dev.off()

pred<-cbind(c(PCA_results[,1],predict_comp[,1]),c(PCA_results[,2],predict_comp[,2]),c(PCA_results[,3],predict_comp[,3]))
library(scatterplot3d)
pdf("All_miseq_Llorens_PCA_Llorens_projectionMyCells_conesnsusgenes_3d.pdf",height=8,width=8)
scatterplot3d(pred,color=alpha(col,0.7),pch=16,cex.symbols=2.0,cex.axis=1.5,cex.lab=2.0,xlab=paste("PC1    (",round(summa$importance[2,1],digits = 3)*100,"% of variance)", sep = ""),
              ylab=paste("PC2    (",round(summa$importance[2,2],digits = 3)*100,"% of variance)", sep = ""),
              zlab=paste("PC3    (",round(summa$importance[2,3],digits = 3)*100,"% of variance)", sep = ""),
              angle=150)
dev.off()