#Generation of PCAs with all single cells 
#Figure 1B

rm(list=ls())
library(ggplot2)
library(edgeR)


#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_allcells>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_genefilt<-allcounts_allcells[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt <- rpkm(glm.norm, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)


#Plotting PCAs all genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm<-fpkm_glm_genefilt

comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"

#PCA with variable genes
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = T, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-comb_fpkm_col
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm),col=col)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y),color=data$col,size=5,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#FIGURE 1B
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA.pdf",height=8,width=8)
print(p)
dev.off()

library(scatterplot3d)
pdf("All_miseq_PCA.pdf_3d.pdf",height=8,width=8)
scatterplot3d(PCA_results[,1:3],color=alpha(comb_fpkm_col,0.7),pch=16,cex.symbols=2.0,cex.axis=1.5,cex.lab=2.0,xlab=paste("PC1    (",round(summa$importance[2,1],digits = 3)*100,"% of variance)", sep = ""),
              ylab=paste("PC2    (",round(summa$importance[2,2],digits = 3)*100,"% of variance)", sep = ""),
              zlab=paste("PC3    (",round(summa$importance[2,3],digits = 3)*100,"% of variance)", sep = ""))
dev.off()


#Determining genes that define the oligodendrocyte distinction.
oligo_genes<-sort(PCA_int$rotation[,2][PCA_int$rotation[,2]<(-0.05)],decreasing=F)
write.table(oligo_genes,"Oligo_genes.txt")
write.table(rownames(comb_fpkm),"Background.txt")
