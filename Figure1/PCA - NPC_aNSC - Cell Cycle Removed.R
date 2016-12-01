#PCA with only aNSCs and NPCs, with cell cycle genes removed

rm(list=ls())
library(ggplot2)
library(edgeR)
library(ConsensusClusterPlus)
library(monocle)

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

#Removing Oligodendrocytes and Outliers
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_allcells_nooligo>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_nooligo_genefilt<-allcounts_allcells_nooligo[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_nooligo_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)

#Removing cell cycle genes
setwd("/Volumes/guacamole/scLMV/")
cellcycle<-as.vector(read.table("cellcyclesymbols.txt")[,1])
fpkm_glm_genefilt_nocycle<-fpkm_glm_genefilt_nooligo[-na.omit(match(cellcycle,rownames(fpkm_glm_genefilt_nooligo))),]

#=====================================================================================================================
#=====================================================================================================================
#=====================================================================================================================
#Just NPCs aNSCs
fpkm_aNSC<-fpkm_glm_genefilt_nocycle[,grepl("aNSC",colnames(fpkm_glm_genefilt_nocycle))]
fpkm_NPC<-fpkm_glm_genefilt_nocycle[,grepl("NPC",colnames(fpkm_glm_genefilt_nocycle))]

comb_fpkm<-cbind(fpkm_aNSC,fpkm_NPC)

#All Cells - No Oligos - All Genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"

#PCA
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = F, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-comb_fpkm_col
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm),col=col)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y),color=data$col,size=7,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#With cell cycle
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_aNSC_NPC_only_nocellcycle.pdf",height=8,width=8)
print(p)
dev.off()




#=====================================================================================================================
#=====================================================================================================================
#With cell cycle
fpkm_aNSC<-fpkm_glm_genefilt_nooligo[,grepl("aNSC",colnames(fpkm_glm_genefilt_nooligo))]
fpkm_NPC<-fpkm_glm_genefilt_nooligo[,grepl("NPC",colnames(fpkm_glm_genefilt_nooligo))]


comb_fpkm<-cbind(fpkm_aNSC,fpkm_NPC)

#All Cells - No Oligos - All Genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"

#PCA
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = F, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-comb_fpkm_col
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm),col=col)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y),color=data$col,size=7,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#FIGURE 1E
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_aNSC_NPC_only_withcellcycle.pdf",height=8,width=8)
print(p)
dev.off()





#=====================================================================================================================
#=====================================================================================================================
#All Cells
comb_fpkm<-fpkm_glm_genefilt_nocycle

#All Cells - No Oligos - All Genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"

#PCA
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = F, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-comb_fpkm_col
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm),col=col)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y),color=data$col,size=7,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#With cell cycle
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_allcells_nocellcycle.pdf",height=8,width=8)
print(p)
dev.off()




#=====================================================================================================================
#=====================================================================================================================
#With cell cycle
comb_fpkm<-fpkm_glm_genefilt_nooligo

#All Cells - No Oligos - All Genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC",colnames(comb_fpkm))]<-"#FF3300"
comb_fpkm_col[grepl("NPC",colnames(comb_fpkm))]<-"#00CCFF"

#PCA
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = F, center = T, retx=T)
PCA_results<-PCA_int$x
summa<-summary(PCA_int)
col<-comb_fpkm_col
data<-data.frame(x=PCA_results[,1],y=PCA_results[,2],factors=colnames(comb_fpkm),col=col)
data$factors<-as.character(data$factors)
data$factors <- factor(data$factors, levels=unique(data$factors), ordered = T)
p<-ggplot(data)
#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-p+geom_point(aes(x=data$x,y=data$y),color=data$col,size=7,alpha=0.65)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=24))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=24))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
#FIGURE 1E
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_allcells_withcellcycle.pdf",height=8,width=8)
print(p)
dev.off()