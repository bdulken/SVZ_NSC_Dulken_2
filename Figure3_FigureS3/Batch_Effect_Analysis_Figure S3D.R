#============================================================================================================================================#
#============================================================================================================================================#
#Batch effect analysis
#FIGURE S3D
rm(list=ls())
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")

group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus.txt")[,1])
group2_aNSC<-group2[grepl("aNSC",group2)]
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus.txt")[,1])
group3_aNSC<-group3[grepl("aNSC",group3)]
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus.txt")[,1])
group4_aNSC<-group4[grepl("aNSC",group4)]


group2_aNSC_percents<-c((length(group2_aNSC[grepl("aNSC_A",group2_aNSC)])/length(group2_aNSC)),(length(group2_aNSC[grepl("aNSC_B",group2_aNSC)])/length(group2_aNSC)),(length(group2_aNSC[grepl("aNSC_C",group2_aNSC)])/length(group2_aNSC)),(length(group2_aNSC[grepl("aNSC_D",group2_aNSC)])/length(group2_aNSC)))

group3_aNSC_percents<-c((length(group3_aNSC[grepl("aNSC_A",group3_aNSC)])/length(group3_aNSC)),(length(group3_aNSC[grepl("aNSC_B",group3_aNSC)])/length(group3_aNSC)),(length(group3_aNSC[grepl("aNSC_C",group3_aNSC)])/length(group3_aNSC)),(length(group3_aNSC[grepl("aNSC_D",group3_aNSC)])/length(group3_aNSC)))

group4_aNSC_percents<-c((length(group4_aNSC[grepl("aNSC_A",group4_aNSC)])/length(group4_aNSC)),(length(group4_aNSC[grepl("aNSC_B",group4_aNSC)])/length(group4_aNSC)),(length(group4_aNSC[grepl("aNSC_C",group4_aNSC)])/length(group4_aNSC)),(length(group4_aNSC[grepl("aNSC_D",group4_aNSC)])/length(group4_aNSC)))

groups<-rbind(group2_aNSC_percents,
              group3_aNSC_percents,
              group4_aNSC_percents)

names<-c("aNSC_early","aNSC_mid","aNSC_late")

groups_data<-data.frame(names,groups)

library(reshape)
groups_melt<-melt(groups_data,id.var="names")


p<-ggplot(groups_melt)

p<-p+geom_bar(aes(x=factor(groups_melt$names,levels=names<-c("aNSC_early","aNSC_mid","aNSC_late"),ordered=T),
                  y=as.numeric(groups_melt$value),fill=groups_melt$variable),stat="identity",width=0.6)
p<-p+labs(x="Cell grouping by monocle ordering",y="Cumulative fraction of aNSC batch")
p<-p+theme_classic()
p<-p+theme(legend.position="none")
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+scale_fill_manual(values=c("#663300","#cc0099","#cc9900","#cc0000"))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
pdf("GroupComposition_5groups_50minModels_11052015_aNSC_batch.pdf",height=9,width=11)
print(p)
dev.off()



#============================================================================================================================================#
#============================================================================================================================================#
#PCA batch effect analysis 


#None of this is used for the paper currently DO NOT CHECK BEYOND HERE

#Figure 1E initial PCA of NSC only cells - with oligodendrocytes removed
rm(list=ls())
library(ggplot2)
library(edgeR)


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






#All Cells - No Oligos - All Genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm<-fpkm_glm_genefilt_nooligo
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC_A",colnames(comb_fpkm))]<-"#663300"
comb_fpkm_col[grepl("aNSC_B",colnames(comb_fpkm))]<-"#cc0099"
comb_fpkm_col[grepl("aNSC_C",colnames(comb_fpkm))]<-"#cc9900"
comb_fpkm_col[grepl("aNSC_D",colnames(comb_fpkm))]<-"#cc0000"
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

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_aNSC_Batch.pdf",height=8,width=8)
print(p)
dev.off()


comb_fpkm<-fpkm_glm_genefilt_nooligo[,grepl("aNSC",colnames(fpkm_glm_genefilt_nooligo))]
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC",colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[grepl("aNSC_A",colnames(comb_fpkm))]<-"#663300"
comb_fpkm_col[grepl("aNSC_B",colnames(comb_fpkm))]<-"#cc0099"
comb_fpkm_col[grepl("aNSC_C",colnames(comb_fpkm))]<-"#cc9900"
comb_fpkm_col[grepl("aNSC_D",colnames(comb_fpkm))]<-"#cc0000"
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

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_aNSC_Batch_aNSConly.pdf",height=8,width=8)
print(p)
dev.off()











#All Cells - No Oligos - All Genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
comb_fpkm<-fpkm_glm_genefilt_nooligo
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")

comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC_A",colnames(comb_fpkm))]<-"#663300"
comb_fpkm_col[grepl("qNSC_B",colnames(comb_fpkm))]<-"#cc0099"
comb_fpkm_col[grepl("qNSC_C",colnames(comb_fpkm))]<-"#cc9900"
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

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_qNSC_Batch.pdf",height=8,width=8)
print(p)
dev.off()


comb_fpkm<-fpkm_glm_genefilt_nooligo[,grepl("qNSC",colnames(fpkm_glm_genefilt_nooligo))]
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[grepl("Ast",colnames(comb_fpkm))]<-"#009900"
comb_fpkm_col[grepl("qNSC_A",colnames(comb_fpkm))]<-"#663300"
comb_fpkm_col[grepl("qNSC_B",colnames(comb_fpkm))]<-"#cc0099"
comb_fpkm_col[grepl("qNSC_C",colnames(comb_fpkm))]<-"#cc9900"
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

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure1_Output")
pdf("All_miseq_PCA_nooligos_aNSC_Batch_qNSConly.pdf",height=8,width=8)
print(p)
dev.off()


