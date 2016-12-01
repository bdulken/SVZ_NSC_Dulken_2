
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


#With cell cycle
fpkm_Ast<-fpkm_glm_genefilt_nooligo[,grepl("Ast",colnames(fpkm_glm_genefilt_nooligo))]
fpkm_qNSC<-fpkm_glm_genefilt_nooligo[,grepl("qNSC",colnames(fpkm_glm_genefilt_nooligo))]


comb_fpkm<-cbind(fpkm_Ast,fpkm_qNSC)

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
p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
#p<-p+geom_point(aes(x=data$x,y=data$y),color=data$col,size=7,alpha=0.65)
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
pdf("All_miseq_PCA_nooligos_Ast_qNSC_only_text.pdf",height=8,width=8)
print(p)
dev.off()



#With cell cycle
fpkm_Ast<-fpkm_glm_genefilt_nooligo[,grepl("Ast",colnames(fpkm_glm_genefilt_nooligo))]
fpkm_qNSC<-fpkm_glm_genefilt_nooligo[,grepl("qNSC",colnames(fpkm_glm_genefilt_nooligo))]

comb_fpkm<-cbind(fpkm_Ast,fpkm_qNSC)
comb_fpkm<-comb_fpkm[,-match(c("Ast_92_1g","qNSC_A_37_1g","qNSC_B_27_1g","qNSC_B_38_1g"),colnames(comb_fpkm))]
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
pdf("All_miseq_PCA_nooligos_Ast_qNSC_only_qNSC_outliersremoved.pdf",height=8,width=8)
print(p)
dev.off()



#Astrocyte, qNSC differential expression
#Remove weird cells (qNSCs that are actually activated, one strange Astrocyte)
allcounts_allcells_nooligo_genefilt_noqNSCOutlier<-allcounts_allcells_nooligo_genefilt[,-match(c("Ast_92_1g","qNSC_A_37_1g","qNSC_B_27_1g","qNSC_B_38_1g"),colnames(allcounts_allcells_nooligo_genefilt))]
allcounts_Ast<-allcounts_allcells_nooligo_genefilt_noqNSCOutlier[,grepl("Ast",colnames(allcounts_allcells_nooligo_genefilt_noqNSCOutlier))]
allcounts_qNSC<-allcounts_allcells_nooligo_genefilt_noqNSCOutlier[,grepl("qNSC",colnames(allcounts_allcells_nooligo_genefilt_noqNSCOutlier))]
comb_counts<-cbind(allcounts_Ast,allcounts_qNSC)
comb_fac<-factor(c(rep("Ast",length(allcounts_Ast[1,])),
                   rep("qNSC",length(allcounts_qNSC[1,]))))


n.cores <- 10;
# calculate models
o.ifm <- scde.error.models(counts=comb_counts,groups=comb_fac,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1);


valid.cells <- o.ifm$corr.a >0;
table(valid.cells)
o.ifm <- o.ifm[valid.cells,];

o.prior <- scde.expression.prior(models=o.ifm,counts=comb_counts,length.out=400,show.plot=F)

#Group1

comb_fac<-factor(c(rep("Ast",length(allcounts_Ast[1,])),
                   rep("qNSC",length(allcounts_qNSC[1,]))))

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_Ast_qNSC_diff.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"Ast_qNSC_diff.rnk",quote=F,row.names=F,col.names=F,sep="\t")