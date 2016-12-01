
#Differential Expression - calculation by SCDE
#Enrichment plots for GSEA enrichments (weighted, sorted by Z-score)
#Figure 3B-E


rm(list=ls())
library(edgeR)
library(scde)
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


#============================================================================================================================================#
#Establishing groups
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])


group1_counts<-allcounts_allcells_nooligo_noast_genefilt[,match(group1,colnames(allcounts_allcells_nooligo_noast_genefilt))]
group2_counts<-allcounts_allcells_nooligo_noast_genefilt[,match(group2,colnames(allcounts_allcells_nooligo_noast_genefilt))]
group3_counts<-allcounts_allcells_nooligo_noast_genefilt[,match(group3,colnames(allcounts_allcells_nooligo_noast_genefilt))]
group4_counts<-allcounts_allcells_nooligo_noast_genefilt[,match(group4,colnames(allcounts_allcells_nooligo_noast_genefilt))]
group5_counts<-allcounts_allcells_nooligo_noast_genefilt[,match(group5,colnames(allcounts_allcells_nooligo_noast_genefilt))]

comb_counts<-cbind(group1_counts,group2_counts,group3_counts,group4_counts,group5_counts)
comb_fac<-factor(c(rep("group1",length(group1_counts[1,])),
                   rep("group2",length(group2_counts[1,])),
                   rep("group3",length(group3_counts[1,])),
                   rep("group4",length(group4_counts[1,])),
                   rep("group5",length(group5_counts[1,]))))

#============================================================================================================================================#
#Setting up SCDE
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
n.cores <- 10;
# calculate models
o.ifm <- scde.error.models(counts=comb_counts,groups=comb_fac,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1);


valid.cells <- o.ifm$corr.a >0;
table(valid.cells)
o.ifm <- o.ifm[valid.cells,];

o.prior <- scde.expression.prior(models=o.ifm,counts=comb_counts,length.out=400,show.plot=F)



#Group1 v Group2

comb_fac<-factor(c(rep("group1",length(group1_counts[1,])),rep("group2",length(group2_counts[1,])),rep(NA,length(group3_counts[1,])),rep(NA,length(group4_counts[1,])),rep(NA,length(group5_counts[1,]))))

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group1_group2_ediff_full_10282015_round2_4groups_allgenesINIT_table50.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group1_group2_10282015_round2_4groups_allgenesINIT_table50_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")



#Group2 vs. Group3

comb_fac<-factor(c(rep(NA,length(group1_counts[1,])),rep("group2",length(group2_counts[1,])),rep("group3",length(group3_counts[1,])),rep(NA,length(group4_counts[1,])),rep(NA,length(group5_counts[1,]))))

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group2_group3_ediff_full_10282015_round2_4groups_allgenesINIT_table50.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group2_group3_10282015_round2_4groups_allgenesINIT_table50_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")




#Group3 vs. Group4

comb_fac<-factor(c(rep(NA,length(group1_counts[1,])),rep(NA,length(group2_counts[1,])),rep("group3",length(group3_counts[1,])),rep("group4",length(group4_counts[1,])),rep(NA,length(group5_counts[1,]))))

# define two groups of cells

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group3_group4_ediff_full_10282015_round2_4groups_allgenesINIT_table50.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group3_group4_10282015_round2_4groups_allgenesINIT_table50_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")



#Group4 vs. Group5

comb_fac<-factor(c(rep(NA,length(group1_counts[1,])),rep(NA,length(group2_counts[1,])),rep(NA,length(group3_counts[1,])),rep("group4",length(group4_counts[1,])),rep("group5",length(group5_counts[1,]))))

# define two groups of cells

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_group4_group5_ediff_full_10282015_round2_4groups_allgenesINIT_table50.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"ML_iterative_group4_group5_10282015_round2_4groups_allgenesINIT_table50_5groups.rnk",quote=F,row.names=F,col.names=F,sep="\t")




#============================================================================================================================================#
#============================================================================================================================================#
library(ggplot2)
#Group Enrichments - GSEA was given ranked lists of genes by z score for each transition.

#Figure 3B,C,D,E

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Enrichments/min")
enrichment<-read.table("group1_group2_GSEA_enrichments_11032015_min.txt")

#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
color<-vector()
enrich<-vector()
for(i in 1:length(enrichment[,2])){
  if(enrichment[i,2]<0){
    enrich_temp<-log10(abs(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#9966FF")
  }else{
    enrich_temp<-(-log10(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#c9cc00")
  }
}

data1<-data.frame(names=enrichment[length(enrichment[,1]):1,1],enrich=enrich[length(enrichment[,1]):1],color=color[length(color):1])
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
p<-ggplot(data1)
p<-p+geom_bar(aes(x=data1$names,y=data1$enrich),stat="identity",fill=data1$color)
p<-p+theme_classic()
p<-p+coord_flip()+labs(y="-log10(FDR)",x=NULL)
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=15,face="bold"))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_y_continuous(limits=c(-4,4),breaks=c(-4,-2,0,2,4), labels=c("4", "2", "0","2","4"))

p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
pdf("11032015_Group1_Group2_enrichments.pdf",height=3,width=8)
print(p)
dev.off()

#group2 vs group3
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Enrichments/min")
enrichment<-read.table("group2_group3_GSEA_enrichments_11032015_min.txt")

#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
color<-vector()
enrich<-vector()
for(i in 1:length(enrichment[,2])){
  if(enrichment[i,2]<0){
    enrich_temp<-log10(abs(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#c9cc00")
  }else{
    enrich_temp<-(-log10(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#ff9933")
  }
}

data1<-data.frame(names=enrichment[length(enrichment[,1]):1,1],enrich=enrich[length(enrichment[,1]):1],color=color[length(color):1])
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
p<-ggplot(data1)
p<-p+geom_bar(aes(x=data1$names,y=data1$enrich),stat="identity",fill=data1$color)
p<-p+theme_classic()
p<-p+coord_flip()+labs(y="-log10(FDR)",x=NULL)
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=15,face="bold"))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_y_continuous(limits=c(-4,4),breaks=c(-4,-2,0,2,4), labels=c("4", "2", "0","2","4"))

p

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
pdf("11032015_Group2_Group3_enrichments.pdf",height=3,width=8)
print(p)
dev.off()




#group3 vs group4
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Enrichments/min")
enrichment<-read.table("group3_group4_GSEA_enrichments_11032015_min.txt")

#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
color<-vector()
enrich<-vector()
for(i in 1:length(enrichment[,2])){
  if(enrichment[i,2]<0){
    enrich_temp<-log10(abs(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#ff9933")
  }else{
    enrich_temp<-(-log10(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#ff3333")
  }
}

data1<-data.frame(names=enrichment[length(enrichment[,1]):1,1],enrich=enrich[length(enrichment[,1]):1],color=color[length(color):1])
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
p<-ggplot(data1)
p<-p+geom_bar(aes(x=data1$names,y=data1$enrich),stat="identity",fill=data1$color)
p<-p+theme_classic()
p<-p+coord_flip()+labs(y="-log10(FDR)",x=NULL)
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=15,face="bold"))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_y_continuous(limits=c(-4,4),breaks=c(-4,-2,0,2,4), labels=c("4", "2", "0","2","4"))

p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
pdf("11032015_Group3_Group4_enrichments.pdf",height=3,width=8)
print(p)
dev.off()


#group4 vs group5
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Enrichments/min")
enrichment<-read.table("group4_group5_GSEA_enrichments_11032015_min.txt")
#c("#9966FF","#c9cc00","#ff9933","#ff3333","")
color<-vector()
enrich<-vector()
for(i in 1:length(enrichment[,2])){
  if(enrichment[i,2]<0){
    enrich_temp<-log10(abs(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#ff3333")
  }else{
    enrich_temp<-(-log10(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#00CCFF")
  }
}

data1<-data.frame(names=enrichment[length(enrichment[,1]):1,1],enrich=enrich[length(enrichment[,1]):1],color=color[length(color):1])
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
p<-ggplot(data1)
p<-p+geom_bar(aes(x=data1$names,y=data1$enrich),stat="identity",fill=data1$color)
p<-p+theme_classic()
p<-p+coord_flip()+labs(y="-log10(FDR)",x=NULL)
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=15,face="bold"))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+scale_y_continuous(limits=c(-4,4),breaks=c(-4,-2,0,2,4), labels=c("4", "2", "0","2","4"))

p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
pdf("11032015_Group4_Group5_enrichments.pdf",height=3,width=8)
print(p)
dev.off()
