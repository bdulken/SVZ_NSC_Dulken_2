#4/26/16 - Neurosphere and in vivo NSC differential expression and enrichment plotting

#This script performs differential expresion between all in vitro single NS cells and all
#in vivo NSCs

# The comparison that was used for the paper was the first one which is calculated which is 
#that between all in vitro NS and all aNSC and NPC cells (aNSC-early, aNSC-mid, aNSC-late, and NPC)



#Reading in NS and specific population cells
rm(list=ls())
library(ggplot2)
library(edgeR)
library(scde)

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
allcounts_allcells<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]
#remove astrocytes for regression learning
allcounts_allcells_nooligo_noast<-allcounts_allcells_nooligo[,!grepl("Ast",colnames(allcounts_allcells_nooligo))]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])

allcounts_group2<-allcounts_allcells_nooligo_noast[,match(group2,colnames(allcounts_allcells_nooligo_noast))]
allcounts_group3<-allcounts_allcells_nooligo_noast[,match(group3,colnames(allcounts_allcells_nooligo_noast))]
allcounts_group4<-allcounts_allcells_nooligo_noast[,match(group4,colnames(allcounts_allcells_nooligo_noast))]
allcounts_group5<-allcounts_allcells_nooligo_noast[,match(group5,colnames(allcounts_allcells_nooligo_noast))]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
neurosphere<-read.table("AllCounts_Neurosphere_read_gene_filt_FINAL.txt")

allcounts_comb<-cbind(allcounts_group2,allcounts_group3,allcounts_group4,allcounts_group5,neurosphere)

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_comb>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_comb_genefilt<-allcounts_comb[greaterthan0sum>=5,]


comb_counts<-allcounts_comb_genefilt
# comb_fac<-factor(c(rep("aNSC-mid",length(allcounts_group3[1,])),
#                    rep("aNSC-late",length(allcounts_group4[1,])),
#                    rep("NPC-like",length(allcounts_group5[1,])),
#                    rep("Neurosphere",length(neurosphere[1,]))))

comb_fac<-factor(c(rep("InVivo",length(allcounts_group2[1,])),
                   rep("InVivo",length(allcounts_group3[1,])),
                   rep("InVivo",length(allcounts_group4[1,])),
                   rep("InVivo",length(allcounts_group5[1,])),
                   rep("InVitro",length(neurosphere[1,]))))

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



#This was the comparison that was used for the enrichments.
#Comparing all aNSCs (aNSC-early,-mid, and -late) to the in vitro NS cells

comb_fac<-factor(c(rep("InVivo",length(allcounts_group2[1,])),
                   rep("InVivo",length(allcounts_group3[1,])),
                   rep("InVivo",length(allcounts_group4[1,])),
                   rep("InVivo",length(allcounts_group5[1,])),
                   rep("InVitro",length(neurosphere[1,]))))


# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_InVivo_G2345_vs_Neurosphere.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"InVivo_G2345_vs_Neurosphere.rnk",quote=F,row.names=F,col.names=F,sep="\t")



#only groups 3 and 4
comb_fac<-factor(c(rep("InVivo",length(allcounts_group2[1,])),
                   rep("InVivo",length(allcounts_group3[1,])),
                   rep("InVivo",length(allcounts_group4[1,])),
                   rep(NA,length(allcounts_group5[1,])),
                   rep("InVitro",length(neurosphere[1,]))))


# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_InVivo_G234_vs_Neurosphere.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"InVivo_G234_vs_Neurosphere.rnk",quote=F,row.names=F,col.names=F,sep="\t")



comb_fac<-factor(c(rep("InVivo",length(allcounts_group2[1,])),
                   rep("InVivo",length(allcounts_group3[1,])),
                   rep(NA,length(allcounts_group4[1,])),
                   rep(NA,length(allcounts_group5[1,])),
                   rep("InVitro",length(neurosphere[1,]))))


# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,comb_counts,o.prior,groups=comb_fac,n.randomizations=100,n.cores=n.cores,verbose=1)

head(ediff[order(ediff$Z,decreasing=F),])

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
sort_ediff<-ediff[order(ediff$Z,decreasing=F),]
write.table(sort_ediff,"sort_InVivo_G23_vs_Neurosphere.txt")
sort_ediff_rnk<-cbind(rownames(sort_ediff),sort_ediff[,5])
write.table(sort_ediff_rnk,"InVivo_G23_vs_Neurosphere.rnk",quote=F,row.names=F,col.names=F,sep="\t")















#Enrichments

#Figure 5C

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Enrichments")
enrichment<-read.table("inVitro_inVivo_enrichments.txt")

#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
color<-vector()
enrich<-vector()
for(i in 1:length(enrichment[,2])){
  if(enrichment[i,2]<0){
    enrich_temp<-log10(abs(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#ff704d")
  }else{
    enrich_temp<-(-log10(enrichment[i,2]))
    enrich<-c(enrich,enrich_temp)
    color<-c(color,"#666666")
  }
}

data1<-data.frame(names=enrichment[length(enrichment[,1]):1,1],enrich=enrich[length(enrichment[,1]):1],color=color[length(color):1])
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
p<-ggplot(data1)
p<-p+geom_bar(aes(x=data1$names,y=data1$enrich),stat="identity",fill=data1$color)
p<-p+theme_classic()
p<-p+coord_flip()+labs(y="log10(FDR)",x="pathways")
p<-p+theme(axis.text.x=element_text(size=15))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=15,face="bold"))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))

p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
pdf("Neurosphere_Enrichments.pdf",height=10,width=11)
print(p)
dev.off()
