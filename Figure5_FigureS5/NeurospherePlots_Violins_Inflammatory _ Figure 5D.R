#4/26/16 - Ben Dulken

#This script plots the expression of inflammatory genes that are associated with neurospheres.
#This was made a separate script to plot genes that were detected in neurospheres but were not detected
#in any in vivo cells.

#The code is essentially identical to that found in Neurosphere Plots _ PCA and Violin _ Figure 4B,C,D,E, Figure S4B

#Reading in NS and specific population cells
rm(list=ls())
library(ggplot2)
library(edgeR)

#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
allcounts_allcells<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]
#remove astrocytes for regression learning
allcounts_allcells_nooligo_noast<-allcounts_allcells_nooligo[,!grepl("Ast",colnames(allcounts_allcells_nooligo))]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
neurosphere<-read.table("AllCounts_Neurosphere_read_gene_filt_FINAL.txt")

allcounts_comb<-cbind(allcounts_allcells_nooligo_noast,neurosphere)

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_comb>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_comb_genefilt<-allcounts_comb[greaterthan0sum>=5,]


glm <- DGEList(counts=allcounts_comb_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)



setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])

group1_fpkm<-fpkm_glm_genefilt_nooligo[,match(group1,colnames(fpkm_glm_genefilt_nooligo))]
group2_fpkm<-fpkm_glm_genefilt_nooligo[,match(group2,colnames(fpkm_glm_genefilt_nooligo))]
group3_fpkm<-fpkm_glm_genefilt_nooligo[,match(group3,colnames(fpkm_glm_genefilt_nooligo))]
group4_fpkm<-fpkm_glm_genefilt_nooligo[,match(group4,colnames(fpkm_glm_genefilt_nooligo))]
group5_fpkm<-fpkm_glm_genefilt_nooligo[,match(group5,colnames(fpkm_glm_genefilt_nooligo))]
fpkm_NS<-fpkm_glm_genefilt_nooligo[,grepl("NS_",colnames(fpkm_glm_genefilt_nooligo))]
comb_fpkm<-cbind(group1_fpkm,group2_fpkm,group3_fpkm,group4_fpkm,group5_fpkm,fpkm_NS)



#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"#ff9933"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"#ff3333"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"#00CCFF"
comb_fpkm_col[grep("NS_",colnames(comb_fpkm))]<-"#000000"

comb_fac_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fac_col[match(group1,colnames(comb_fpkm))]<-"qNSC-like"
comb_fac_col[match(group2,colnames(comb_fpkm))]<-"aNSC-early"
comb_fac_col[match(group3,colnames(comb_fpkm))]<-"aNSC-mid"
comb_fac_col[match(group4,colnames(comb_fpkm))]<-"aNSC-late"
comb_fac_col[match(group5,colnames(comb_fpkm))]<-"NPC-like"
comb_fac_col[grep("NS_",colnames(comb_fpkm))]<-"NS"

comb_fac_2<-c(rep("#000000",length(comb_fpkm[1,])))
comb_fac_2[grepl("NS_",colnames(comb_fpkm))]<-"#ff4d4d"

#Larger plots used for Figure 4E, 4F
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
int_genes<-c("Ifitm3","Fas","Cx3cl1","Ccl2","Ifngr2","Rhog","Csf1")
for(i in 1:length(int_genes)){
  int_fpkm_glm<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]
  #comb_fpkm_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=log2(int_fpkm_glm[i,]+1),col=comb_fpkm_col,fac=comb_fac_col,colors=comb_fac_2)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+
    geom_jitter(aes(x=dataframe$fac,y=dataframe$data),position = position_jitter(width = .3),color=dataframe$colors,alpha=0.7)
  p<-p+theme_classic()
  p<-p + theme(legend.position="none")
  p<-p+scale_x_discrete(labels=c("qNSC-\nlike","aNSC-\nearly","aNSC-\nmid","aNSC-\nlate","NPC-\nlike","NS"))
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="log2(FPKM+1)")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+scale_fill_manual(values=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF","#000000"))
  
  pdf(paste(rownames(int_fpkm_glm)[i],"_expression_violin_withNS_large.pdf",sep=""),width=10,height=4)
  print(p)
  dev.off()
}