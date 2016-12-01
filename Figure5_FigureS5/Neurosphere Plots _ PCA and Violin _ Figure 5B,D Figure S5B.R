#============================================================================================================================================#
#============================================================================================================================================#
#PCA analysis with NS projection


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

#Filtering for expressed by 5 cells at 10 counts
greaterthan0<-allcounts_allcells_nooligo>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_nooligo_genefilt<-allcounts_allcells_nooligo[greaterthan0sum>=5,]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
neurosphere<-read.table("AllCounts_Neurosphere_read_gene_filt_FINAL.txt")
neurosphere_genefilt<-neurosphere[match(rownames(allcounts_allcells_nooligo_genefilt),rownames(neurosphere)),]

#Filtering genes for neurospheres based on detected genes in in vivo cells
allcounts_allcells_neurosphere_nooligo_genefilt<-cbind(allcounts_allcells_nooligo_genefilt,neurosphere_genefilt)

glm <- DGEList(counts=allcounts_allcells_neurosphere_nooligo_genefilt)
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
comb_fpkm<-cbind(group1_fpkm,group2_fpkm,group3_fpkm,group4_fpkm,group5_fpkm)

#============================================================================================================================================#
#============================================================================================================================================#
#Projections with NS cells using consensus ordering genes (Figure 5B)
#============================================================================================================================================#
#============================================================================================================================================#

var_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])

#Extract FPKM of consensus ordering genes to be used for PCA and projection
comb_fpkm<-comb_fpkm[na.omit(match(var_genes,rownames(comb_fpkm))),]
comb_fpkm_NS<-fpkm_NS[na.omit(match(var_genes,rownames(fpkm_NS))),]


#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"#ff9933"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"#ff3333"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"#00CCFF"

comb_fpkm_col<-c(comb_fpkm_col,rep("#000000",length(fpkm_NS[1,])))


#PCA
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = F, center = T, retx=T)
summa<-summary(PCA_int)
PCA_results<-PCA_int$x
predict_mydata<-predict(PCA_int,t(log2(comb_fpkm_NS+1)))

col<-comb_fpkm_col
data_2<-data.frame(x=c(PCA_results[,1],predict_mydata[,1]),y=c(PCA_results[,2],predict_mydata[,2]), col=col)

#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-ggplot()
p<-p+geom_point(data=data_2,aes(x=data_2$x,y=data_2$y),color=data_2$col,size=7,alpha=0.6)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2]-summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
#FIGURE 4B
pdf("All_miseq_PCA_projection_consensusgenes_allcells_withNS.pdf",height=8,width=8)
print(p)
dev.off()


#============================================================================================================================================#
#============================================================================================================================================#
#NS projection with all genes #currently not used in paper
#============================================================================================================================================#
#============================================================================================================================================#

#Perform PCA with all cells together
comb_fpkm<-cbind(group1_fpkm,group2_fpkm,group3_fpkm,group4_fpkm,group5_fpkm)
comb_fpkm_NS<-fpkm_NS
#============================================================================================================================================#
#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"#ff9933"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"#ff3333"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"#00CCFF"

comb_fpkm_col<-c(comb_fpkm_col,rep("#000000",length(fpkm_NS[1,])))


#PCA with variable genes
#PCA
PCA_int<-prcomp(t(log2(comb_fpkm+1)), scale = F, center = T, retx=T)
summa<-summary(PCA_int)
PCA_results<-PCA_int$x
predict_mydata<-predict(PCA_int,t(log2(comb_fpkm_NS+1)))

col<-comb_fpkm_col
data_2<-data.frame(x=c(PCA_results[,1],predict_mydata[,1]),y=c(PCA_results[,2],predict_mydata[,2]), col=col)

#p<-p+geom_text(aes(x=data$x,y=data$y),label=data$factors,color=data$col)
p<-ggplot()
p<-p+geom_point(data=data_2,aes(x=data_2$x,y=data_2$y),color=data_2$col,size=7,alpha=0.65)
p<-p+theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<- p+ labs(y = paste("PC2    (",(round(summa$importance[3,2], digits = 2)-round(summa$importance[3,1],digits = 2))*100,"% of variance)", sep = ""), x =paste("PC1    (",round(summa$importance[3,1],digits = 2)*100,"% of variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=22))
p<-p+theme(axis.title.x=element_text(size=26))
p<-p+theme(axis.text.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(size=26))
p<-p+theme(plot.title=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
#Currently not used with paper
pdf("All_miseq_PCA_Projection_allgenes_allcells_withNS.pdf",height=8,width=8)
print(p)
dev.off()
#============================================================================================================================================#
#============================================================================================================================================#


#Violin Plots with NS cells
#FIGURE 5D   FIGURE S5B
#============================================================================================================================================#
#============================================================================================================================================#



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

#Larger plots used for Figure 5D
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
int_genes<-c("Dlx2","Dlx1","Jag1","Notch2","Atp1a2","Ntsr2","Gja1","Ascl1","Cdk1","Notch2","Notch1","Dcx","Tubb3","Lgr4","Serpine2","Fgfr3","Fgfr2","Dlx6as1","Egfr","Nrxn3")
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
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
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

#Smaller plots used for Figure S5B
int_genes<-c("Dlx2","Dlx1","Jag1","Notch2","Atp1a2","Ntsr2","Gja1","Ascl1","Cdk1","Notch2","Notch1","Dcx","Tubb3","Lgr4","Serpine2","Fgfr3","Fgfr2","Dlx6as1","Ifitm3","Nrxn3","Cdk4","Mki67","Ccnb1","Ccna2")
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
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
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
  
  pdf(paste(rownames(int_fpkm_glm)[i],"_expression_violin_withNS_small.pdf",sep=""),width=10,height=3)
  print(p)
  dev.off()
}

#============================================================================================================================================#
#============================================================================================================================================#
#============================================================================================================================================#
#============================================================================================================================================#


