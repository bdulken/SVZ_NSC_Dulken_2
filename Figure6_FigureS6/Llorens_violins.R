rm(list=ls())
library(edgeR)
library(monocle)
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")

Llorens_allcounts<-read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Llorens_counts_allgenes.txt")

allcounts_allcells<-Llorens_allcounts

#Remove neuroblasts from Llorenss data
#allcounts_allcells_notaps<-allcounts_allcells[!grepl("tap",colnames(allcounts_allcells))]
#allcounts_allcells_noblasts<-allcounts_allcells[!grepl("PSA",colnames(allcounts_allcells))]

#Remove oligodendrocytes from Llorenss data
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("comp_oligos.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]
allcounts_allcells_nooligo_noERCC<-allcounts_allcells_nooligo[!grepl("ERCC-",rownames(allcounts_allcells_nooligo)),]

#Genes detected in our dataset
detected_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/detected_genes_mydat.txt")[,1])

#Use only genes which were detected in our dataset
allcounts_allcells_nooligo_noERCC_genefilt<-allcounts_allcells_nooligo_noERCC[match(detected_genes,rownames(allcounts_allcells_nooligo_noERCC)),]

glm <- DGEList(counts=allcounts_allcells_nooligo_noERCC_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_comp_cells <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                     normalized.lib.sizes=T, log=F)


#Ordering based on genes which appear in at least 50 models
comb_fpkm_comp_tap<-fpkm_glm_genefilt_comp_cells[,grepl("tap",colnames(fpkm_glm_genefilt_comp_cells))]
comb_fpkm_comp_nb<-fpkm_glm_genefilt_comp_cells[,grepl("PSA",colnames(fpkm_glm_genefilt_comp_cells))]
comb_fpkm_comp<-cbind(comb_fpkm_comp_tap,comb_fpkm_comp_nb)
#Monocle ordering with Llorenss cells
comb_fpkm_comp_col<-c(rep("NSC",length(comb_fpkm_comp[1,])))
comb_fpkm_comp_col[grepl("tap",colnames(comb_fpkm_comp))]<-"TAP\nLlorens-Bobadilla"
comb_fpkm_comp_col[grepl("PSA",colnames(comb_fpkm_comp))]<-"NB\nLlorens-Bobadilla"

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
comb_fpkm<-cbind(group5_fpkm)

#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"qNSC-like"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"aNSC-early"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"aNSC-mid"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"aNSC-late"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"NPC-like\nOur Study"


int_genes<-c("Cdk1","Cdk4","Dcx","Nrxn3","Dlx6as1","Dlx1","Dlx2","Prc1")
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure6_Output")
int_fpkm_glm_comp<-comb_fpkm_comp[match(int_genes,rownames(comb_fpkm_comp)),]
int_fpkm_glm<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]

int_fpkm_glm_combined<-cbind(int_fpkm_glm,int_fpkm_glm_comp)
comb_fpkm_combined_col<-c(comb_fpkm_col,comb_fpkm_comp_col)

for(i in 1:length(int_genes)){

  #comb_counts_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=log2(int_fpkm_glm_combined[i,]+1),col=comb_fpkm_combined_col,fac=comb_fpkm_combined_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  #p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data,color=dataframe$dot),position = position_jitter(width = .3),alpha=0.9)
  p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$fac),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data),position = position_jitter(width = .3,height=0),alpha=0.9)
  p<-p+ theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p + theme(legend.position="none")
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=15))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(plot.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="log2(FPKM+1)")
  p<-p+theme(axis.title.x=element_blank())
  p<-p+labs(title=int_genes[i])
  p<-p+scale_fill_manual(values=c("#00CCFF","#0066ff","#575757"))
  
  pdf(paste(rownames(int_fpkm_glm_combined)[i],"_neuroblast_comp.pdf",sep=""),width=7,height=3)
  print(p)
  dev.off()
}