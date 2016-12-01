rm(list=ls())
library(edgeR)
library(monocle)
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")

# ====================================================================================================
#Repeat analysis with HiSeq cells only.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/")
run1<-read.table("AllCounts_HiSeq_run1.txt")
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2/")
run2<-read.table("AllCounts_HiSeq_run2.txt")

run1_NPC<-run1[,grepl("NPC",colnames(run1))]
run1_noNPC<-run1[,!grepl("NPC",colnames(run1))]

comb<-run1_noNPC+run2

comb_all<-cbind(comb,run1_NPC)


#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-spec_pops

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells[,-na.omit(match(oligos,colnames(allcounts_allcells)))]

comb_all_cellfilt<-comb_all[,!is.na(match(colnames(comb_all),colnames(allcounts_allcells_nooligo)))]

comb_all_cellfilt_noast<-comb_all_cellfilt[,!grepl("Ast",colnames(comb_all_cellfilt))]
#Filtering for expressed by 5 cells at 10 counts

greaterthan0<-comb_all_cellfilt_noast>10
greaterthan0sum<-rowSums(greaterthan0)
comb_all_cellfilt_noast_genefilt<-comb_all_cellfilt_noast[greaterthan0sum>=5,]


glm <- DGEList(counts=comb_all_cellfilt_noast_genefilt)
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


group1_counts<-fpkm_glm_genefilt_nooligo[,na.omit(match(group1,colnames(fpkm_glm_genefilt_nooligo)))]
group2_counts<-fpkm_glm_genefilt_nooligo[,na.omit(match(group2,colnames(fpkm_glm_genefilt_nooligo)))]
group3_counts<-fpkm_glm_genefilt_nooligo[,na.omit(match(group3,colnames(fpkm_glm_genefilt_nooligo)))]
group4_counts<-fpkm_glm_genefilt_nooligo[,na.omit(match(group4,colnames(fpkm_glm_genefilt_nooligo)))]
group5_counts<-fpkm_glm_genefilt_nooligo[,na.omit(match(group5,colnames(fpkm_glm_genefilt_nooligo)))]


comb_fpkm<-cbind(group1_counts,group2_counts,group3_counts,group4_counts,group5_counts)


#c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF")
comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"#9966FF"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"#c9cc00"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"#ff9933"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"#ff3333"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"#00CCFF"


comb_fac_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fac_col[match(group1,colnames(comb_fpkm))]<-"qNSC-like"
comb_fac_col[match(group2,colnames(comb_fpkm))]<-"aNSC-early"
comb_fac_col[match(group3,colnames(comb_fpkm))]<-"aNSC-mid"
comb_fac_col[match(group4,colnames(comb_fpkm))]<-"aNSC-late"
comb_fac_col[match(group5,colnames(comb_fpkm))]<-"NPC-like"

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure3_Output")
int_genes<-c("Dlx2","Dlx1","Atp1a2","Ntsr2","Gja1","Ascl1","Cdk1","Fgfr3","Fgfr2","Dlx6as1","Egfr","Prc1","Tpx2","Ccna2","Ccnb1","Nono","Mki67")
for(i in 1:length(int_genes)){
  int_fpkm_glm<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]
  #comb_fpkm_col<-c("#009900","#3366FF","#990000","#CC6600","#FF0066","#000000")
  dataframe<-data.frame(data=log2(int_fpkm_glm[i,]+1),col=comb_fpkm_col,fac=comb_fac_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe[,1],fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+
    geom_jitter(aes(x=dataframe$fac,y=dataframe[,1]),position = position_jitter(width = .3),alpha=0.7)
  p<-p+theme_classic()
  p<-p + theme(legend.position="none")
  p<-p+scale_x_discrete(labels=c("qNSC-\nlike","aNSC-\nearly","aNSC-\nmid","aNSC-\nlate","NPC-\nlike"))
  p<-p+theme(axis.text.y=element_text(size=15))
  p<-p+theme(axis.text.x=element_text(size=20))
  p<-p+theme(axis.title=element_text(size=20))
  p<-p+theme(axis.title.y=element_text(vjust=1))
  p<-p+theme(axis.title.x=element_text(vjust=-0.10))
  p<-p+labs(y="log2(FPKM+1)")
  p<-p+theme(axis.title.x=element_blank())
  
  p<-p+scale_fill_manual(values=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"))
  
  pdf(paste(rownames(int_fpkm_glm)[i],"_expression_violin_hiseq_group3_group4_large.pdf",sep=""),width=10,height=4)
  print(p)
  dev.off()
}