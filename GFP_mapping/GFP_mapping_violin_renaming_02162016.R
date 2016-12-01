#GFP mapping violin plots





#Reading in GFP mapped reads

#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_GFP/")
fileList=list.files(pattern=".counts.txt$")

#Combine the counts from each of the individual files to make one large counts matrix and save the counts matrix to
#a file.
allcounts<-read.table(fileList[1],header=F)
colnames(allcounts)[2]<-strsplit(fileList[1],split="[.]")[[1]][1]

for (i in 2:length(fileList)){
  temp = read.table(fileList[i],header=F)
  allcounts=cbind(allcounts,temp[,2])
  colnames(allcounts)[i+1]<-strsplit(fileList[i],split="[.]")[[1]][1]
}

allcounts<-allcounts[-c((length(allcounts[,1])-4):length(allcounts[,1])),] #removes HTseq QC values from bottom of counts matrix
allcounts_fin<-allcounts[,2:length(allcounts[1,])]
rownames(allcounts_fin)<-allcounts[,1]

allcounts_allcells<-cbind(allcounts_fin)

#The following code changes the names to match the naming approach which I used to name the cells in previous analyses.
str_vec<-vector()
for(i in 1:length(colnames(allcounts_allcells))){
  str<-""
  if(strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))[[1]][1]=="PGE120"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))
    num<-strsplit(filename[[1]][2],split=c("[-]"))[[1]][2]
    qual<-strsplit(filename[[1]][2],split=c("[-]"))[[1]][3]
    str<-paste(str,"aNSC_C_",num,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))[[1]][1]=="PGE210"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))
    num<-filename[[1]][3]
    qual<-filename[[1]][4]
    str<-paste(str,"aNSC_D_",num,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))[[1]][1]=="PG210"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))
    num<-filename[[1]][3]
    qual<-filename[[1]][4]
    str<-paste(str,"qNSC_C_",num,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="Na"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    str<-paste(str,"aNSC_A_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="Nq"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"qNSC_A_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="PGE"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"aNSC_B_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="PG"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"qNSC_B_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="G"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"Ast_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="E"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"NPC_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="N"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"NS_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }
}

check<-cbind(str_vec,colnames(allcounts_allcells))

allcounts_allcells_names<-allcounts_allcells
colnames(allcounts_allcells_names)<-str_vec
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
write.table(allcounts_allcells_names,"Allcounts_GFP_mapped.txt")


#The following code changes the names to match the naming approach which I used to name the cells in previous analyses. With the correct quality value
str_vec<-vector()
for(i in 1:length(colnames(allcounts_allcells))){
  str<-""
  if(strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))[[1]][1]=="PGE120"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))
    num<-strsplit(filename[[1]][2],split=c("[-]"))[[1]][2]
    qual<-strsplit(filename[[1]][2],split=c("[-]"))[[1]][3]
    str<-paste(str,"aNSC_C_",num,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))[[1]][1]=="PGE210"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))
    num<-filename[[1]][3]
    qual<-filename[[1]][4]
    str<-paste(str,"aNSC_D_",num,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))[[1]][1]=="PG210"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[_]"))
    num<-filename[[1]][3]
    qual<-filename[[1]][4]
    str<-paste(str,"qNSC_C_",num,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="Na"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    str<-paste(str,"aNSC_A_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="Nq"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"qNSC_A_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="PGE"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"aNSC_B_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="PG"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"qNSC_B_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="G"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"Ast_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="E"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"NPC_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][2]=="N"){
    filename<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))
    num<-gsub("C","",strsplit(filename[[1]][3],split=c("[_]"))[[1]][1])
    qual<-strsplit(colnames(allcounts_allcells)[i],split=c("[-]"))[[1]][1]
    str<-paste(str,"NS_",num,sep="")
    str<-paste(str,"_",qual,sep="")
    str_vec<-c(str_vec,str)
  }
}

allcounts_allcells_names<-allcounts_allcells
colnames(allcounts_allcells_names)<-str_vec
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
write.table(allcounts_allcells_names,"Allcounts_GFP_mapped_qualfromfile.txt")









#Reading in NS and specific population cells
rm(list=ls())
library(ggplot2)
library(edgeR)
library(ConsensusClusterPlus)
library(monocle)

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

gfp_mapped<-read.table("Allcounts_GFP_mapped.txt")

gfp_mapped_nooligo<-gfp_mapped[,match(colnames(allcounts_allcells_nooligo),colnames(gfp_mapped))]
gfp_mapped_nooligo_noERCC<-gfp_mapped_nooligo[!grepl("ERCC",rownames(gfp_mapped_nooligo)),]


greaterthan0<-gfp_mapped_nooligo_noERCC>10
greaterthan0sum<-rowSums(greaterthan0)
gfp_mapped_nooligo_noERCC_genefilt<-gfp_mapped_nooligo_noERCC[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_nooligo_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
exonic.gene.sizes.2<-c(exonic.gene.sizes,720)
names(exonic.gene.sizes.2)<-c(names(exonic.gene.sizes),'gfp')
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes.2))
exonic.gene.sizes.ord <- exonic.gene.sizes.2[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                      normalized.lib.sizes=T, log=F)

glm <- DGEList(counts=gfp_mapped_nooligo_noERCC_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
exonic.gene.sizes.2<-c(exonic.gene.sizes,720)
names(exonic.gene.sizes.2)<-c(names(exonic.gene.sizes),'gfp')
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes.2))
exonic.gene.sizes.ord <- exonic.gene.sizes.2[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt_nooligo_gfp <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])

group1_fpkm<-fpkm_glm_genefilt_nooligo_gfp[,match(group1,colnames(fpkm_glm_genefilt_nooligo_gfp))]
group2_fpkm<-fpkm_glm_genefilt_nooligo_gfp[,match(group2,colnames(fpkm_glm_genefilt_nooligo_gfp))]
group3_fpkm<-fpkm_glm_genefilt_nooligo_gfp[,match(group3,colnames(fpkm_glm_genefilt_nooligo_gfp))]
group4_fpkm<-fpkm_glm_genefilt_nooligo_gfp[,match(group4,colnames(fpkm_glm_genefilt_nooligo_gfp))]
group5_fpkm<-fpkm_glm_genefilt_nooligo_gfp[,match(group5,colnames(fpkm_glm_genefilt_nooligo_gfp))]
comb_fpkm<-cbind(group1_fpkm,group2_fpkm,group3_fpkm,group4_fpkm,group5_fpkm)



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

#Larger plots used for Figure 4E, 4F
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/GFP_mapping")
int_genes<-c("Dlx2","Dlx1","Jag1","Notch2","Atp1a2","Ntsr2","Gja1","Ascl1","Cdk1","Notch2","Notch1","Dcx","Tubb3","Lgr4","Serpine2","Fgfr3","Fgfr2","Dlx6as1","Egfr","GFP","Ccna2","Ccnb1")
for(i in 1:length(int_genes)){
  int_fpkm_glm<-comb_fpkm[match(int_genes,rownames(comb_fpkm)),]
  dataframe<-data.frame(data=log2(int_fpkm_glm[i,]+1),col=comb_fpkm_col,fac=comb_fac_col)
  dataframe$col<-as.character(dataframe$col)
  dataframe$col <- factor(dataframe$col, levels=unique(dataframe$col), ordered = T)
  dataframe$fac<-as.character(dataframe$fac)
  dataframe$fac <- factor(dataframe$fac, levels=unique(dataframe$fac), ordered = T)
  p<-ggplot(dataframe)
  p<-p+geom_violin(aes(x=dataframe$fac,y=dataframe$data,fill=dataframe$col),trim=F,scale="width",alpha=0.5,width=0.5)+geom_jitter(aes(x=dataframe$fac,y=dataframe$data),position = position_jitter(width = .3))
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
  
  pdf(paste(rownames(int_fpkm_glm)[i],"_gfp_mapping.pdf",sep=""),width=10,height=4)
  print(p)
  dev.off()
}



#Use only aNSC-mid and aNSC-late cells
comb_fpkm<-cbind(group3_fpkm,group4_fpkm)

#Genes to be used for correlation
comb_fpkm_curr<-comb_fpkm[match(c("Atp1a2","Ntsr2","Gja1","Jag1","Fgfr3","Dlx1","Dlx2","GFP","Prom1","Egfr"),rownames(comb_fpkm)),]

library(Hmisc)

#Calculate correlation
correlation_cellcycle<-rcorr(t(log(comb_fpkm_curr+1)),type="spearman")

my_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100)

#Plot correlations for Figure 3D
library(gplots)
pdf("Correlation_select_genes_heatmap_2_cluster_withgfp.pdf")
heatmap.2(correlation_cellcycle$r,col=my_palette,density.info="none",trace="none", keysize = 1)

dev.off()








comb_fpkm<-cbind(group1_fpkm,group2_fpkm,group3_fpkm,group4_fpkm,group5_fpkm)

comb_fpkm_col<-vector(length=length(comb_fpkm[1,]),mode="character")
comb_fpkm_col[match(group1,colnames(comb_fpkm))]<-"Group1"
comb_fpkm_col[match(group2,colnames(comb_fpkm))]<-"Group2"
comb_fpkm_col[match(group3,colnames(comb_fpkm))]<-"Group3"
comb_fpkm_col[match(group4,colnames(comb_fpkm))]<-"Group4"
comb_fpkm_col[match(group5,colnames(comb_fpkm))]<-"Group5"


#This pheno vector can be referenced by calling "color_by='type'" in any of the plotting 
#functions in monocle to color the cells plotted by their type.
pheno<-comb_fpkm_col

pheno.data.df <- data.frame(type=pheno)
rownames(pheno.data.df) <- colnames(comb_fpkm)

feature<-data.frame(rownames(comb_fpkm))


pd <- new("AnnotatedDataFrame", data = pheno.data.df)
fd <- new("AnnotatedDataFrame", data = feature)
HSMM_q <- newCellDataSet(comb_fpkm, phenoData = pd)

#This is the list of ordering genes which is used to order the qNSCs, all have known roles in the differentiation of neurons,
#and because we have the qNSCs and Astrocytes in this as well I will also use some genes associated with quiescence

curated_ordering_genes<-as.vector(read.table("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/10282015_names_in50models_3.txt")[,1])

HSMM_q <- setOrderingFilter(HSMM_q, curated_ordering_genes) # Set list of genes for ordering
HSMM_q <- reduceDimension(HSMM_q, use_irlba = F) # Reduce dimensionality
HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = T) # Order cells
if(as.character(pData(HSMM_q)[,1])[match(max(pData(HSMM_q)[,2]),pData(HSMM_q)[,2])]=="Group1"){
  HSMM_q <- orderCells(HSMM_q, num_paths = 1, reverse = F) # Order cells
}

library(scales)
pseudotime<-pData(HSMM_q)[,2]
pseudotime_mod<-rescale(pseudotime,to=c(0,100))
pData(HSMM_q)[,2]<-pseudotime_mod

source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")

#Plotting genes with respect to pseudotime, colored by putative population (qNSC-like, aNSC-early, aNSC-mid, aNSC-late, NPC-like)
my_genes<-c("Dlx2")
my_genes_ind <- match(my_genes,row.names(subset(fData(HSMM_q))))
cds_subset <- HSMM_q[my_genes_ind, ]
pdf("Dlx2_aggregated_ordering_top100_Figure2_groups.pdf",height=1.3,width=5.5)
p<-plot_genes_in_pseudotime_mod(cds_subset,color_by="type",cell_size=2,color_custom=c("#9966FF","#c9cc00","#ff9933","#ff3333","#00CCFF"),ncol=1,trend_formula=NULL)
p<-p+theme(axis.line=element_blank()) 
p<-p+theme(axis.title.y=element_blank()) 
p<-p+theme(axis.title.x=element_blank())
p<-p+theme(legend.position="none")
print(p)
dev.off()
