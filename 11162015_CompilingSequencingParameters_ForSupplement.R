rm(list=ls())

setwd("/Volumes/guacamole/Analyses/10132014_SingleCell_inVivo_rep2_PG_PGE/10222014_QNSC/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

qNSC_A<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(qNSC_A,"qNSC_A_mapping.txt")

#==============================================================================================================================#
#==============================================================================================================================#
#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/10132014_SingleCell_inVivo_rep2_PG_PGE/10222014_ANSC/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

aNSC_A<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(aNSC_A,"aNSC_A_mapping.txt")

#==============================================================================================================================#
#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E/G_RawReads/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

G<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(G,"G_mapping.txt")


#==============================================================================================================================#
#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E/PG_RawReads/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PG<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PG,"PG_mapping.txt")

#==============================================================================================================================#
#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E/PGE_RawReads/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PGE<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PGE,"PGE_mapping.txt")

#==============================================================================================================================#
#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E/E_RawReads/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

E<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(E,"E_mapping.txt")

#==============================================================================================================================#
#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/01202015_SingleCell_inVivo_PGE_spikeIns/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="_")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PGE_C<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PGE_C,"PGE_C_mapping.txt")

#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/02102015_SingleCell_inVivo_PG_PGE_Mid_spikeIns/QNSC_spikeIn/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="[.]")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PG_C<-rbind(name,input,mapped)
setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PG_C,"PG_C_mapping.txt")


#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/02102015_SingleCell_inVivo_PG_PGE_Mid_spikeIns/ANSC_spikeIn/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="[.]")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PGE_D<-rbind(name,input,mapped)

setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PGE_D,"PGE_D_mapping.txt")


#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/06302014_Single_Cell_RNA_seq_neurosphere/Data/STAR/")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="[.]")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

NS<-rbind(name,input,mapped)

setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(NS,"NS_mapping.txt")

#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E_HighThroughput/Run_1/PGE/STAR")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="[.]")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PGE_1<-rbind(name,input,mapped)

setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PGE_1,"PGE_Hiseq_run1_mapping.txt")


#==============================================================================================================================#

setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E_HighThroughput/Run_2/PGE/STAR")

fileList=list.files(pattern=".out$")

name<-c()
input<-c()
mapped<-c()

for(i in 1:length(fileList)){
  stats<-read.table(fileList[i],fill=T,sep="\t",stringsAsFactors=F)
  name_t<-strsplit(fileList[i],split="[.]")[[1]][1]
  name<-c(name,name_t)
  input_t<-as.numeric(stats[5,2])
  input<-c(input,input_t)
  mapped_t<-as.numeric(stats[8,2])
  mapped<-c(mapped,mapped_t)
}

PGE_2<-rbind(name,input,mapped)

setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
write.table(PGE_2,"PGE_Hiseq_run2_mapping.txt")







#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
counts_NS<-read.table("AllCounts_Neurosphere_read_gene_filt_FINAL.txt")

allcounts_allcells<-spec_pops

counts_aNSC_A<-allcounts_allcells[grep("aNSC_A",colnames(allcounts_allcells))]
counts_aNSC_B<-allcounts_allcells[grep("aNSC_B",colnames(allcounts_allcells))]
counts_aNSC_C<-allcounts_allcells[grep("aNSC_C",colnames(allcounts_allcells))]
counts_aNSC_D<-allcounts_allcells[grep("aNSC_D",colnames(allcounts_allcells))]
counts_qNSC_A<-allcounts_allcells[grep("qNSC_A",colnames(allcounts_allcells))]
counts_qNSC_B<-allcounts_allcells[grep("qNSC_B",colnames(allcounts_allcells))]
counts_qNSC_C<-allcounts_allcells[grep("qNSC_C",colnames(allcounts_allcells))]
counts_Ast<-allcounts_allcells[grep("Ast",colnames(allcounts_allcells))]
counts_NPC<-allcounts_allcells[grep("NPC",colnames(allcounts_allcells))]
counts_aNSC_A_B<-cbind(counts_aNSC_A,counts_aNSC_B)

setwd("/Volumes/guacamole/Analyses/09012015_Figure1/")
qNSC_A<-read.table("qNSC_A_mapping.txt",stringsAsFactors=F)
aNSC_A<-read.table("aNSC_A_mapping.txt",stringsAsFactors=F)
Ast<-read.table("G_mapping.txt",stringsAsFactors=F)
qNSC_B<-read.table("PG_mapping.txt",stringsAsFactors=F)
aNSC_B<-read.table("PGE_mapping.txt",stringsAsFactors=F)
NPC<-read.table("E_mapping.txt",stringsAsFactors=F)
qNSC_C<-read.table("PG_C_mapping.txt",stringsAsFactors=F)
qNSC_C<-qNSC_C[,!grepl("M",as.vector(qNSC_C[1,],mode="character"))]
aNSC_C<-read.table("PGE_C_mapping.txt",stringsAsFactors=F)
aNSC_D<-read.table("PGE_D_mapping.txt",stringsAsFactors=F)
aNSC_D<-aNSC_D[,!grepl("M",as.vector(aNSC_D[1,],mode="character"))]
NS<-read.table("NS_mapping.txt",stringsAsFactors=F)
PGE_1<-read.table("PGE_Hiseq_run1_mapping.txt",stringsAsFactors=F)
PGE_2<-read.table("PGE_Hiseq_run2_mapping.txt",stringsAsFactors=F)

#aNSC_A
split<-strsplit(colnames(counts_aNSC_A),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
aNSC_A_nums<-c()
aNSC_A_names<-as.vector(aNSC_A[1,],mode="character")
for(i in 1:length(aNSC_A_names)){
  num<-gsub("C","",strsplit(aNSC_A_names[i],split="[-]")[[1]][3])
  aNSC_A_nums<-c(aNSC_A_nums,num)
}
aNSC_A_specs<-aNSC_A[,match(nums,aNSC_A_nums)]
aNSC_A_final_specs<-rbind(aNSC_A_specs,colnames(counts_aNSC_A))
write.table(aNSC_A_final_specs,"aNSC_A_specs.txt")



#qNSC_A
split<-strsplit(colnames(counts_qNSC_A),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
qNSC_A_nums<-c()
qNSC_A_names<-as.vector(qNSC_A[1,],mode="character")
for(i in 1:length(qNSC_A_names)){
  num<-gsub("C","",strsplit(qNSC_A_names[i],split="[-]")[[1]][3])
  qNSC_A_nums<-c(qNSC_A_nums,num)
}
qNSC_A_specs<-qNSC_A[,match(nums,qNSC_A_nums)]
qNSC_A_final_specs<-rbind(qNSC_A_specs,colnames(counts_qNSC_A))
write.table(qNSC_A_final_specs,"qNSC_A_specs.txt")




#Ast
split<-strsplit(colnames(counts_Ast),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][2])
}
Ast_nums<-c()
Ast_names<-as.vector(Ast[1,],mode="character")
for(i in 1:length(Ast_names)){
  num<-gsub("C","",strsplit(Ast_names[i],split="[-]")[[1]][3])
  Ast_nums<-c(Ast_nums,num)
}
Ast_specs<-Ast[,match(nums,Ast_nums)]
Ast_final_specs<-rbind(Ast_specs,colnames(counts_Ast))
write.table(Ast_final_specs,"Ast_specs.txt")



#qNSC_B
split<-strsplit(colnames(counts_qNSC_B),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
qNSC_B_nums<-c()
qNSC_B_names<-as.vector(qNSC_B[1,],mode="character")
for(i in 1:length(qNSC_B_names)){
  num<-gsub("C","",strsplit(qNSC_B_names[i],split="[-]")[[1]][3])
  qNSC_B_nums<-c(qNSC_B_nums,num)
}
qNSC_B_specs<-qNSC_B[,match(nums,qNSC_B_nums)]
qNSC_B_final_specs<-rbind(qNSC_B_specs,colnames(counts_qNSC_B))
write.table(qNSC_B_final_specs,"qNSC_B_specs.txt")

#aNSC_B
split<-strsplit(colnames(counts_aNSC_B),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
aNSC_B_nums<-c()
aNSC_B_names<-as.vector(aNSC_B[1,],mode="character")
for(i in 1:length(aNSC_B_names)){
  num<-gsub("C","",strsplit(aNSC_B_names[i],split="[-]")[[1]][3])
  aNSC_B_nums<-c(aNSC_B_nums,num)
}
aNSC_B_specs<-aNSC_B[,match(nums,aNSC_B_nums)]
aNSC_B_final_specs<-rbind(aNSC_B_specs,colnames(counts_aNSC_B))
write.table(aNSC_B_final_specs,"aNSC_B_specs.txt")

#NPC
split<-strsplit(colnames(counts_NPC),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][2])
}
NPC_nums<-c()
NPC_names<-as.vector(NPC[1,],mode="character")
for(i in 1:length(NPC_names)){
  num<-gsub("C","",strsplit(NPC_names[i],split="[-]")[[1]][3])
  NPC_nums<-c(NPC_nums,num)
}
NPC_specs<-NPC[,match(nums,NPC_nums)]
NPC_final_specs<-rbind(NPC_specs,colnames(counts_NPC))
write.table(NPC_final_specs,"NPC_specs.txt")


#aNSC_C
split<-strsplit(colnames(counts_aNSC_C),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
aNSC_C_nums<-c()
aNSC_C_names<-as.vector(aNSC_C[1,],mode="character")
for(i in 1:length(aNSC_C_names)){
  num<-strsplit(aNSC_C_names[i],split="[-]")[[1]][2]
  aNSC_C_nums<-c(aNSC_C_nums,num)
}
aNSC_C_specs<-aNSC_C[,match(nums,aNSC_C_nums)]
aNSC_C_final_specs<-rbind(aNSC_C_specs,colnames(counts_aNSC_C))
write.table(aNSC_C_final_specs,"aNSC_C_specs.txt")

#qNSC_C
split<-strsplit(colnames(counts_qNSC_C),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
qNSC_C_nums<-c()
qNSC_C_names<-as.vector(qNSC_C[1,],mode="character")
for(i in 1:length(qNSC_C_names)){
  num<-strsplit(qNSC_C_names[i],split="[_]")[[1]][2]
  qNSC_C_nums<-c(qNSC_C_nums,num)
}
qNSC_C_specs<-qNSC_C[,match(nums,qNSC_C_nums)]
qNSC_C_final_specs<-rbind(qNSC_C_specs,colnames(counts_qNSC_C))
write.table(qNSC_C_final_specs,"qNSC_C_specs.txt")


#aNSC_D
split<-strsplit(colnames(counts_aNSC_D),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][3])
}
aNSC_D_nums<-c()
aNSC_D_names<-as.vector(aNSC_D[1,],mode="character")
for(i in 1:length(aNSC_D_names)){
  num<-strsplit(aNSC_D_names[i],split="[_]")[[1]][2]
  aNSC_D_nums<-c(aNSC_D_nums,num)
}
aNSC_D_specs<-aNSC_D[,match(nums,aNSC_D_nums)]
aNSC_D_final_specs<-rbind(aNSC_D_specs,colnames(counts_aNSC_D))
write.table(aNSC_D_final_specs,"aNSC_D_specs.txt")

#NS
split<-strsplit(colnames(counts_NS),split="_")
nums<-c()
for(i in 1:length(split)){
  nums<-c(nums,split[[i]][2])
}
NS_nums<-c()
NS_names<-as.vector(NS[1,],mode="character")
for(i in 1:length(NS_names)){
  num<-gsub("C","",strsplit(strsplit(NS_names[i],split="[_]")[[1]][1],split="[-]")[[1]][3])
  NS_nums<-c(NS_nums,num)
}
NS_specs<-NS[,match(nums,NS_nums)]
NS_final_specs<-rbind(NS_specs,colnames(counts_NS))
write.table(NS_final_specs,"NS_specs.txt")

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
final_specs<-cbind(Ast_final_specs,qNSC_A_final_specs,qNSC_B_final_specs,qNSC_C_final_specs,aNSC_A_final_specs,aNSC_B_final_specs,aNSC_C_final_specs,aNSC_D_final_specs,NPC_final_specs,NS_final_specs)
write.table(t(final_specs),"final_specs_withOligos_withNS.txt")




#PGE1
PGE1_newnames<-c()
PGE1_names<-as.vector(PGE_1[1,],mode="character")
for(i in 1:length(PGE1_names)){
  if(strsplit(PGE1_names[i],split="[_]")[[1]][2]=="PGE1"){
    num<-gsub("C","",strsplit(PGE1_names[i],split="[_]")[[1]][1])
    temp_name<-paste("aNSC_A_",num,"_1g",sep="")
    PGE1_newnames<-c(PGE1_newnames,temp_name)
  }else if(strsplit(PGE1_names[i],split="[_]")[[1]][2]=="PGE2"){
    num<-gsub("C","",strsplit(PGE1_names[i],split="[_]")[[1]][1])
    temp_name<-paste("aNSC_B_",num,"_1g",sep="")
    PGE1_newnames<-c(PGE1_newnames,temp_name)
  }else{
    PGE1_newnames<-c(PGE1_newnames,"Empty")
  }
}

PGE1_hiseq_newnames<-paste(PGE1_newnames,"_hiseq_run1",sep="")
PGE2_hiseq_newnames<-paste(PGE1_newnames,"_hiseq_run2",sep="")

PGE_hiseq1_quals<-PGE_1[,!is.na(match(PGE1_newnames,colnames(counts_aNSC_A_B)))]
PGE_hiseq2_quals<-PGE_2[,!is.na(match(PGE1_newnames,colnames(counts_aNSC_A_B)))]

PGE_hiseq1_quals_names<-PGE1_hiseq_newnames[!is.na(match(PGE1_newnames,colnames(counts_aNSC_A_B)))]
PGE_hiseq2_quals_names<-PGE2_hiseq_newnames[!is.na(match(PGE1_newnames,colnames(counts_aNSC_A_B)))]

PGE_hiseq1_quals<-rbind(PGE_hiseq1_quals,PGE_hiseq1_quals_names)
PGE_hiseq2_quals<-rbind(PGE_hiseq2_quals,PGE_hiseq2_quals_names)

PGE_hiseq1_quals_final<-cbind(PGE_hiseq1_quals[,grep("aNSC_A",PGE_hiseq1_quals[4,])],PGE_hiseq1_quals[,grep("aNSC_B",PGE_hiseq1_quals[4,])])
PGE_hiseq2_quals_final<-cbind(PGE_hiseq2_quals[,grep("aNSC_A",PGE_hiseq2_quals[4,])],PGE_hiseq2_quals[,grep("aNSC_B",PGE_hiseq2_quals[4,])])

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
final_specs<-cbind(Ast_final_specs,qNSC_A_final_specs,qNSC_B_final_specs,qNSC_C_final_specs,aNSC_A_final_specs,aNSC_B_final_specs,aNSC_C_final_specs,aNSC_D_final_specs,NPC_final_specs,NS_final_specs,PGE_hiseq1_quals_final,PGE_hiseq2_quals_final)
write.table(t(final_specs),"final_specs_withOligos_withNS_withHiSeq.txt")






#Getting Groups
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10282015_50minModels_quiescent1_consensus_2.txt")[,1])
group2<-as.vector(read.table("group2_10282015_50minModels_active_nondividing_consensus_2.txt")[,1])
group3<-as.vector(read.table("group3_10282015_50minModels_active_dividing_early_consensus_2.txt")[,1])
group4<-as.vector(read.table("group4_10282015_50minModels_active_dividing_late_consensus_2.txt")[,1])
group5<-as.vector(read.table("group5_10282015_50minModels_active_dividing_consensus_2.txt")[,1])

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])

final_specs_miseq<-cbind(Ast_final_specs,qNSC_A_final_specs,qNSC_B_final_specs,qNSC_C_final_specs,
                         aNSC_A_final_specs,aNSC_B_final_specs,aNSC_C_final_specs,aNSC_D_final_specs,NPC_final_specs,NS_final_specs)

group_type<-c(rep("Outlier",length(final_specs_miseq[1,])))
group_type[grepl("Ast",final_specs_miseq[4,])]<-"Astrocyte"
group_type[charmatch(oligos,final_specs_miseq[4,])]<-"Oligodendrocyte - Like Cell"
group_type[charmatch(group1,final_specs_miseq[4,])]<-"qNSC-like"
group_type[charmatch(group2,final_specs_miseq[4,])]<-"aNSC-early"
group_type[charmatch(group3,final_specs_miseq[4,])]<-"aNSC-mid"
group_type[charmatch(group4,final_specs_miseq[4,])]<-"aNSC-late"
group_type[charmatch(group5,final_specs_miseq[4,])]<-"NPC-like"
group_type[grepl("NS_",final_specs_miseq[4,])]<-"Passage 3 - Neurosphere"

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
final_specs_miseq_2<-rbind(final_specs_miseq,group_type)
write.table(t(final_specs_miseq_2),"final_specs_withOligos_withNS_miseq_withgroups.txt")

final_specs_hiseq<-cbind(PGE_hiseq1_quals_final)

group_type<-c(rep("Outlier",length(final_specs_hiseq[1,])))
group_type[grepl("Ast",final_specs_hiseq[4,])]<-"Astrocyte"
group_type[charmatch(oligos,final_specs_hiseq[4,])]<-"Oligodendrocyte - Like Cell"
group_type[charmatch(group1,final_specs_hiseq[4,])]<-"qNSC-like"
group_type[charmatch(group2,final_specs_hiseq[4,])]<-"aNSC-early"
group_type[charmatch(group3,final_specs_hiseq[4,])]<-"aNSC-mid"
group_type[charmatch(group4,final_specs_hiseq[4,])]<-"aNSC-late"
group_type[charmatch(group5,final_specs_hiseq[4,])]<-"NPC-like"
group_type[grepl("NS_",final_specs_hiseq[4,])]<-"Passage 3 - Neurosphere"

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
final_specs_hiseq_2<-rbind(final_specs_hiseq,group_type)
write.table(t(final_specs_hiseq_2),"final_specs_withOligos_withNS_hiseq_withgroups.txt")

#GeneratingFPKM tables


#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
counts_NS<-read.table("AllCounts_Neurosphere_read_gene_filt_FINAL.txt")


allcounts_allcells<-cbind(spec_pops,counts_NS)
final_specs<-cbind(Ast_final_specs,qNSC_A_final_specs,qNSC_B_final_specs,qNSC_C_final_specs,aNSC_A_final_specs,aNSC_B_final_specs,aNSC_C_final_specs,aNSC_D_final_specs,NPC_final_specs,NS_final_specs)

allcounts_allcells_ordered<-allcounts_allcells[,match(final_specs[4,],colnames(allcounts_allcells))]

write.table(allcounts_allcells_ordered,"AllCounts_ORDERED_forSupplement.txt")

greaterthan0<-allcounts_allcells_ordered>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_ordered_genefilt<-allcounts_allcells_ordered[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_ordered_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Annotations/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt <- rpkm(glm.norm, gene.length=genes$exonic.size,
                                  normalized.lib.sizes=T, log=F)

write.table(fpkm_glm_genefilt,"AllFPKM_ORDERED_forSupplement.txt")





#Loading all high quality cells and filtering for lowly expressed genes and cell cycle genes
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")

allcounts_allcells<-cbind(spec_pops,counts_NS)
final_specs<-cbind(Ast_final_specs,qNSC_A_final_specs,qNSC_B_final_specs,qNSC_C_final_specs,aNSC_A_final_specs,aNSC_B_final_specs,aNSC_C_final_specs,aNSC_D_final_specs,NPC_final_specs)

allcounts_allcells_ordered<-allcounts_allcells[,match(final_specs[4,],colnames(allcounts_allcells))]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
oligos<-as.vector(read.table("STAR_oligos_updated_09232015.txt")[,1])
allcounts_allcells_nooligo<-allcounts_allcells_ordered[,-na.omit(match(oligos,colnames(allcounts_allcells_ordered)))]

greaterthan0<-allcounts_allcells_nooligo>10
greaterthan0sum<-rowSums(greaterthan0)
allcounts_allcells_nooligo_genefilt<-allcounts_allcells_nooligo[greaterthan0sum>=5,]

glm <- DGEList(counts=allcounts_allcells_nooligo_genefilt)
glm.norm <- calcNormFactors(glm,method="TMM")
load("/Volumes/guacamole/Software/Annotations/Katja_exonic_gene_sizes_ERCC92.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(tolower(rownames(glm.norm)), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
genes <- data.frame(gene.symbol=rownames(glm.norm),
                    exonic.size=exonic.gene.sizes.ord)
fpkm_glm_genefilt <- rpkm(glm.norm, gene.length=genes$exonic.size,
                          normalized.lib.sizes=T, log=F)

write.table(fpkm_glm_genefilt,"AllFPKM_ORDERED_forSupplement_noOutliers.txt")






#==============================================================================================================================#
#==============================================================================================================================#
#==============================================================================================================================#
#==============================================================================================================================#
#Counts matrices for hiseq cells

#Writing the table of all counts for miseq 
setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E_HighThroughput/Run_1/PGE/HTseqcounts")
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

run1_counts<-allcounts_fin[,charmatch(PGE_hiseq1_quals_final[1,],colnames(allcounts_fin))]

PGE_hiseq1_quals_final[4,]

colnames(run1_counts)<-PGE_hiseq1_quals_final[4,]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
write.table(run1_counts,"12082015_hiseq_run1_PGE_hiqual_counts.txt")


#==============================================================================================================================#
#Counts matrices for hiseq cells

#Writing the table of all counts for miseq 
setwd("/Volumes/guacamole/Analyses/11132014_SingleCell_inVivo_G_PG_PGE_E_HighThroughput/Run_2/PGE/HTseqcounts")
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

run2_counts<-allcounts_fin[,charmatch(PGE_hiseq2_quals_final[1,],colnames(allcounts_fin))]

PGE_hiseq2_quals_final[4,]

colnames(run2_counts)<-PGE_hiseq2_quals_final[4,]

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
write.table(run2_counts,"12082015_hiseq_run2_PGE_hiqual_counts.txt")
