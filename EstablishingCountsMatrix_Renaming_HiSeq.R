# All hiseq pops ---------------------------------------------------------------------

#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/G")
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

G_run1_counts<-allcounts_fin[!grepl("ERCC-",rownames(allcounts_fin)),]

#PG
#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/PG")
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

PG_run1_counts<-allcounts_fin[!grepl("ERCC-",rownames(allcounts_fin)),]


#PGE
#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/PGE")
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

PGE_run1_counts<-allcounts_fin[!grepl("ERCC-",rownames(allcounts_fin)),]





#Reading in E cells for run1
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1/E")
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

E_run1_counts<-allcounts_fin

All_run1_counts<-cbind(G_run1_counts,PG_run1_counts,PGE_run1_counts,E_run1_counts)


str_vec<-vector()
for(i in 1:length(colnames(All_run1_counts))){
  str<-""
 if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="PGE1"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"aNSC_A_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="PG1"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_A_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="PGE2"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"aNSC_B_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="PG2"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_B_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="PG"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_",num,sep="")
    str<-paste(str,"_Empty",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="PGE"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"aNSC_",num,sep="")
    str<-paste(str,"_Empty",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="G"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"Ast_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_run1_counts)[i],split=c("[_]"))[[1]][2]=="E"){
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"NPC_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else{
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_",num,sep="")
    str<-paste(str,"_Empty",sep="")
    str_vec<-c(str_vec,str)
  }
}
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1")
colnames(All_run1_counts)<-str_vec
write.table(All_run1_counts,"AllCounts_HiSeq_run1.txt")


#===========================================================================================================================================================
#===========================================================================================================================================================
#===========================================================================================================================================================
#===========================================================================================================================================================
#Run2

# All hiseq pops ---------------------------------------------------------------------

#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2/G")
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

G_Run2_counts<-allcounts_fin[!grepl("ERCC-",rownames(allcounts_fin)),]

#PG
#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2/PG")
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

PG_Run2_counts<-allcounts_fin[!grepl("ERCC-",rownames(allcounts_fin)),]

write.table(PG_Run2_counts,"PG_Run2_counts_test.txt")
#PGE
#Reading in all HTseq counts files.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2/PGE")
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

PGE_Run2_counts<-allcounts_fin[!grepl("ERCC-",rownames(allcounts_fin)),]





All_Run2_counts<-cbind(G_Run2_counts,PG_Run2_counts,PGE_Run2_counts)


str_vec<-vector()
for(i in 1:length(colnames(All_Run2_counts))){
  str<-""
  if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="PGE1"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"aNSC_A_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="PG1"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_A_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="PGE2"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"aNSC_B_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="PG2"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_B_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="PG"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_",num,sep="")
    str<-paste(str,"_Empty",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="PGE"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"aNSC_",num,sep="")
    str<-paste(str,"_Empty",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="G"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"Ast_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else if(strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))[[1]][2]=="E"){
    filename<-strsplit(colnames(All_Run2_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"NPC_",num,sep="")
    str<-paste(str,"_1g",sep="")
    str_vec<-c(str_vec,str)
  }else{
    filename<-strsplit(colnames(All_run1_counts)[i],split=c("[_]"))
    num<-gsub("C","",filename[[1]][1])
    str<-paste(str,"qNSC_",num,sep="")
    str<-paste(str,"_Empty",sep="")
    str_vec<-c(str_vec,str)
  }
}

test<-rbind(str_vec,colnames(All_Run2_counts))

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2")
colnames(All_Run2_counts)<-str_vec
write.table(All_Run2_counts,"AllCounts_HiSeq_Run2.txt")






#High Quality HiSeq Cells

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run1")
allcounts_run1<-read.table('AllCounts_HiSeq_Run1.txt')
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/AllCounts_HiSeq/Run2")
allcounts_run2<-read.table('AllCounts_HiSeq_Run2.txt')
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/")
spec_pops<-read.table("AllCounts_specPops_read_gene_ERCC_filt_FINAL.txt")


allcounts_run1_highqual<-allcounts_run1[,!is.na(match(colnames(allcounts_run1),colnames(spec_pops)))]
allcounts_run2_highqual<-allcounts_run2[,!is.na(match(colnames(allcounts_run2),colnames(spec_pops)))]

write.table(allcounts_run1_highqual,"AllCounts_HiSeq_Run1_highqual.txt")

write.table(allcounts_run2_highqual,"AllCounts_HiSeq_Run2_highqual.txt")
