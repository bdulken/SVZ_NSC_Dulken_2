#Generation of Counts Matrices
#This file establishes the initial counts matrix for the data from Llorens-Bobadilla, et al.

rm(list=ls())
library(ggplot2)
library(edgeR)


# All spec pops ---------------------------------------------------------------------

#Writing the table of all counts for miseq 
setwd("/Volumes/guacamole/Analyses/Competitors_Dataset/HTseqcounts/new_name")
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

write.table(allcounts_fin,"/Volumes/guacamole/Analyses/Competitors_Dataset/Competitors_counts_allgenes.txt")
