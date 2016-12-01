

# ====================================================================================================
# initial set up
# ====================================================================================================
rm(list=ls())
library(flux)
library(gbm)
#library(scales)
#library(fields)
library(caret)
library(doMC)
library(e1071)
library(edgeR)
library(monocle)
library(scde)
library("Hmisc")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_spanning_tree_mod.R")
source("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Functions/09282015_plot_genes_in_pseudotime_mod.txt")
registerDoMC(cores = 7)
set.seed(123)

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

#Reading in cells that are in the aNSC - cell cycle low and aNSC - cell cycle high states.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files/Ordering")
group1<-as.vector(read.table("group1_10252015_quiescent1_all_genes.txt")[,1])
group2<-as.vector(read.table("group2_10252015_active_nondividing_all_genes.txt")[,1])
group3<-as.vector(read.table("group3_10252015_active_dividing_all_genes.txt")[,1])

group1_fpkm<-fpkm_glm_genefilt_nooligo[,na.omit(match(group1,colnames(fpkm_glm_genefilt_nooligo)))]
group1_fpkm_col<-c(rep("group1",length(group1_fpkm[1,])))
group2_fpkm<-fpkm_glm_genefilt_nooligo[,na.omit(match(group2,colnames(fpkm_glm_genefilt_nooligo)))]
group2_fpkm_col<-c(rep("group2",length(group2_fpkm[1,])))
group3_fpkm<-fpkm_glm_genefilt_nooligo[,na.omit(match(group3,colnames(fpkm_glm_genefilt_nooligo)))]
group3_aNSC_fpkm<-group3_fpkm[,grepl("aNSC",colnames(group3_fpkm))]
group3_NPC_fpkm<-group3_fpkm[,grepl("NPC",colnames(group3_fpkm))]
group3_fpkm_aNSC_col<-c(rep("group3",length(group3_aNSC_fpkm[1,])))
group3_NPC_fpkm_col<-c(rep("group4",length(group3_NPC_fpkm[1,])))
test_fpkm<-cbind(group1_fpkm,group2_fpkm,group3_aNSC_fpkm,group3_NPC_fpkm)
for(i in 1:100){
  inTraining_group1<-sample(c(1:length(group1_fpkm[1,])),20)
  group1_fpkm_training<-group1_fpkm[,inTraining_group1]
  group1_fpkm_testing<-group1_fpkm[,-inTraining_group1]
  group1_fpkm_training_col<-group1_fpkm_col[inTraining_group1]
  group1_fpkm_testing_col<-group1_fpkm_col[-inTraining_group1]
  
  inTraining_group2<-sample(c(1:length(group2_fpkm[1,])),20)
  group2_fpkm_training<-group2_fpkm[,inTraining_group2]
  group2_fpkm_testing<-group2_fpkm[,-inTraining_group2]
  group2_fpkm_training_col<-group2_fpkm_col[inTraining_group2]
  group2_fpkm_testing_col<-group2_fpkm_col[-inTraining_group2]
  
  inTraining_group3<-sample(c(1:length(group3_aNSC_fpkm[1,])),20)
  group3_fpkm_training<-group3_aNSC_fpkm[,inTraining_group3]
  group3_fpkm_testing<-group3_aNSC_fpkm[,-inTraining_group3]
  group3_fpkm_training_col<-group3_fpkm_aNSC_col[inTraining_group3]
  group3_fpkm_testing_col<-group3_fpkm_aNSC_col[-inTraining_group3]
  
  inTraining_group4<-sample(c(1:length(group3_NPC_fpkm[1,])),20)
  group4_fpkm_training<-group3_NPC_fpkm[,inTraining_group4]
  group4_fpkm_testing<-group3_NPC_fpkm[,-inTraining_group4]
  group4_fpkm_training_col<-group3_NPC_fpkm_col[inTraining_group4]
  group4_fpkm_testing_col<-group3_NPC_fpkm_col[-inTraining_group4]
  
  
  Training_all<-cbind(group1_fpkm_training,group2_fpkm_training,group3_fpkm_training,group4_fpkm_training)
  Testing_all<-cbind(group1_fpkm_testing,group2_fpkm_testing,group3_fpkm_testing,group4_fpkm_testing)
  Training_all_col<-c(group1_fpkm_training_col,group2_fpkm_training_col,group3_fpkm_training_col,group4_fpkm_training_col)
  Testing_all_col<-c(group1_fpkm_testing_col,group2_fpkm_testing_col,group3_fpkm_testing_col,group4_fpkm_testing_col)
  
  #   Training_all<-cbind(group1_fpkm_training,group2_fpkm_training,group3_fpkm_training,group4_fpkm_training,group5_fpkm_training)
  #   Testing_all<-cbind(group1_fpkm_testing,group2_fpkm_testing,group3_fpkm_testing,group4_fpkm_testing,,group5_fpkm_testing)
  #   Training_all_col<-c(group1_fpkm_training_col,group2_fpkm_training_col,group3_fpkm_training_col,group4_fpkm_training_col,group5_fpkm_training_col)
  #   Testing_all_col<-c(group1_fpkm_testing_col,group2_fpkm_testing_col,group3_fpkm_testing_col,group4_fpkm_testing_col,group5_fpkm_testing_col)
  
  setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Machine Learning")
  training_data<-data.frame(t(Training_all),Training_all_col)
  testing_data<-data.frame(t(Testing_all),Testing_all_col)
  #   training_data<-data.frame(t(log(Training_all+1)),Training_all_col)
  #   testing_data<-data.frame(t(log(Testing_all+1)),Testing_all_col)
  write.table(rownames(training_data),paste("training_names_",i,".txt",sep=""))
  write.table(rownames(testing_data),paste("testing_names_",i,".txt",sep=""))
  
  #Train Model
  #my.ctrl.opt <-trainControl(method = "cv", number = 10, classProbs=TRUE, allowParallel=TRUE, verbose=T, summaryFunction=getMultiClassConfusionMatrixBalancedAccuracy)
  my.ctrl.opt <-trainControl(method = "cv", number = 5, classProbs=TRUE, allowParallel=TRUE, verbose=F)
  courseGbmGrid <- expand.grid(.interaction.depth = 5, .n.trees = 2500, .shrinkage = 0.001, .n.minobsinnode = 5 )
  
  # print("Print beginning to train GBM")
  # myOptdGbm <- train(x = as.data.frame(classifTraining[,-ncol(classifTraining)])
  #                      , y = classifTraining[,ncol(classifTraining)], method = "gbm", trControl = my.ctrl.opt
  #                      , tuneGrid = courseGbmGrid, verbose=FALSE)
  #   
  
  newBaseName<-paste("Classification_BalancedGroups_Accuracy_4groups_all_based_",i,sep="")
  myOptdGbm <- train(x = as.data.frame(training_data[,-ncol(training_data)]), y = as.factor(training_data[,ncol(training_data)])
                     , method = "gbm", trControl = my.ctrl.opt
                     , tuneGrid = courseGbmGrid)
  
  testingPr<-  predict(myOptdGbm$finalModel, type='response', newdata=testing_data[,-ncol(testing_data)], n.trees= myOptdGbm$finalModel$n.trees)
  
  assignmentPreds <- apply(testingPr, 1, function(x){ colnames(testingPr)[which.max(x)]})
  forConfMat <- as.factor(assignmentPreds)
  testConfMat <- confusionMatrix(forConfMat, testing_data[,ncol(testing_data)])
  for(ind in 1:length(testConfMat)){write.table(testConfMat[[ind]], file=paste(newBaseName, "testSetConfusionMatrixAndStats.txt", sep="_"), append=TRUE, quote=F, sep="\t", col.names=F)}
  
  caretedClassificationGbmFile<-paste(newBaseName,"_model.R",sep="")
  saveRDS(myOptdGbm, file = caretedClassificationGbmFile)
  
  importance <- summary(myOptdGbm$finalModel, n.trees=myOptdGbm$finalModel$n.trees)
  write.table(importance, file=paste(newBaseName,"featureImportance.tsv", sep="_"), quote=F, sep="\t")
}

importance<-as.vector(read.table("Classification_BalancedGroups_Accuracy_4groups_all_based_1_featureImportance.tsv")[,1])
for(i in 2:100){
  temp<-as.vector(read.table(paste("Classification_BalancedGroups_Accuracy_4groups_all_based_",i,"_featureImportance.tsv",sep=""))[,1])
  importance<-cbind(importance,temp)
}

setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Machine Learning")
write.table(importance,"importance_matrix_original_all_based.txt")