
#Extracting Accuracy from model
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Machine Learning/NewOutput_5/")

library(ggplot2)

accuracy_training_resample<-vector()
accuracy_testing<-vector()
for(i in 1:100){
  myOptdGbm<-readRDS(paste("Classification_BalancedGroups_Accuracy_4groups_all_based_",i,"_model.R",sep=""))
  accuracy_training<-mean(myOptdGbm$resample[,1])
  accuracy_training_resample<-c(accuracy_training_resample,accuracy_training)
  training_conf<-as.numeric(read.table(paste("Classification_BalancedGroups_Accuracy_4groups_all_based_",i,"_testSetConfusionMatrixAndStats.txt",sep=""),fill=T,stringsAsFactors=F)[5,2])
  accuracy_testing<-c(accuracy_testing,training_conf)
}

names<-c(rep("Random",100),rep("Cross-vaildation\nTraining",100),rep("Testing",100))
values<-c(rep(0.25,100),accuracy_training_resample,accuracy_testing)

data1<-data.frame(names=names,values=values)
data1$names<-as.character(data1$names)
data1$names <- factor(data1$names, levels=unique(data1$names), ordered = T)
p<-ggplot(data1)
p<-p+geom_boxplot(aes(x=data1$names,y=data1$values),fill="green")
p<-p+theme_classic()
p<-p+labs(y="Accuracy",x="Condition")
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=18))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(plot.title=element_text(size=20))
p<-p+scale_y_continuous(limits=c(0,1))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Machine Learning")
pdf("Accuracy_GBM_original_10222015.pdf")
print(p)
dev.off()







#Class wise accuracy
#not currently used in paper
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure2_Output/Machine Learning/NewOutput_5/")
truevals<-c(rep("group1",20),
            rep("group2",20),
            rep("group3",20),
            rep("group4",20))
accuracy_training<-c()
for(i in 1:100){
myOptdGbm<-readRDS(paste("Classification_BalancedGroups_Accuracy_4groups_all_based_",i,"_model.R",sep=""))
trainingAssignments <- apply(myOptdGbm$finalModel$fit, 1, function(x){ colnames(myOptdGbm$finalModel$fit)[which.max(x)]})
trainingaccuracy<-sum(truevals==trainingAssignments)/80
accuracy_training<-c(accuracy_training,trainingaccuracy)
}