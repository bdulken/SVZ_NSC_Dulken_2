

#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
#Quality control for NS cells
#All NS cells - Miseq

#Read in counts matrix of all cells and extract neurosphere cells.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
allcounts<-read.table("All_spec_pops_counts_09222015_correctNames.txt",header=T)
allcounts_single<-allcounts[,grepl("NS_",colnames(allcounts))]

allcells_geneNum<-colSums(allcounts_single!= 0)   #Number of nonzero FPKM values for each cell
d=data.frame(Genes_Detected=allcells_geneNum[2:length(allcells_geneNum)])
p<-ggplot(d)
p<-p+geom_histogram(aes(x=d$Genes_Detected),smooth=T,binwidth=100, color="black",fill="red")
p<-p+geom_density(aes(x=d$Genes_Detected))
p<-p+theme_classic()
p<-p+labs(x="# of Genes Detected", title="Total Number of Genes Detected_AllCells")
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=18))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S4A
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
pdf("NumGenesDetect_neurosphere_miseq.pdf", width = 8, height=3.5)
print(p)
dev.off()

#Filter based on cells expressing 500 genes
allcells_counts_genes_single<-allcounts_single[,allcells_geneNum>500]

#----------------------------------------------------------------------------------------------------------------------------#
#Number of reads mapped to genes

#Read in counts matrix of all cells and extract neurosphere cells.
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Files")
allcounts<-read.table("All_spec_pops_counts_09222015_correctNames.txt",header=T)
allcounts_single<-allcounts[,grepl("NS_",colnames(allcounts))]

allcells_readsMap<-colSums(allcounts_single)   #Number of nonzero FPKM values for each cell
d=data.frame(Genes_Detected=allcells_readsMap[2:length(allcells_readsMap)])
p<-ggplot(d)
p<-p+geom_histogram(aes(x=d$Genes_Detected),smooth=T,binwidth=5000, color="black",fill="red")
p<-p+geom_density(aes(x=d$Genes_Detected))
p<-p+labs(x="# of Reads Mapped to Transcriptome", title="Total Number of Reads Mapped to transcriptome_AllCells")
p<-p+theme_classic()
p<-p+theme(axis.text.x=element_text(size=18))
p<-p+theme(axis.title.x=element_text(size=22))
p<-p+theme(axis.text.y=element_text(size=18))
p<-p+theme(axis.title.y=element_text(size=22))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.1))
p
#FIGURE S5A
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
pdf("NumReads_in_genes_Detect_neurosphere_miseq.pdf", width = 8, height=3.5)
print(p)
dev.off()


#Filter based on greater than 20,000 reads mapping to transcriptome
allcells_counts_reads_single<-allcounts_single[,allcells_readsMap>20000]

intersect_names<-intersect(colnames(allcells_counts_reads_single),colnames(allcells_counts_genes_single))

#filtering cells based on reads mapping to transcriptome and genes detected
allcounts_filtered<-allcounts[match(intersect_names,colnames(allcounts))]

#Write final counts matrix that is used for all subsequent analyses
setwd("/Volumes/guacamole/Software/Final Scripts for Paper _ 01152015/Figure5_Output")
#write.table(allcounts_filtered,"AllCounts_Neurosphere_read_gene_filt_FINAL.txt")
