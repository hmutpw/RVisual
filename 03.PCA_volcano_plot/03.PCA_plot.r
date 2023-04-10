rm(list=ls())
################################################
#plot PCA_result
################################################

PCA_res = read.table("./03.PCA_volcano_plot/01.PCA_data.txt",sep="\t",header=T,check.names=F)
row.names(PCA_res) = PCA_res$sample

#------setting my theme
library(ggplot2)
line_size = 0.5
text_size = 8
mytheme <- theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

#------plot PCA results.

p <- ggplot(PCA_res,aes(x=PC1,y=PC2,color = stage))+
geom_point(alpha=1,size=1)+
geom_text(aes(label=sample))+
mytheme

p

################################################
#how to do PCA
################################################

gene_counts = read.table("./01.calculate_TPM/gene_level_counts.tab",sep="\t",header=T,row.names=1,check.names=F)

#------get gene expression matrix
gene_mat = gene_counts[,-1]
gene_expressed = gene_mat[apply(gene_mat,1,function(x){length(x[which(x>1)])>2}),]

#------PCA analysis
pr = prcomp(t(gene_expressed), scale. = T)
summary(pr)
pca <- predict(pr)
