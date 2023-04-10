rm(list=ls())

data = read.table("./05.violin_plot/violin_data.txt",sep="\t",header=T,row.names=1)
sample2stage = read.table("./05.violin_plot/sample_infor.txt",sep="\t",header=T,row.names=1)

#------get expressed gene number
expressed_gene_num = apply(data,2,function(x){length(which(x>1))})
expressed_gene_num_tab = data.frame(sample2stage,num = expressed_gene_num[row.names(sample2stage)])

stages = c(paste(6:9,"w",sep=""),"11w","14w","16w","19w","21w","24w","25w")
expressed_gene_num_tab$stage = factor(expressed_gene_num_tab$stage,levels = stages)
expressed_gene_num_tab$sample = factor(expressed_gene_num_tab$sample,levels = c("esophagus","stomach","SI","LI"))

#------get my theme
text_size = 8
line_size = 0.5/1.07
library(ggplot2)
#------set plot theme
mytheme = theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA,color=NA),
strip.background = element_rect(fill=NA,color=NA),
strip.text = element_text(size=text_size,color="black"),
panel.grid.major = element_line(size=line_size/2, linetype = 8, color = "lightgray"),
panel.grid.minor = element_line(size=NA, color = NA),
axis.ticks = element_line(size=line_size,color="black"),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45, hjust = 1, vjust =1,color="black"), 
axis.text.y= element_text(size=text_size,color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

p <- ggplot(plot_frm,mapping = aes(x=Description,y=-LogP))+
geom_bar(stat="identity",width=0.7,fill="blue")+
mytheme+
coord_flip()
p



plot_frm = expressed_gene_num_tab

p <- ggplot(plot_frm,mapping = aes(x=stage,y=num, color=stage))+
geom_violin(position = "dodge", na.rm = TRUE, scale = "width",width = 0.9)+
geom_boxplot(width=0.25,fill="white", size=0.2,outlier.colour=NA, alpha = 0.8)+
geom_jitter(width=0.25,size=0.8,shape=21,fill="white", alpha = 0.8)+
labs(x="",y='number of expressed genes in each single cells')+
mytheme+
facet_wrap(.~sample,nrow=4)

p

