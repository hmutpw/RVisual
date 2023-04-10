rm(list=ls())
#------import data
data = read.table("./06.density_plot/violin_data.txt",sep="\t",header=T,row.names=1)
sample2stage = read.table("./06.density_plot/sample_infor.txt",sep="\t",header=T,row.names=1)

#------set mytheme
library(ggplot2)
text_size = 8
line_size = 0.5
mytheme <- theme(line = element_line(size=line_size), text = element_text(size=text_size, color="black"), 
panel.border = element_rect(size = line_size, fill=NA),
panel.background = element_rect(fill=NA),
panel.grid.major = element_line(size=0.25*line_size, linetype = 8, color = "lightgray"),
axis.ticks = element_line(size=line_size),
axis.title.x = element_text(size=text_size),
axis.text.x= element_text(size=text_size, angle=45,hjust=1,vjust=1,color="black"), 
axis.text.y= element_text(size=text_size, color="black"),
legend.text = element_text(size = text_size), 
legend.title = element_text(size=text_size*1.2),
plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5),legend.position='right')

#------plot my data
expressed_gene_num = apply(data,2,function(x){length(which(x>1))})
expressed_gene_num_tab = data.frame(sample_stage,num = expressed_gene_num[row.names(sample_stage)])

stages = c(paste(6:9,"w",sep=""),"11w","14w","16w","19w","21w","24w","25w")
expressed_gene_num_tab$stage = factor(expressed_gene_num_tab$stage,levels = stages)
expressed_gene_num_tab$sample = factor(expressed_gene_num_tab$sample,levels = c("esophagus","stomach","SI","LI"))

library(ggridges)
plot_frm = expressed_gene_num_tab

#------get histogram distribution
p = ggplot(data = plot_frm, aes(x = num))+
geom_histogram(binwidth=10,fill="pink")+
labs(x='gene num',y='cell number',title = '')+
mytheme+
facet_wrap(.~sample,nrow=2)
p

#------get density distribution
p = ggplot(data = plot_frm, aes(x = num,fill=stage))+
geom_density(alpha=0.5)+
#scale_color_manual(values=heavy4color)+
#scale_fill_manual(values=light4color)+
labs(x='gene num',y='',title = '')+
mytheme+
facet_wrap(.~sample,nrow=2)
p

#------get split density distribution

p = ggplot(data = plot_frm, aes(x = num, y = stage, fill = stage, height=..density..))+
geom_density_ridges_gradient(size = line_size, scale = 1.4, panel_scaling = F, alpha=.8,na.rm=T)+
#scale_color_manual(values=heavy4color)+
#scale_fill_manual(values=rev(light4color))+
labs(x='cor',y='',title = '')+
mytheme+
facet_wrap(.~sample,nrow=2)
p

#########################################

#########################################

library(plyr)
plot_ecdf = ddply(plot_frm[,c("sample", "stage","num")], .(sample,stage), transform, ecdf = ecdf(num)(num) )

text_size = 8
line_size = 0.5
p = ggplot(plot_ecdf,aes(x = num, y=ecdf, color = sample))+
stat_ecdf(size = line_size)+
mytheme+
facet_wrap(.~stage,nrow=5)

p

