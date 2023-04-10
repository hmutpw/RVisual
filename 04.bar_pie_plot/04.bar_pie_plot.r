rm(list=ls())

test_res = read.table("./04.bar_pie_plot/bar_plot_data.txt",sep="\t",header=T,row.names=1,stringsAsFactors=F)
diff_res = test_res[test_res$is_diff!="not_diff",]
gene_type = table(diff_res[,c(2,8)])
gene_stat = as.data.frame(gene_type[rowSums(gene_type)>1,])
colnames(gene_stat) = c("gene_type","type","num")

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

plot_frm = gene_stat
#------bar plot
p = ggplot(data = plot_frm,aes(x=type,y=num, fill = gene_type))+
geom_bar(stat = "identity", position = "dodge")+
#geom_bar(stat = "identity", position = "stack")+
#geom_bar(stat = "identity", position = "fill")+
mytheme
p

#------coord flip
p = ggplot(data = plot_frm,aes(x=type,y=num, fill = gene_type))+
geom_bar(stat = "identity", position = "dodge")+
mytheme+
coord_flip()

p

#------pie plot

p = ggplot(data = plot_frm,aes(x="",y=num, fill = gene_type))+
geom_bar(stat="identity",width = 1, position="fill")+
coord_polar(theta = "y")+
mytheme+
facet_wrap(.~type,nrow=1)

p


