rm(list=ls())

test_res = read.table("./03.PCA_volcano_plot/02.volcano_plot_data.txt",sep="\t",header=T,row.names=1)

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

#------construct ggplot data.
plot_frm = test_res[,c("logFC","PValue","is_diff")]

#------plot

demo_gene = c("CD24","PSG3")
names(demo_gene) = c("ENSG00000272398","ENSG00000221826")
demo_gene_tab = plot_frm[names(demo_gene),]

p <- ggplot(plot_frm,aes(x = logFC,y = -log10(PValue),color = is_diff))+
geom_point(size=1)+
geom_text(data= demo_gene_tab, aes(x=logFC, y=-log10(PValue)),label=demo_gene,color = "black")+
xlim(-12,12)+
scale_color_manual(values=c("lightblue","darkgray","pink"))+
mytheme

p


