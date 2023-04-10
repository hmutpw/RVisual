data = read.table("./02.heatmap_line_plot/heatmap_data.txt",sep="\t",header=T,row.names=1)

#------plot heatmap
library(pheatmap)
plot_frm = log2(data+1)
p = pheatmap(plot_frm,scale="row",cutree_row=5,cluster_col=F)

hclust_res<-p$tree_row
cluster = cutree(hclust_res,k=5)
cluster_tab = data.frame(clusters = cluster)
cluster_tab$clusters = paste('cluster',cluster_tab$clusters,sep="")

p = pheatmap(plot_frm,
scale="row",
cutree_row=5,
gaps_col=seq(0,20,by=5),
cluster_col=F,
annotation_row = cluster_tab,
show_rownames=F,
treeheight_row=0)

gene_cluster = data.frame(cluster_tab,plot_frm[row.names(cluster_tab),])
#######################################################
#------line plot
#######################################################
#------prepare for line plot
stage_split = strsplit(colnames(data),split="_")
stage_vec = sapply(stage_split,function(x){x[1]})
names(stage_vec) = colnames(data)
stage = unique(stage_vec)

gene_exp_merge = t(apply(gene_cluster[-1],1,function(x){tapply(x,stage_vec[names(x)],median)}))
gene_merge_frm = data.frame(cluster = gene_cluster[row.names(gene_exp_merge),1],gene_exp_merge)


par(mfrow=c(3,2))
#-----Cluster1
ylim = c(0,12)
yseq = seq(0,12,by=4)
xlim = c(1:5)
plot_mat = gene_merge_frm[gene_cluster$clusters == "cluster1",-1]
for (i in row.names(plot_mat)){
plot(x=xlim, plot_mat[i,],type="l",col="lightgray",xaxt="n", yaxt="n", xlab="", ylab="", ylim=ylim)
par(new=T)
}
plot(x=xlim,y=apply(plot_mat,2,median),ylim=ylim,xaxt="n", yaxt="n",xlab = "Stage", ylab="log2(TPM+1)", col="red",lwd=2,pch=19,type="o")
axis(1,at=xlim,label=unique(stage))
axis(2,at=yseq, label=yseq)

#-----Cluster2
ylim = c(0,12)
yseq = seq(0,12,by=4)
xlim = c(1:5)
plot_mat = gene_merge_frm[gene_cluster$clusters == "cluster2",-1]
for (i in row.names(plot_mat)){
plot(x=xlim, plot_mat[i,],type="l",col="lightgray",xaxt="n", yaxt="n", xlab="", ylab="", ylim=ylim)
par(new=T)
}
plot(x=xlim,y=apply(plot_mat,2,median),ylim=ylim,xaxt="n", yaxt="n",xlab = "Stage", ylab="log2(TPM+1)", col="red",lwd=2,pch=19,type="o")
axis(1,at=xlim,label=unique(stage))
axis(2,at=yseq, label=yseq)

#-----Cluster3
ylim = c(0,12)
yseq = seq(0,12,by=4)
xlim = c(1:5)
plot_mat = gene_merge_frm[gene_cluster$clusters == "cluster3",-1]
for (i in row.names(plot_mat)){
plot(x=xlim, plot_mat[i,],type="l",col="lightgray",xaxt="n", yaxt="n", xlab="", ylab="", ylim=ylim)
par(new=T)
}
plot(x=xlim,y=apply(plot_mat,2,median),ylim=ylim,xaxt="n", yaxt="n",xlab = "Stage", ylab="log2(TPM+1)", col="red",lwd=2,pch=19,type="o")
axis(1,at=xlim,label=unique(stage))
axis(2,at=yseq, label=yseq)

#-----Cluster4
ylim = c(0,12)
yseq = seq(0,12,by=4)
xlim = c(1:5)
plot_mat = gene_merge_frm[gene_cluster$clusters == "cluster4",-1]
for (i in row.names(plot_mat)){
plot(x=xlim, plot_mat[i,],type="l",col="lightgray",xaxt="n", yaxt="n", xlab="", ylab="", ylim=ylim)
par(new=T)
}
plot(x=xlim,y=apply(plot_mat,2,median),ylim=ylim,xaxt="n", yaxt="n",xlab = "Stage", ylab="log2(TPM+1)", col="red",lwd=2,pch=19,type="o")
axis(1,at=xlim,label=unique(stage))
axis(2,at=yseq, label=yseq)

#-----Cluster5
ylim = c(0,12)
yseq = seq(0,12,by=4)
xlim = c(1:5)
plot_mat = gene_merge_frm[gene_cluster$clusters == "cluster5",-1]
for (i in row.names(plot_mat)){
plot(x=xlim, plot_mat[i,],type="l",col="lightgray",xaxt="n", yaxt="n", xlab="", ylab="", ylim=ylim)
par(new=T)
}
plot(x=xlim,y=apply(plot_mat,2,median),ylim=ylim,xaxt="n", yaxt="n",xlab = "Stage", ylab="log2(TPM+1)", col="red",lwd=2,pch=19,type="o")
axis(1,at=xlim,label=unique(stage))
axis(2,at=yseq, label=yseq)

#######################################################

#------how to change cluster order
hclust.result<-p$tree_row
cluster = cutree(hclust.result,k=cluster_num)
cluster_name = paste("C",cluster,sep="")
names(cluster_name) = names(cluster)
hclust.result.order<-p$tree_row$order
cluster_name_ordered = cluster_name[hclust.result.order]
cluster_categroy = unique(cluster_name_ordered)
new_cluster_name = paste("cluster",1:length(cluster_categroy),sep="")
names(new_cluster_name) = cluster_categroy



#------how to change colorbar
ramp.red <- colorRampPalette(c("#FF0000","#FF4500"))
ramp.white <- colorRampPalette(c("#FF4500","#FFFFFF","#1E90FF"))
ramp.blue <- colorRampPalette(c("#1E90FF","#0000FF"))
ramp.color = rev(c(ramp.red(80),ramp.white(320),ramp.blue(80)))



