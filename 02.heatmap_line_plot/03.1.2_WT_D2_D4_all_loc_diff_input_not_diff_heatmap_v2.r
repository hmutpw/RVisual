rm(list=ls())
gc()
library(data.table)
#######################################
#get stage diff genelist.
#######################################
library(data.table)
#------read a file and get diff or non-diff results.
read_diff <- function(stagepare, inPath, return = c("diff", "not_diff"), 
                      diff_fdr = 0.05, diff_logfc = 1, not_diff_fdr = 0.2){
  return <- match.arg(return)
  file = file.path(inPath,paste(stagepare,"all_gene_edgeR_res.txt",sep="_"))
  tab <- fread(file,sep="\t",header=T,data.table=F)
  row.names(tab) <- tab$gene_id
  if(return=="diff"){
    diff_res <- tab[tab$FDR<diff_fdr,]
    out_res <- diff_res[abs(diff_res$logFC)>diff_logfc,]
  }else if(return=="not_diff"){
    out_res <- tab[tab$FDR>not_diff_fdr,]
  }
  out_res
}
#---get differential expressed gene
get_diff_gene <- function(inPath, diff_fdr = 0.05, diff_logfc = 1, not_diff_fdr = 0.2){
  #------input
  input_dir <- file.path(inPath,"input","all_res")
  WT_D2_input_not_diff = read_diff(stagepare="WT_D2", inPath=input_dir, return="not_diff", 
                                   diff_fdr=diff_fdr, diff_logfc=diff_logfc, 
                                   not_diff_fdr=not_diff_fdr)
  D2_D4_input_not_diff = read_diff(stagepare = "D2_D4", inPath=input_dir, return="not_diff", 
                                   diff_fdr=diff_fdr, diff_logfc=diff_logfc, 
                                   not_diff_fdr=not_diff_fdr)
  WT_D4_input_not_diff = read_diff(stagepare = "WT_D4", inPath=input_dir, return="not_diff", 
                                   diff_fdr=diff_fdr, diff_logfc=diff_logfc, 
                                   not_diff_fdr=not_diff_fdr)
  WT_D2_D4_input_not_diff <- Reduce(intersect,list(row.names(WT_D2_input_not_diff),
                                                   row.names(D2_D4_input_not_diff),
                                                   row.names(WT_D4_input_not_diff)))
  #------CE
  CE_dir <- file.path(inPath,"CE","all_res")
  WT_D2_CE_diff = read_diff(stagepare="WT_D2", inPath=CE_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  D2_D4_CE_diff = read_diff(stagepare="D2_D4", inPath=CE_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D4_CE_diff = read_diff(stagepare="WT_D4", inPath=CE_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D2_D4_CE_diff <- Reduce(union,list(row.names(WT_D2_CE_diff), 
                                        row.names(D2_D4_CE_diff),
                                        row.names(WT_D4_CE_diff)))
  #------ME
  ME_dir <- file.path(inPath,"ME","all_res")
  WT_D2_ME_diff = read_diff(stagepare="WT_D2", inPath=ME_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  D2_D4_ME_diff = read_diff(stagepare="D2_D4", inPath=ME_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D4_ME_diff = read_diff(stagepare="WT_D4", inPath=ME_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D2_D4_ME_diff <- Reduce(union,list(row.names(WT_D2_ME_diff), 
                                        row.names(D2_D4_ME_diff),
                                        row.names(WT_D4_ME_diff)))
  
  #------SNE
  SNE_dir <- file.path(inPath,"SNE","all_res")
  WT_D2_SNE_diff = read_diff(stagepare="WT_D2", inPath=SNE_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  D2_D4_SNE_diff = read_diff(stagepare="D2_D4", inPath=SNE_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D4_SNE_diff = read_diff(stagepare="WT_D4", inPath=SNE_dir, return="diff", 
                            diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D2_D4_SNE_diff <- Reduce(union,list(row.names(WT_D2_SNE_diff), 
                                         row.names(D2_D4_SNE_diff),
                                         row.names(WT_D4_SNE_diff)))
  
  #------cbNE
  cbNE_dir <- file.path(inPath,"cbNE","all_res")
  WT_D2_cbNE_diff = read_diff(stagepare="WT_D2", inPath=cbNE_dir, return="diff", 
                              diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  D2_D4_cbNE_diff = read_diff(stagepare="D2_D4", inPath=cbNE_dir, return="diff", 
                              diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D4_cbNE_diff = read_diff(stagepare="WT_D4", inPath=cbNE_dir, return="diff", 
                              diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D2_D4_cbNE_diff <- Reduce(union,list(row.names(WT_D2_cbNE_diff), 
                                          row.names(D2_D4_cbNE_diff),
                                          row.names(WT_D4_cbNE_diff)))
  #------PE
  PE_dir <- file.path(inPath,"PE","all_res")
  WT_D2_PE_diff = read_diff(stagepare="WT_D2", inPath=PE_dir, return="diff", 
                              diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  D2_D4_PE_diff = read_diff(stagepare="D2_D4", inPath=PE_dir, return="diff", 
                              diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D4_PE_diff = read_diff(stagepare="WT_D4", inPath=PE_dir, return="diff", 
                              diff_fdr=diff_fdr, diff_logfc=diff_logfc, not_diff_fdr=not_diff_fdr)
  WT_D2_D4_PE_diff <- Reduce(union,list(row.names(WT_D2_PE_diff), 
                                        row.names(D2_D4_PE_diff),
                                        row.names(WT_D4_PE_diff)))
  
  locs_diff_gene = Reduce(union,list(WT_D2_D4_CE_diff,WT_D2_D4_ME_diff, WT_D2_D4_SNE_diff, 
                                     WT_D2_D4_cbNE_diff,WT_D2_D4_PE_diff))
  intersect(locs_diff_gene,  WT_D2_D4_input_not_diff)
}

edgeR_dir = "./edgeR_test"

WT_D2_D4_candidate_gene = get_diff_gene(inPath = edgeR_dir, not_diff_fdr = 0.2)

###################
# merge all logFC
##################
locs <- list.files(edgeR_dir)
all_tab <- list()
for(loc in locs){
  inPath <- file.path(edgeR_dir,loc,"all_res")
  all_res <- list.files(inPath)
  stagepare <- stringr::str_replace(string = all_res, pattern = "_all_gene_edgeR_res.txt", replacement = "")
  for(sp in stagepare){
    tab <- data.table::fread(file.path(inPath, paste(sp,"_all_gene_edgeR_res.txt",sep="")),sep="\t")
    tab <- tab[,c("gene_id","logFC")]
    sp_loc <- paste(sp,loc,sep="_")
    tab[,stage:=rep(sp_loc,nrow(tab))]
    all_tab[[sp_loc]] <- tab
  }
}
all_tab_out <- data.table::rbindlist(all_tab)
all_tab_dc <- data.table::dcast.data.table(all_tab_out, gene_id ~ stage, value.var = "logFC")
all_tab_frm <- as.data.frame(all_tab_dc)[,-1]
row.names(all_tab_frm) <- all_tab_dc$gene_id
samples <- grep(colnames(all_tab_frm),pattern = "^WT",value = T)
all_tab_filter <- all_tab_frm[,samples]

#######################################
#heatmap.
#######################################
gene_tpm = fread("../../expression_data/merged_exp/exon_level_TPM.tab",sep="\t",header=T,data.table=F)
row.names(gene_tpm) = gene_tpm$gene_id
gene_exp = gene_tpm[,-c(1:3)]
gene_ratio = fread("../../expression_data/SL_ratio/exon_level_sl_ratio.txt",sep="\t",header=T,data.table=F)
row.names(gene_ratio) = gene_ratio$gene_id
sample2stage = read.table("../../sampleInfor/sample2stage.txt",sep="\t",header=T)
row.names(sample2stage) = sample2stage[,2]
sample2stage$stage_loc = apply(sample2stage,1,function(x){paste(x[3],x[4],sep="_")})

stage_loc_order = c("WT_input","D2_input","D4_input",
                    "WT_CE","D2_CE","D4_CE",
                    "WT_ME","D2_ME","D4_ME",
                    "WT_SNE","D2_SNE","D4_SNE",
                    "WT_cbNE","D2_cbNE","D4_cbNE",
                    "WT_PE","D2_PE","D4_PE")

#------gene TPM>5 at least in one sample in input.
candidate_gene_exp = gene_exp[WT_D2_D4_candidate_gene,]
candidate_gene_median = t(apply(candidate_gene_exp,1,function(x){
  tapply(x,sample2stage[names(x),"stage_loc"],median)}))

input_gene_median <- candidate_gene_median[,c("WT_input","D2_input","D4_input")]
input_gene_median_high <- input_gene_median[apply(input_gene_median,1,function(x){length(x[which(x>=5)])>0}),]



candidate_gene_logFC <- all_tab_filter[row.names(input_gene_median_high),]
colnames(candidate_gene_logFC) <- stringr::str_replace(colnames(candidate_gene_logFC),pattern = "^WT_",replacement = "")
candidate_gene_logFC[is.na(candidate_gene_logFC)] <- 0
WT_logFC <- matrix(0,nrow = nrow(candidate_gene_logFC),ncol=6, 
                   dimnames = list(row.names(candidate_gene_logFC),
                                   paste("WT",c("input","CE","ME","SNE","cbNE","PE"),sep="_")))
WT_logFC <- as.data.frame(WT_logFC)
candidate_gene_logFC <- cbind(candidate_gene_logFC, WT_logFC)

candidate_gene_type <- gene_tpm[row.names(candidate_gene_logFC),1:3]
candidate_gene_pc <- row.names(candidate_gene_type[candidate_gene_type$gene_type=="protein_coding",])
candidate_gene_npc <- setdiff(row.names(candidate_gene_type), candidate_gene_pc)

pc_gene_frm <- candidate_gene_logFC[,stage_loc_order]
#npc_gene_frm <- candidate_gene_logFC[candidate_gene_npc,stage_loc_order]

plot_frm = pc_gene_frm
#plot_frm = npc_gene_frm
plot_frm[plot_frm>=1] =1
plot_frm[plot_frm<= -1] = -1

cluster_num = 5
library(pheatmap)
library(RColorBrewer)
cols <- c("#1565c0","#757575","#c62828")
blue_col <- colorRampPalette(c("#0d47a1","#4fc3f7"))(37.5)
median_col <- colorRampPalette(c("#757575"))(10)
red_col <- colorRampPalette(c("#ef9a9a","#b71c1c"))(52.5)
pal <- c(blue_col, median_col, red_col)

color_breaks <- seq(min(plot_frm),max(plot_frm),length.out = 100)

p = pheatmap(plot_frm[,-c(1:3)],scale="none",cutree_row=cluster_num,cluster_col=F,clustering_method = "ward.D2",
             show_rownames=F,fontsize = 8, gaps_col = seq(3,ncol(plot_frm[,-c(1:3)]),by=3),
             color = pal, breaks = color_breaks)

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
anno_row = data.frame(cluster = cluster_name_ordered)
anno_row$cluster = new_cluster_name[as.character(anno_row$cluster)]
anno_row$cluster = factor(anno_row$cluster,
                          levels = paste("cluster",1:cluster_num,sep=""))

cluster_frm = data.frame(gene_id = row.names(anno_row),
                         gene_tpm[row.names(anno_row),c("gene_name","gene_type")],
                         cluster = anno_row[,"cluster"],
                         candidate_gene_median[row.names(anno_row),colnames(plot_frm)],check.names=F)
out_dir = "./WT_D2_D4_dynamic/heatmap_clusters"
dir.create(out_dir,recursive = T)
#write.table(cluster_frm,file.path(out_dir,"WT_D2_D4_input_not_diff_5_loc_diff_ncRNA_heatmap_clusters.txt"),sep="\t",quote=F,row.names=F)
row_gaps <- cumsum(table(cluster_frm$cluster))
plot_frm_ordered = plot_frm[row.names(cluster_frm),]
p = pheatmap(plot_frm_ordered,scale="none",cluster_row=F,cluster_col=F,
             show_rownames=F,annotation_row = anno_row, treeheight_row=F,fontsize = 8,
             gaps_col = seq(3,ncol(plot_frm),by=3), main="Input_not_diff_other_locs_diff",
             color = pal, gaps_row = row_gaps)

#mass spec data
library(stringr)
msp <- fread("./mass_spectrum/mass_spec_data.txt",sep="\t",header=T,data.table=F) 

gene_name <- str_match(msp$Description, pattern = "GN=(\\w+)")
msp_tab <- msp[,c(1:9)]
msp_tab$Description <- gene_name[,2]
msp_tab <- msp_tab[!is.na(msp_tab$Description),]
colnames(msp_tab) <- c("Confidence","Accession","gene_name","WT_1","D2_1","D4_1","WT_2","D2_2","D4_2")
msp_tab_filter <- msp_tab[,c("gene_name","WT_1","D2_1","D4_1","WT_2","D2_2","D4_2")]
msp_tab_filter <- msp_tab_filter[apply(msp_tab_filter,1,function(x){length(x[is.na(x)])<length(x)}),]
msp_tab_median <- apply(msp_tab_filter[,-1], 2, function(x){tapply(x, msp_tab_filter[,1], function(y){mean(y,na.rm = T)})})
msp_tab_median <- na.omit(msp_tab_median)
stages <- rep(c("WT","D2","D4"),2)
names(stages) <- colnames(msp_tab_median)

msp_tab_merge <- t(apply(msp_tab_median,1,function(x){tapply(x,stages[names(x)],mean)}))
write.table(msp_tab_merge,"./mass_spectrum/mass_spec_merge.txt",sep="\t",quote=F)

candidate_msp <- msp_tab_merge[intersect(cluster_frm$gene_name,row.names(msp_tab_merge)),c("WT","D2","D4")]
candidate_msp_scale <- t(apply(candidate_msp,1,scale))
colnames(candidate_msp_scale) <- c("WT","D2","D4")

plot_frm_new <- plot_frm_ordered
row.names(plot_frm_new) <- cluster_frm[row.names(plot_frm_new),"gene_name"]

others_name <- setdiff(row.names(plot_frm_new), row.names(candidate_msp_scale))
plot_frm_others <- matrix(0,nrow=length(others_name),ncol=3,dimnames = list(others_name,c("WT","D2","D4")))
plot_frm_others <- as.data.frame(plot_frm_others)
plot_frm_all <- rbind(candidate_msp_scale, plot_frm_others)
plot_frm_all <- plot_frm_all[row.names(plot_frm_new),]

cols <- c("#1565c0","white","#c62828")
pal <- colorRampPalette(cols)(100)
pheatmap(plot_frm_all,scale = "none",cluster_cols = F,cluster_row=F,fontsize_row = 6, color = pal, 
         show_rownames = F, gaps_row = row_gaps)

p = pheatmap(plot_frm_ordered,scale="none",cluster_row=F,cluster_col=F,
             show_rownames=F,annotation_row = anno_row, treeheight_row=F,fontsize = 8,
             gaps_col = seq(3,ncol(plot_frm),by=3), main="Input_not_diff_other_locs_diff",
             color = pal, gaps_row = row_gaps)

