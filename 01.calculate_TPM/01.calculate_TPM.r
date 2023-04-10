
#------calculate TPM with expression data.
counts_tab = read.table("./01.calculate_TPM/gene_level_counts.tab",sep="\t",header=T,row.names=1)
#------normalize read counts with gene length.
length_normalize = t(apply(counts_tab,1,function(x){x[-1]/x[1]*1000}))
#------calculate TPM.
gene_tpm = apply(length_normalize,2,function(x){1e6*x/sum(x)})

#############################################################################
#FUNCTION:calculate TPM for each gene.
#############################################################################
#------input requriment.
#------row names is gene name,first column is exon length for a gene,
#------other column is counts for each all samples of a gene.

calculate_tpm <- function(count_mat){
length_normalize = t(apply(count_mat,1,function(x){x[-1]/x[1]*1000}))
tpm = apply(length_normalize,2,function(x){1e6*x/sum(x)})
return(tpm)
}

gene_tpm = calculate_tpm(count_mat = counts_tab)

write.table(gene_tpm,"./01.calculate_TPM/gene_level_TPM.tab",sep="\t",quote=F)
