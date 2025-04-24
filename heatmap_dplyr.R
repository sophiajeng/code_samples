library(anndata)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

rda_dir<-"/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/PPG/cluster_specific_analysis/rda"
adata<-read_h5ad("/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/PPG/cluster_specific_analysis/h5ad_files/harmony_cmv_clusters_refined_04042025.h5ad")
#get unique clusters and samples
bc_leiden_batch<-adata$obs %>% rownames_to_column() %>% group_by(leiden,batch) %>% select(rowname,leiden,batch) %>% group_split()


#create a list of clusters
#each element will be a list of batch
#each element of batch will be a list of samples/barcodes for that batch in that cluster

#subset X with just hcmv genes
gtf<-read.delim("/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/PPG/scRNA_20250128/genomes/genes_sorted.gtf",stringsAsFactors=F,header=F)

hcmv_gtf<-gtf[which(gtf[,"V2"]=="Geneious"),]
hcmv_gene<-gsub(";.*","",gsub("gene_id ","",hcmv_gtf[,"V9"]))
hcmv_gtf[,"gene"]<-hcmv_gene
hcmv_gene_unique<-hcmv_gtf[!duplicated(hcmv_gtf[,"gene"]),"gene"]
#what is not in the exprset
setdiff(hcmv_gene_unique,rownames(adata$var))
#[1] "US34A"               "AmpR\\(B-Lactamase)" "UL2"                
#[4] "UL19"                "UL48A"               "UL80.5"             
#[7] "UL147A"             
hcmv_idx<-which(rownames(adata$var) %in% hcmv_gene_unique)

#subset counts
hcmv_X<-adata$X[,hcmv_idx]

#get the mean for each hcmv gene
#by batch and leiden
get_hcmv_mean<-function(bc_df,expr) {
	bcs<-pull(bc_df,rowname)
	leiden<-bc_df %>% select(leiden) %>% distinct()
	batch<-bc_df %>% select(batch) %>% distinct()
	if (length(bcs) > 1) {
         	cluster_batch_mean<-Matrix::colMeans(expr[bcs,])
                return(cbind(leiden,batch,cluster_batch_mean))
        }
}

cluster_matrix_list<-map(bc_leiden_batch,get_hcmv_mean,hcmv_X)
create_wider_tibble<-function(long_t) {
	if (!is.null(long_t)) {
        	wider_t<-long_t %>% rownames_to_column(var="gene") %>% pivot_wider(id_cols=c(batch,leiden),names_from=gene,values_from=cluster_batch_mean)
        	return(wider_t)
	}

}
cluster_row<-map(cluster_matrix_list,create_wider_tibble)
cluster_matrix<-do.call(rbind,cluster_row)
hcmv_mean_cluster_list<-cluster_matrix %>% group_by(leiden) 

#create a heatmap, do not cluster batch/gene for easy comparison
#create consistent color gradient
col_fun = colorRamp2(c(0, 3, 6), c("grey", "blue", "red"))
col_fun(seq(-3, 3))
hcmv_gene_order<-hcmv_gene_unique[which(hcmv_gene_unique %in% rownames(adata$var))]
dir<-"/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/PPG/cluster_specific_analysis/heatmaps_dplyr"

create_hm<-function(cluster_tibble,outputdir,gene_order,col_ramp) {
	leiden<-cluster_tibble %>% select(leiden) %>% distinct() %>% pull(leiden)
	print(as.character(leiden))
	cluster_df<-cluster_tibble %>% column_to_rownames(var="batch") %>% select(!(leiden))
	cluster_matrix<-as.matrix(cluster_df)
	hm<-ComplexHeatmap::Heatmap(cluster_matrix,col=col_ramp,column_order=gene_order,cluster_rows=F,row_names_gp=gpar(fontsize=7),column_names_gp=gpar(fontsize=8),name=paste0("cluster_",leiden))
	png(file.path(outputdir,paste0("cluster_",leiden,"_heatmap.png")),width=1500)
	draw(hm)
	dev.off()
}
hm_list<-map(hcmv_mean_cluster_list,create_hm,dir,hcmv_gene_order,col_fun)
	
#create row annotation, batch and leiden/cluster 
annot_batch<-cluster_matrix %>% pull(batch) %>% as.character()
annot_leiden<-cluster_matrix %>% pull(leiden) %>% as.character()
row_annotation<-rowAnnotation(batch=annot_batch,cluster=annot_leiden)

#create a cluster matrix with all clusters/batch for heatmap
cluster_matrix_hm<-cluster_matrix %>% mutate(batch_cluster=paste0(batch,"_",leiden)) %>% column_to_rownames(var="batch_cluster") %>% select(!c(leiden,batch))
cluster_matrix_hm<-as.matrix(cluster_matrix_hm)
all_single_hm<-Heatmap(cluster_matrix_hm,column_order=hcmv_gene_order,cluster_rows=F,row_names_gp=gpar(fontsize=2),column_names_gp=gpar(fontsize=5),name="Clusters 04022025",right_annotation=row_annotation)

#save the single heatmap
png(file.path(dir,"all_clusters_04042025_dplyr_heatmap.png"),width=1500)
draw(all_single_hm)
dev.off()
save.image(file.path(rda_dir,"heatmap_dplyr.rda"))

