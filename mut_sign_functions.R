library(MutationalPatterns)
library(viridis)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(targets)

#outputdir<-"/home/groups/mcweeney_lab/jengs/HNSCC/mutational_patterns_output/"
outputdir<-"/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/MKM/HNSCC/mutational_patterns_output/"

load("/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/MKM/HNSCC/PureCN/cl_barcode_pairs.RData")
home_dir<-"/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/MKM/HNSCC/"
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#reads in a text file
#first col is vcf file
#second col is tumour barcode
#thir col is normal barcode
readInFile<-function(x) {
        mapping_file<-read.delim(x,as.is=T,header=F)
        colnames(mapping_file)<-c("vcf_file","tumor","normal")
	print(colnames(mapping_file))
        return(mapping_file)
}

#creates an annotated data frame by merging vcf_file from readInFile
#and samp.bams which is a data.table loaded in cl_barcode_pairs.RData
#samp.bams cols are barcode, sample name, labid, tissue, value, found,sample_id, file, file_type
createAnnotFile<-function(file_df,bam_annot) {
	bam_annot<-as.data.frame(bam_annot)
	print(colnames(file_df))
        mutational_sig_file_sample<-merge(file_df,bam_annot[which(!duplicated(bam_annot[,"barcode"]) & bam_annot[,"tissue"]=="Tumor"),c("barcode","sample")],by.x="tumor",by.y="barcode",all.x=T)
        return(mutational_sig_file_sample)
}

#return genomic ranges based on annotated file from createAnnotFile
getGrl<-function(annot_df) {
	vcf_file<-paste0(home_dir,"geno_mutect2_vcf_files/",basename(annot_df[,"vcf_file"]))
	tum_samp_name<-annot_df[,"sample"]
        return(read_vcfs_as_granges(vcf_file, tum_samp_name, ref_genome,type="all"))
}

#gets SNVs from granges
getSnvGrl<-function(all_grl) {
	return(get_mut_type(all_grl,type="snv"))
}

#get mutational matrix from MutationalPatterns package
getMutMat<-function(grl_list) {
	return(mut_matrix(vcf_list=grl_list,ref_genome=ref_genome))
}

#get known signatures from MutationalPatterns package
getKnownSig<-function(mut_type,source,sig_type,incl_poss_art,genome) {
#	return(get_known_signatures(muttype="sn,source="COSMIC_v3.2",sig_type="reference",incl_poss_artifacts=F,genome="GRCh37"))
	return(get_known_signatures(muttype=mut_type,source=source,sig_type=sig_type,incl_poss_artifacts=incl_poss_art,genome=genome))
}

#return fit_to_signature from Mutational Patterns Package
getFitSig<-function(mut_mat,sign) {
	return(fit_to_signatures(mut_mat,sign))
}

#plots conribution from Mutational Patterns Package
plotContribution<-function(fit_res,view_mode) {
	p<-plot_contribution(fit_res$contribution, coord_flip=F,mode=view_mode) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	return(p)
}

#plots original versus reconstructed
plotOGReconstuct<-function(mut_mat,fit_res) {
	p<-plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed,y_intercept = 0.95)
        return(p)
}
#returns get strict refit
getStrictRefit<-function(mut_mat,sign,md,method) {
	strict_refit<-fit_to_signatures_strict(mut_mat,sign,max_delta=md,method=method)
	return(strict_refit)
}
	
#return fit signatures to bootstrapped
performBootstrap<-function(mut_mat,sign,num_boots,method) {
	return(fit_to_signatures_bootstrapped(mut_mat,sign,n_boots=num_boots,method=method))
}

#plots bootstrapped contribution
plotBootStrap<-function(contribution_boots) {
	p<-plot_bootstrapped_contribution(contribution_boots)
	return(p)
}

#createsa a melted df from the strict refit dataframe
#used for plotting
createMeltDFPropMedian<-function(strict_refit) {
	strict_refit_cont<-data.frame(strict_refit$fit_res$contribution)
	proportion_samples<-rowSums(strict_refit_cont>0)/ncol(strict_refit_cont)
	median_contrib<-apply(strict_refit_cont,1,median)
	o<-gsub("[a-z]","",gsub("SBS","",rownames(strict_refit_cont)))
	o<-as.numeric(o)
	strict_refit_cont$order<-o
	strict_refit_cont$proportion<-proportion_samples
	strict_refit_cont$median<-median_contrib
	strict_refit_cont_melt<-strict_refit_cont %>% rownames_to_column(var = "signature") %>% melt(id.vars=c("signature","order","proportion","median"))
        print(head(strict_refit_cont_melt))
	colnames(strict_refit_cont_melt)<-c("signature","order","proportion","median","patient","contribution")
	strict_refit_cont_melt[,"signature"]<-as.factor(strict_refit_cont_melt[,"signature"])
	strict_refit_cont_melt<-strict_refit_cont_melt[which(strict_refit_cont_melt[,"contribution"]!=0),]
	return(strict_refit_cont_melt)
}

#create a dot plot from the melted dataframe
createDotPlot<-function(melted_df) {
	p<-melted_df %>% mutate(signature = fct_reorder(signature, -order)) %>%
  	ggplot() +
   	geom_point(aes(x=patient,y=signature,color = contribution,size=5,na.rm=T)) +
   	scale_colour_gradientn(colours = terrain.colors(7)) +
   	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text=element_text(size=8)) +
  	guides(size = "none")
	return(p)
}

#creates a dotplot from the median melted dataframe (from createMeltDFPropMedian)
createDotPlotMedian<-function(melted_df) {
	p_a<-melted_df %>% mutate(signature = fct_reorder(signature, -order)) %>%
  	ggplot() +
  	geom_point(aes(x=1,y=signature,color=median,size=proportion,na.rm=T)) +
  	scale_color_viridis(discrete=F,direction=-1) +
   	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text=element_text(size=8)) +
   	xlim(c(1,1)) +
   	theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
	return(p_a)
}








