basedir<-"/Users/jengs/Box/QTL/"
rdatadir<-paste0(basedir,"RData/candidate_prioritization")
outputdir<-paste0(basedir,"Output/candidate_prioritization")
annotationdir<-paste0(basedir,"Annotation")

#read in annotation files
#uniprot to reactome annotation
reactome_uniprot<-read.delim(file.path(annotationdir,"UniProt2Reactome_All_Levels.txt"),header=F,stringsAsFactors=F)
colnames(reactome_uniprot)<-c("Uniprot","Reactome ID","Reactome URL","Reactome Pathway","Reactome Evidence","Species")
reactome_mouse_uniprot<-reactome_uniprot[which(reactome_uniprot[,"Species"]=="Mus musculus"),]
tmp<-paste(reactome_mouse_uniprot[,"Uniprot"],reactome_mouse_uniprot[,"Reactome ID"])
reactome_mouse_uniprot_u<-reactome_mouse_uniprot[!duplicated(tmp),]
rm(tmp)

#uniprot
uniprot_annotation<-read.delim(file.path(annotationdir,"MOUSE_10090_idmapping.dat"),stringsAsFactors=F,header=F)
colnames(uniprot_annotation)<-c("Uniprot","ID_type","ID")
uniprot_annotation_gene<-uniprot_annotation[which(uniprot_annotation[,"ID_type"]=="Gene_Name"),]
uniprot_annotation_mgi<-uniprot_annotation[which(uniprot_annotation[,"ID_type"]=="MGI"),]
uniprot_annotation_gene_mgi<-merge(uniprot_annotation_gene[,c("Uniprot","ID")],uniprot_annotation_mgi[,c("Uniprot","ID")],suffix=c("_GeneName","_MGI"),by="Uniprot")

#mgi annotation
mgi_annotation<-read.delim(file.path(annotationdir,"mgi_annotation.rpt"),stringsAsFactors=F)

#read in the candidate prioritization results files
#candidate prioritization results files from candidate_prioritization_rank.R script
rda_files<-list.files(path=rdatadir,pattern="WNV.*rda",full.names=T)
rda_files_8<-grep("8_mb",rda_files,value=T)

#annotateRegion function below
lapply(rda_files_8,annotateRegion,"X")

#set qtl region of interest
qtl_region_interest<-rbind(c("X",140.539364),c("X",146.207907))
colnames(qtl_region_interest)<-c("chromosome","position.B38.")
qtl_region_interest<-as.data.frame(qtl_region_interest,stringsAsFactors=F)
qtl_region_interest[,2]<-as.numeric(qtl_region_interest[,2])

#identifySnpsIndels funciton defined below
identifySnpsIndels(qtl_region_interest_list,"X","WNV_3_mb_bin_mean_composite_mean") 
getCompositeScore<-function(rda_file) {
	tryCatch({load(rda_file); return(get("composite_score_start_list"))},error=function(e) { return(NULL)})

}
getQTLRegionInterest<-function(idx1, idx2_list,chr,composite_score_df) {
	idx2<-idx2_list[[as.character(idx1)]]
	qtl_region_int<-data.frame(rep(chr,length(c(idx1,idx2))),as.numeric(composite_score_df[c(idx1,idx2),1]))
	colnames(qtl_region_int)<-c("chromosome","position.B38.")
	return(qtl_region_int)
}
annotateRegion<-function(rda_file,chr,threshold=NULL) {
        analysis_name<-gsub(paste0(rdatadir,"/"),"",gsub(".rda","",rda_file))
        increment<-unlist(strsplit(analysis_name,"_"))[2]
        print(increment)

	tmp<-getCompositeScore(rda_file)
	print(class(tmp))
	print(is.null(tmp))
	if (!is.null(tmp)) {
	tmp<-tmp[[chr]]
        if (is.null(threshold)) {
                threshold<-max(as.numeric(tmp[,2]))
	}
	idx1<-which(as.numeric(tmp[,2])>=threshold)
	idx2<-lapply(idx1,function(x, df,inc) { which(round(as.numeric(df[,1]))==round(as.numeric(df[x,1]) + as.numeric(inc)))},tmp,increment)
	names(idx2)<-idx1
	qtl_region_interest_list<-lapply(idx1,getQTLRegionInterest,idx2,chr,tmp)
	names(qtl_region_interest_list)<-as.character(idx1)
	chr_min_max_list<-lapply(qtl_region_interest_list,getChrMinMax)

	#idx2<-which(round(as.numeric(tmp[,1]))==round(as.numeric(tmp[idx1,1]) + as.numeric(increment)))
	#qtl_region_interest<-data.frame(rep(chr,length(c(idx1,idx2))),as.numeric(tmp[c(idx1,idx2),1]))
	#colnames(qtl_region_interest)<-c("chromosome","position.B38.")
	
	#identify mgi ids within region of interest
	#region of interest has to lie within gene
	#region_interest_mgi_uniprot_reactome<-annotatedMgiUniprotReactome(qtl_region_interest,chr,analysis_name) 
	region_interest_mgi_uniprot_reactome<-lapply(qtl_region_interest_list,annotatedMgiUniprotReactome,chr,analysis_name)
	names(region_interest_mgi_uniprot_reactome)<-names(qtl_region_interest_list)

	#identify snps/indels
	#missense_indels_trunc<-identifySnpsIndels(qtl_region_interest,chr,analysis_name) 
	missense_indels_trunc<-lapply(qtl_region_interest_list,identifySnpsIndels,chr,analysis_name)
	names(missense_indels_trunc)<-names(qtl_region_interest_list)

	#merge with reactome annotation
#	missense_reactome<-apply(missense_indels_trunc,1,getReactomeAnnotation,region_interest_mgi_uniprot_reactome)
	lapply(names(missense_indels_trunc),mergeSnpsReactome,region_interest_mgi_uniprot_reactome,missense_indels_trunc,chr_min_max_list)
#	missense_reactome<-do.call(rbind,missense_reactome)
#	colnames(missense_reactome)[1:2]<-c("X.CHROM","POS")
#	missense_indels_trunc_reactome<-merge(missense_indels_trunc,missense_reactome,by=c("X.CHROM","POS"),all.x=T)
#	print(file.path(outputdir,paste0(analysis_name,"_chr",chr,"_missense_indels_trunc_reactome.txt")))
#	write.table(missense_indels_trunc_reactome,file=file.path(outputdir,paste0(analysis_name,"_chr",chr,"_missense_indels_trunc_reactome.txt")),sep="\t",col.names=T,row.names=F)
}
}

#min_tmp_idx<-lapply(tmp_idx,min)
#min_tmp_idx<-unlist(min_tmp_idx)
#min_tmp_idx<-min_tmp_idx[is.finite(min_tmp_idx)]

getChrMinMax<-function(qtl_region) {
	return(paste0(qtl_region[1,"chromosome"],"_",min(qtl_region[,"position.B38."]),"_",max(qtl_region[,"position.B38."])))
}
mergeSnpsReactome<-function(idx_txt,reactome,missense,chr_min_max_list) {
	chr_min_max<-chr_min_max_list[[idx_txt]]
	missense_df<-missense[[idx_txt]]
	reactome_df<-reactome[[idx_txt]]
	missense_reactome<-apply(missense_df,1,getReactomeAnnotation,reactome_df)
	missense_reactome<-do.call(rbind,missense_reactome)
        colnames(missense_reactome)[1:2]<-c("X.CHROM","POS")
        missense_indels_trunc_reactome<-merge(missense_df,missense_reactome,by=c("X.CHROM","POS"),all.x=T)
        write.table(missense_indels_trunc_reactome,file=file.path(outputdir,paste0(analysis_name,"_",chr_min_max,"_missense_indels_trunc_reactome.txt")),sep="\t",col.names=T,row.names=F)

}
#function to create tabix files
identifySnpsIndels<-function(qtl_region_interest,chr,analysis_name) {
        high_impact<-c("transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","transcript_amplification")
        moderate_impact<-c("inframe_insertion","inframe_deletion", "missense_variant","regulatory_region_ablation")
        filters<-paste(c(high_impact,moderate_impact),collapse="|")
        chr_min_max<-data.frame(chr,min(qtl_region_interest[,"position.B38."])*1000000,max(qtl_region_interest[,"position.B38."])*1000000)
	chr_min_max_text<-paste0(chr_min_max[1,1],"_",chr_min_max[1,2],"_",chr_min_max[1,3])

        tabix_snps<-createTabixFile("ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz",chr_min_max,filters,"snps")
        write(tabix_snps,file=file.path(outputdir,paste0(analysis_name,"_snps.txt")))
        system(paste0("chmod u+x ",file.path(outputdir,paste0(analysis_name,"_snps.txt"))), intern=T)
        system(file.path(outputdir,paste0(analysis_name,"_snps.txt")),intern=T)

        tabix_indels<-createTabixFile("ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz", chr_min_max,filters,"indels")
        write(tabix_indels,file=file.path(outputdir,paste0(analysis_name,"_indels.txt")))
        system(paste0("chmod u+x ",file.path(outputdir,paste0(analysis_name,"_indels.txt"))), intern=T)
        system(file.path(outputdir,paste0(analysis_name,"_indels.txt")),intern=T)

        tabix_files<-list.files(path=outputdir,pattern=paste0(chr_min_max[1,1],"_",chr_min_max[1,2],"_",chr_min_max[1,3],"_[indels|snps]"),full.names=T)
	print(tabix_files)
        tabix_output<-lapply(tabix_files,function(x) { if (file.info(x)[1,"size"]>0) { read.delim(x,stringsAsFactors=F,header=F) }})
	print(tabix_output)
	if (!is.null(unlist(tabix_output))) {
        tabix_output<-do.call(rbind,tabix_output)

        header<-read.delim(file.path(annotationdir,"header.out"),stringsAsFactors=F,skip=103)
        colnames(tabix_output)<-colnames(header)

        short_genotype<-apply(tabix_output[,10:ncol(tabix_output)],2,setGenotypeShort)
        missense_indels<-data.frame(tabix_output, short_genotype)
        colnames(missense_indels)<-gsub("_","/",colnames(missense_indels))
        f_interest<-c("CAST/EiJ","PWK/PhJ","WSB/EiJ","C57BL/6J","NZO/HlLtJ","129S1/SvImJ","A/J","NOD/ShiLtJ")
        f_interest[which(f_interest=="129S1/SvImJ")]<-"X129S1/SvImJ"
        f_idx<-unlist(lapply(f_interest,function(x,y) { return (grep(paste("^",x,sep=""),colnames(y)))},missense_indels))
        missense_indels_trunc<-missense_indels[,c(1:9,f_idx)]
        write.table(missense_indels_trunc,file=file.path(outputdir,paste0(analysis_name,"_",chr_min_max_text,"_founder_snps_indels.txt")),sep="\t",col.names=T,row.names=F)
	return(missense_indels_trunc)
	} 
}
annotatedMgiUniprotReactome<-function(qtl_region_interest,chr,analysis_name) {
        #identify mgi ids within region of interest
        #region of interest has to lie within gene
        region_interest_mgi<-mgi_annotation[which(mgi_annotation[,"Chr"]==chr & mgi_annotation[,"genome.coordinate.start"]>= (min(qtl_region_interest[,"position.B38."])*1000000) & mgi_annotation[,"genome.coordinate.end"]<= (max(qtl_region_interest[,"position.B38."])*1000000)),]

	region_int_txt<-paste0("_",chr,"_",(min(qtl_region_interest[,"position.B38."])*1000000),"_",(max(qtl_region_interest[,"position.B38."])*1000000))

        #identify uniprots given mgi id for region of interest
        region_interest_mgi_uniprot<-merge(region_interest_mgi,uniprot_annotation_gene_mgi,by.x="MGI.Accession.ID",by.y="ID_MGI",all.x=T)

        #identify reactome pathways for region of interest
        region_interest_mgi_uniprot_reactome<-merge(region_interest_mgi_uniprot,reactome_mouse_uniprot_u,by.x="Uniprot",by.y="Uniprot",all.x=T)

        #save output
        write.table(region_interest_mgi_uniprot_reactome,file=file.path(outputdir,paste0(analysis_name,region_int_txt,"_reactome.txt")),sep="\t",col.names=T,row.names=F)
	return(region_interest_mgi_uniprot_reactome)

}
createTabixFile<-function(ftp_url,chr_min_max,filters,type) {
	return(paste0("tabix -f ",ftp_url," ", chr_min_max[1,1],":",chr_min_max[1,2],"-",chr_min_max[1,3]," | awk '$8~/",filters,"/ { print $0 }' > ", file.path(outputdir,paste0(chr_min_max[1,1],"_",chr_min_max[1,2],"_",chr_min_max[1,3],"_",type,".txt"))))

}

setGenotypeShort<-function(x) {
        y<-rep("Alt",length(x))
        y[grep("0/0",x)]<-"Ref"
        y[grep("\\./\\.",x)]<-"No Call"
        return(y)
}



getReactomeAnnotation<-function(missense_indels_trunc_row,reactome_annotation) {
	chrom<-missense_indels_trunc_row["X.CHROM"]
	pos<-missense_indels_trunc_row["POS"]
	missense_reactome<-reactome_annotation[which(reactome_annotation[,"genome.coordinate.start"]<=missense_indels_trunc_row["POS"] & reactome_annotation[,"genome.coordinate.end"]>=missense_indels_trunc_row["POS"]),]
	return(data.frame(rep(chrom,nrow(missense_reactome)),rep(pos,nrow(missense_reactome)),missense_reactome))
}







