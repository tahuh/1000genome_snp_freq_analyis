#!/usr/bin/Rscript

startswith <- function(string1,string2) {
	v = c()
	substr_size <- nchar(string2)
	for(s in string1){
		sub.s <- substr(s, 1, substr_size)
		v <- c(identical(sub.s, string2), v)
	}
	return (v)
}

genes.fobj <- read.table(proper.file.name1, header=T) # file format should be one column of genes with header name GENE fixed
genes <- genes.fobj$GENE

#RS.DF <-data.frame(gene=c(), total=c(), total_vaf1=c(), total_vaf5=c(),total_vaf10=c(), exon_total=c(), exon_vaf1=c(), exon_vaf5=c(), exon_vaf10=c()) 

MDF <- data.frame(gene=c(), num_patient=c(), total_patient_vars = c() , total_vaf1=c(), total_vaf5=c(),total_vaf10=c(), exon_vaf1=c(), exon_vaf5=c(), exon_vaf10=c()) 

# below input file must be TSV file with header names fixed as
# chrom	pos	rsid	ref	alt	gene	pid
# The pid column is patient id in numeric or identifier for each patient
patient.df <- read.table(proper.cohort.germline.info.file.name , sep="\t", header=T, stringsAsFactors=F)



for(gene in genes){
	filename <- paste0(gene, ".snps")
	DF <- read.table(filename, sep="\t",header=T,stringsAsFactors=F)
	
	n <- ncol(DF)	## RS annotated
	DF$VAF1 <- rowSums(DF[,7:n] >= 1)
	DF$VAF5 <- rowSums(DF[,7:n] >= 5)
	DF$VAF10 <- rowSums(DF[,7:n] >= 10)
	
	RS <- DF[startswith(DF$rsid, "rs"),]
	
	patient <- patient.df[patient.df$gene == gene,]
	npatient <- length(as.vector(unique(patient$pid)))
	
	RS <- merge(RS, patient, by="rsid")
	
	nsnps <- nrow(RS)
	vaf1 <- nrow(RS[RS$VAF1 >= 1,])
	vaf5 <- nrow(RS[RS$VAF5 >= 1,])
	vaf10 <- nrow(RS[RS$VAF10 >= 1,])
	EXON.RS <- RS[RS$isExonic == "YES",]
	nsnps.exon <- nrow(EXON.RS)
	vaf1.exon <- nrow(EXON.RS[EXON.RS$VAF1 >= 1,])
	vaf5.exon <- nrow(EXON.RS[EXON.RS$VAF1 >= 5,])
	vaf10.exon <- nrow(EXON.RS[EXON.RS$VAF1 >= 10,])
	new.df <- data.frame(gene=c(gene), num_patient=c(npatient), total_patient_vars=c(nrow(patient)), 
	total_vaf1=c(vaf1), total_vaf5=c(vaf5), total_vaf10=c(vaf10),
	exon_total=c(nsnps.exon), exon_vaf1=c(vaf1.exon), exon_vaf5=c(vaf5.exon), exon_vaf10=c(vaf10.exon))
	
	MDF <- rbind(MDF, new.df)
}


write.table(MDF, "ALColon_30pts_cover.txt", sep="\t", quote=F, row.names=F)
