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

genes.fobj <- read.table("genes.txt", header=T)
genes <- genes.fobj$GENE

O.DF <- data.frame(gene=c(), total=c(), total_vaf1=c(), total_vaf5=c(),total_vaf10=c(), exon_total=c(), exon_vaf1=c(), exon_vaf5=c(), exon_vaf10=c())
RS.DF <-data.frame(gene=c(), total=c(), total_vaf1=c(), total_vaf5=c(),total_vaf10=c(), exon_total=c(), exon_vaf1=c(), exon_vaf5=c(), exon_vaf10=c()) 

for(gene in genes){
	filename <- paste0(gene, ".snps")
	DF <- read.table(filename, sep="\t",header=T,stringsAsFactors=F)
	
	
	
	n <- ncol(DF)
	DF$VAF1 <- rowSums(DF[,7:n] >= 1)
	DF$VAF5 <- rowSums(DF[,7:n] >= 5)
	DF$VAF10 <- rowSums(DF[,7:n] >= 10)
	nsnps <- nrow(DF)
	vaf1 <- nrow(DF[DF$VAF1 >= 1,])
	vaf5 <- nrow(DF[DF$VAF5 >= 1,])
	vaf10 <- nrow(DF[DF$VAF10 >= 1,])
	EXON.DF <- DF[DF$isExonic == "YES",]
	nsnps.exon <- nrow(EXON.DF)
	vaf1.exon <- nrow(EXON.DF[EXON.DF$VAF1 >= 1,])
	vaf5.exon <- nrow(EXON.DF[EXON.DF$VAF1 >= 5,])
	vaf10.exon <- nrow(EXON.DF[EXON.DF$VAF1 >= 10,])
	new.df <- data.frame(gene=c(gene), total=c(nsnps), total_vaf1=c(vaf1), total_vaf5=c(vaf5), total_vaf10=c(vaf10),
	exon_total=c(nsnps.exon), exon_vaf1=c(vaf1.exon), exon_vaf5=c(vaf5.exon), exon_vaf10=c(vaf10.exon))
	
	O.DF <- rbind(O.DF, new.df)
	
	
	## RS annotated
	RS <- DF[startswith(DF$rsid, "rs"),]
	
	
	nsnps <- nrow(RS)
	vaf1 <- nrow(RS[RS$VAF1 >= 1,])
	vaf5 <- nrow(RS[RS$VAF5 >= 1,])
	vaf10 <- nrow(RS[RS$VAF10 >= 1,])
	EXON.RS <- RS[RS$isExonic == "YES",]
	nsnps.exon <- nrow(EXON.RS)
	vaf1.exon <- nrow(EXON.RS[EXON.RS$VAF1 >= 1,])
	vaf5.exon <- nrow(EXON.RS[EXON.RS$VAF1 >= 5,])
	vaf10.exon <- nrow(EXON.RS[EXON.RS$VAF1 >= 10,])
	new.df <- data.frame(gene=c(gene), total=c(nsnps), total_vaf1=c(vaf1), total_vaf5=c(vaf5), total_vaf10=c(vaf10),
	exon_total=c(nsnps.exon), exon_vaf1=c(vaf1.exon), exon_vaf5=c(vaf5.exon), exon_vaf10=c(vaf10.exon))
	
	RS.DF <- rbind(RS.DF, new.df)
}

write.table(O.DF, "SNP.stats.txt", quote=F, row.names=F, sep="\t")

write.table(RS.DF, "SNP.RefSeq.stats.txt", quote=F, row.names=F, sep="\t")