#!/usr/bin/Rscript


is_utr <- function(x) {
	is_five_prime <- ((x$exon_chrom_start == x$"5_utr_start") & ( x$exon_chrom_end == x$"5_utr_end") )
	is_three_prime <- ((x$exon_chrom_start == x$"3_utr_start") & ( x$exon_chrom_end == x$"3_utr_end") )
	return (is_five_prime | is_three_prime)
}

get_cds <- function(x) {
	if(unique(x$strand)[1] == "-1"){
		# Negative
		# 5' side
		x[!is.na(x$"5_utr_start"),]$exon_chrom_end <- x[!is.na(x$"5_utr_start"),]$"5_utr_start"
		# 3' side
		x[!is.na(x$"3_utr_start"),]$exon_chrom_start <- x[!is.na(x$"3_utr_start"),]$"3_utr_end"
	}else {
		# Potitive
		x[!is.na(x$"5_utr_start"),]$exon_chrom_start <- x[!is.na(x$"5_utr_start"),]$"5_utr_end"
		x[!is.na(x$"3_utr_start"),]$exon_chrom_end <- x[!is.na(x$"3_utr_start"),]$"3_utr_start"
	}
	return(x)
}

library(biomaRt)

genes.df <- read.table("genes.txt", sep="\t",header=T, stringsAsFactors=F)

genes <- genes.df$GENE


mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

bm <- getBM(
attributes = c("chromosome_name", "exon_chrom_start","exon_chrom_end",
 "rank", "ensembl_exon_id",
"ensembl_transcript_id", "transcript_length", "external_gene_name", 
"transcript_biotype", "5_utr_start","5_utr_end","3_utr_start","3_utr_end", "strand"),
filters ="external_gene_name", 
values=genes, mart=mart
)

# If 5' and 3' are same then it is UTR

bm <- bm[bm$transcript_biotype == "protein_coding",]
#bm.seq <- getSequence(id = bm$ensembl_transcript_id, seqType="peptide", mart=mart,type="ensembl_transcript_id")
# Now we select the longest transcript as canonical transcrip

headers <- colnames(bm)
empty.df <- as.data.frame(matrix(0, nrow=0, ncol=length(headers)))
colnames(empty.df) <- headers

for(gene in genes) {
	gene_df <- bm[bm$external_gene_name == gene,]
	gene.prot.seq <- getSequence(id=gene_df$ensembl_transcript_id, seqType="peptide", mart=mart, type="ensembl_transcript_id")
	# Select the longest protein coding sequence
	gene.prot.seq$length <- nchar(gene.prot.seq$peptide)
	#print(head(gene.prot.seq))
	maxl <- gene.prot.seq[gene.prot.seq$length == max(gene.prot.seq$length),]$ensembl_transcript_id
	if(length(maxl) >= 2){
		# If there is a tie
		maxl.tmp <- gene_df[gene_df$ensembl_transcript_id %in% maxl,]
		maxl <- maxl.tmp[which(maxl.tmp$transcript_length == max(maxl.tmp$transcript_length)),]$ensembl_transcript_id
	}
	print(paste0(gene, " ", maxl))
	#tx.select <- gene_df[gene_df$transcript_length == max(gene_df$transcript_length),]
	tx.select = gene_df[gene_df$ensembl_transcript_id %in% maxl,]
	utrs <- is_utr(tx.select)
	utrs[is.na(utrs)] <- FALSE
	tx.select <- tx.select[!utrs,]
	tx.select <- get_cds(tx.select)
	write.table(tx.select, paste0(gene, ".exoninfo.txt"), sep="\t", row.names=F, quote=F)
	empty.df <- rbind(empty.df, tx.select)
}

write.table(empty.df, "whole_genes_in_same_file.txlength.selected.txt", sep="\t",quote=F,row.names=F)
write.table(bm, "whole_genes_in_same_file.original.biomart.result.txt", sep="\t", quote=F, row.names=F)
