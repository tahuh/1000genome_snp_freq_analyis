# 1000genome_snp_freq_analysis
1000 Genome SNP frequency analysis scripts


# 1000 genome Amplification Allele Frequency
From 1000 Genomes Allele Frequency Calculator [http://grch37.ensembl.org/info/docs/tools/allelefrequency/index.html] result for some genes

## Example gene list

- AR ( X:66,764,465-66,950,461 )
- EGFR ( 7:55,086,714-55,324,313 )
- HER2(ERBB2) ( 17:37,844,167-37,886,679 )
- FGFR1 ( 8:38,268,656-38,326,352 )
- FGFR2 ( 10:123,237,848-123,357,972 )
- MET ( 7:116,312,444-116,438,440 )
- KRAS ( 12:25,357,723-25,403,870 )
- MYC ( 8:128,747,680-128,753,674 )

Those are usually amplification target genes in clinical samples

## Extraction method
First use Allele Frequency Calculator using link above to extract all population's allele frequency from 1000 Genome DB.

Then use R script calles extractExonRegion.R to extract exon information for each genes from bioMart. This script will select the longest transcript for a given gene if multiple transcript is supporting that.

This will create 3 files.
1. genename.exoninfo.txt - Only include exon information for given gene
2. whole_genes_in_same_file.original.biomart.result.txt - Integrated files
3. whole_genes_in_same_file.txlength.selected.txt - From file 2, the longest 1) protein sequence, 2)transcript length is chosen

Use Ensembl's liftOver (Assembly Converter[http://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter]) to change genomic coordinate from GRCh38 to GRCh37 since as of 20-01-2020, the 1000 genome result is GRCh37 and bioMart uses GRCh38.

After that use separate_by_liftover_gene.py to split all data

Then run snpAlleleFreqStat.py for stat calculation

Finally run postProcess.R for table build

Note. snpCover.R is a script to analyze how many SNPs are covered in cohort. Please change to proper filename at *line 13* and *line 23*

## Genome build

As of current (2020-01-20), 1000 genome offers this tool in with GRCh37.p13 assembly result

## License
MIT

## Author
Thomas Sunghoon Heo
