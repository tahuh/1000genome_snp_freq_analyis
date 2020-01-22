#!/usr/bin/python


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, required=True, help="Input file.")
parser.add_argument("-e", "--exon_region", type=str, required=True, help="Ensembl Exon region.")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file.")
parser.add_argument("--only-refseq", default=False, help="Analyze RefSeq SNPs(rs) only.", action="store_true", dest="only_refseq")

args = parser.parse_args()

# As of current (20-01-2020), 1000 genome offers GRCh37.p13 as assembly genome
# and calculates population frequency over listed below

POPULATIONS = [
"IBS", "CDX" ,"STU",
"PEL", "PUR", "GBR",
"CHB", "ACB", "ASW",
"TSI", "MXL", "GWD",
"GIH", "PJL", "MSL",
"BEB", "ESN", "CHS",
"FIN", "CEU", "JPT",
"KHV", "YRI", "CLM",
"LWK", "ITU"
]

SUFFIX1='_TOTAL_CNT'
SUFFIX2='_ALT_CNT'

infile = open(args.infile)
header = infile.readline()
hdr2idx = { h : i for i, h in enumerate(header.rstrip().split("\t")) }

total = {}
exon = {}
exon_region = {}
with open(args.exon_region) as exon_file:
	exon_file.readline()
	for line in exon_file:
		data = line.rstrip().split("\t")
		chrom, start,end=data[0],data[1],data[2]
		try:
			exon_region[chrom]
		except KeyError:
			exon_region[chrom] = {}
		for i in range(int(start), int(end) + 1):
			exon_region[chrom][i] = False

n = 0
nanal = 0
nexon = 0
snp_list = []
for line in infile:
	n += 1
	data = line.rstrip().split("\t")
	chrom,pos,rs,ref,malt=data[0],int(data[1]),data[2],data[3],data[4]
	if args.only_refseq:
		if rs == '.':
			continue
	#total["_".join(data[:5])] = {}
	malts = malt.split(",")
	nanal += 1
	try:
		exon_region[chrom][pos]
		nexon += 1
	except KeyError:
		pass
	for pop in POPULATIONS:
		tot = int(data[hdr2idx[pop + SUFFIX1]])
		alts = [int(x) for x in data[hdr2idx[pop + SUFFIX2]].split(",")]
		for i, alt in enumerate(alts):
			dx = data[:4]
			dx.append(malts[i])
			k = "|".join(dx)
			try:
				total[k]
			except KeyError:
				total[k] = {}
			#nanal += 1
			freq = 100 * float(alt) / tot
			total[k][pop] = freq
			#exon[pop] = None
			try:
				exon_region[chrom][pos]
				try:
					exon[k]
				except KeyError:
					exon[k] = {}
				exon[k][pop] = freq
				#nexon += 1
			except KeyError:
				# not exon
				continue

outfile = open(args.outfile, "w")
outfile.write("#TotalSNPs=%d\n"%(n))
outfile.write("#TotalAnalyzed=%d\n"%(nanal))
outfile.write("#TotalExonSNPs=%d\n"%(nexon))
outfile.write("#Total\tTotalVAF1\tTotalVAF5\tTotalVAF10\tExon\tExonVAF1\tExonVAF5\tExonVAF10\n")
vaf1 = 0
vaf5 = 0
vaf10= 0
exon_vaf1 = 0
exon_vaf5 = 0
exon_vaf10 = 0
for k in total.keys():
	nvaf1 = 0
	nvaf5 = 0 
	nvaf10 = 0
	for pop in POPULATIONS:
		vaf = total[k][pop]
		if (vaf >= 1) and (vaf >= 5) and (vaf >= 10):
			nvaf10 += 1
			nvaf1 += 1
			nvaf5 += 1
			continue
		elif (vaf >= 1) and (vaf >= 5) and (vaf < 10):
			nvaf5 += 1
			nvaf1 += 1
			continue
		elif (vaf >= 1) and (vaf < 5) and (vaf < 10):
			nvaf1 += 1
	counted = False
	if nvaf10 >= 1:
		vaf10 += 1
		vaf5 += 1
		vaf1 += 1
		counted = True
	if nvaf5 >= 1 and (not counted):
		vaf5 += 1
		vaf1 += 1
		counted = True
	if nvaf1 >= 1 and (not counted):
		vaf1 += 1
	nvaf1 = 0
	nvaf5 = 0
	nvaf10 = 0
	try:
		pops =exon[k]
		for pop in POPULATIONS:
			vaf = exon[k][pop]
			if (vaf >= 1) and (vaf >= 5) and (vaf >= 10):
				nvaf10 += 1
				nvaf1 += 1
				nvaf5 += 1
				continue
			elif (vaf >= 1) and (vaf >= 5) and (vaf < 10):
				nvaf5 == 1
				nvaf1 += 1
				continue
			elif (vaf >= 1) and (vaf < 5) and (vaf < 10):
				nvaf1 += 1
				continue
		counted = False
		if nvaf10 >= 1:
			exon_vaf10 += 1
			exon_vaf5 += 1
			exon_vaf1 += 1
			counted = True
		if nvaf5 >= 1 and (not counted):
			exon_vaf5 += 1
			exon_vaf1 += 1
			counted = True
		if nvaf1 >= 1 and (not counted):
			exon_vaf1 += 1
			counted = True
	except KeyError:
		continue

outfile.write("#"+ str(nanal) + '\t' + str(vaf1) + '\t' + str(vaf5) + '\t' + str(vaf10) + '\t' + str(nexon) + '\t' + str(exon_vaf1) + '\t' + str(exon_vaf5) + '\t' + str(exon_vaf10) + '\n')
outfile.write("#AnalyzedSnpFreqTable\n")
outfile.write("chrom\tpos\trsid\tref\talt\tisExonic")
for pop in POPULATIONS:
	outfile.write("\t" + pop)
outfile.write("\n")
print args.infile
for k in total.keys():
	chrom,pos,rs,ref,alt=k.split('|')
	outfile.write(str(chrom) + '\t' + str(pos) + '\t' + str(rs) + '\t' + str(ref) + '\t' + str(alt))
	try:
		exon_region[chrom][int(pos)]
		is_exonic = 'YES'
	except KeyError:
		is_exonic = 'NO'
	outfile.write('\t' + is_exonic)
	for pop in POPULATIONS:
		frq = total[k][pop]
		outfile.write('\t' + str(frq))
	outfile.write("\n")
	
