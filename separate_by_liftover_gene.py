#!/usr/bin/python

import sys
infile = open(sys.argv[1])
d = {}
#header = infile.readline()
for line in infile:
	data = line.rstrip().split("\t")
	gene = data[-2]
	try:
		d[gene].append(line)
	except KeyError:
		d[gene] = [line]
#print header
for gene in d.keys():
	outfile = open(gene + ".liftover.converted.txt", "w")
	#outfile.write(header)
	for line in d[gene]:
		outfile.write(line)
	outfile.close()
