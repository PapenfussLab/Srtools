#!/usr/bin/env python

"""
bin_reads_into_genes.py

Author: Tony Papenfuss
Date: Tue 23 Feb 2010 21:30:42 EST

"""

import os, sys
from pysam import Samfile


# Input data
annoFilename = "/Users/papenfuss/databases/platypus/annotation/ensembl/Release54/mart_gene_location.txt"
fields = ["geneId", "chrom", "start", "end", "strand"]
converters = [None,None,int,int,strandConverter]

alnFilenames = [
    "/Users/papenfuss/platy/venom/aln/inseason_filtered.txt",
    "/Users/papenfuss/platy/venom/aln/outseason_filtered.txt",
]

oFilename = "counts.txt"


# Construct intersectors
geneIntersectors,annotation = loadAnnotationIntersectors(filename, columns, pad=100, returnAnnotation=True)
counts = {}
for gene in genes:
    count[gene.geneId] = [0,0]


# Parse alignments
# Quick and dirty: just take read 1
for j,alnFilename in enumerate(alnFilenames):
    i = 0
    for line in open(alnFilename):
        i += 1
        if i % 10000==0:
            progressMessage("# lines %s", i)
        
        tokens = line.strip().split("\t")
        if tokens[0][-1]=="2": continue
        strand = tokens[1]
        chrom = tokens[2]
        start = int(tokens[3])
        end = start+len(tokens[4])-1
        
        try:
            overlap = geneIntersectors[(chrom,strand)].find(start, end)
        except KeyError:
            overlap = []
        
        if overlap:
            gene = overlap[0].value
            counts[gene.geneId][j] += 1
    print


# Output counts
oFile = open(oFilename, "w")
for geneId in counts:
    output = (geneId, counts[geneId][0], counts[geneId][1])
    oFile.write("%s\t%i\t%i\n" % output)
oFile.close()
