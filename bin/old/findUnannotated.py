#!/usr/bin/env python

"""
findUnannotated.py

Author: Tony Papenfuss
Date: Fri Aug 15 12:19:24 EST 2008

"""

import os, sys
from bx.intervals.intersection import *
from fasta import FastaFile
from blast import BlastFile
from useful import progressMessage


print "Load Solexa contigs & store as Intervals in an Intersecter object"
contigData = {}
for h,seq in FastaFile('../solexa/solexa_contigs.fa'):
    tokens = h.split()
    name = tokens[0]
    chrom,se = tokens[1].split(':')
    start,end = [int(x) for x in se.split('-')]
    
    contig = Interval(start, end, value=(name,chrom,start,end,seq))
    try:
        contigData[chrom].add_interval(contig)
    except KeyError:
        contigData[chrom] = Intersecter()
        contigData[chrom].add_interval(contig)


# print "Load 454 contig HSPs & store"
# for b in BlastFile('../454/blastn_contigs_v_genome.txt'):
#     b.convertBlockToGenomeCoords()
#     contig = Interval(b.sStart, b.sEnd, value=(name,b.subjectId,b.sStart,b.sEnd,''))
#     try:
#         contigData[chrom].add_interval(contig)
#     except KeyError:
#         contigData[chrom] = Intersecter()
#         contigData[chrom].add_interval(contig)



print 'Parse genes'
iFilename = '/Users/papenfuss/databases/platypus/ensembl/Release50/mart_names_locations.txt'
iFile = open(iFilename)
headers = iFile.readline()

annotated = set()
for i,line in enumerate(iFile):
    if (i % 1000)==0:
        progressMessage('# genes %s', i)
    
    tokens = line.strip().split('\t')
    geneId = tokens[0]
    transId = tokens[1]
    name = tokens[3]
    chrom = tokens[5]
    start = int(tokens[6])
    end = int(tokens[7])
    strand = {'1': '+', '-1': '-'}[tokens[8]]
    
    try:
        for contig in contigData[chrom].find(start-500, end+500):
            annotated.add(contig.value[0])
    except:
        pass


print 'Parse toxprot alignments'
iFilename = '../toxprot/tblastn_toxprot_v_genome.txt'
for b in BlastFile(iFilename):
    chrom = b.subjectId.split(':')[0]
    try:
        for contig in contigData[chrom].find(b.sStart, b.sEnd):
            annotated.add(contig.value[0])
    except:
        pass


print "Write out what's left over"
writer = FastaFile('unannotated_contigs.fa', 'w')
for chrom in contigData:
    for contig in contigData[chrom].intervals:
        if not contig.value[0] in annotated and len(contig.value[4])>60:
            name,chrom,start,end,seq = contig.value
            writer('%s %s:%i-%i' % (name,chrom,start,end), seq)
writer.close()

    