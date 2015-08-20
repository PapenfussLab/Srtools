#!/usr/bin/env python

"""
filterMapview.py <input mapview file> <output filename>

Author: Tony Papenfuss
Date: Mon Jun 16 14:31:42 EST 2008

"""

import os, sys
from useful import progressMessage


iFilename = sys.argv[1]
oFilename = sys.argv[2]

mQ_cutoff = 40
nSeqs = 280000000

oFile = open(oFilename, 'w')
headers = ['name','chrom','start','strand','mQ','numTied','score','numZeroMismatches']
format = '\t'.join(['%s','%s','%i','%s','%i','%i','%i','%i'])

print >> oFile, '\t'.join(headers)

for i,line in enumerate(open(iFilename)):
    tokens = line.strip().split('\t')
    mQ = int(tokens[7])
    if mQ>=mQ_cutoff:
        name = tokens[0]
        chrom = tokens[1]
        start = int(tokens[2])
        strand = tokens[3]
        numTied = int(tokens[10])
        score = int(tokens[11])
        numZeroMismatches = int(tokens[12])
        print >> oFile, format % (name,chrom,start,strand,mQ,numTied,score,numZeroMismatches)
        if (i % 1000)==0: progressMessage('# maq hits %s', i, nSeqs)
oFile.close()
