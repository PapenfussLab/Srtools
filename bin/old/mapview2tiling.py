#!/usr/bin/env python

"""
mapview2tiling.py

Convert mapview format to a tiling

Author: Tony Papenfuss
Date: Mon Jun 16 15:17:20 EST 2008

"""

import os, sys, copy
from useful import progressMessage
from gff import Feature


if '-h' in sys.argv:
    sys.exit(__doc__)


def loadChrSizes(filename):
    data = {}
    for line in open(filename):
        tokens = line.strip().split('\t')
        data[tokens[0]] = int(tokens[1])
    return data


iFilename = '/Users/papenfuss/databases/platypus/venom/solexa/mapview_filtered.txt'  # sys.argv[1]
oFilename = 'tiling.txt'
tileSize = 35
chrSizeFilename = '/Users/papenfuss/databases/chromSizes/ornAna5.txt'

chrSizes = loadChrSizes(chrSizeFilename)

iFile = open(iFilename)
headers = iFile.readline().strip().split('\t')
oFile = open(oFilename, 'w')
format = "%s\t%i\t%i"

chrom = None
lastChrom = None
countDict = {}
for i,line in enumerate(iFile):
    if (i % 1000)==0:
        progressMessage('# reads %s', i, 28000000)
    tokens = line.strip().split('\t')
    d = dict(zip(headers, tokens))
    chrom = d['chrom']
    if chrom=='MT':
        continue
    elif 'Ultra' in chrom or 'Contig' in chrom:
        pass
    else:
        chrom = 'chr%s' % chrom
    start = int(d['start'])
    
    if chrom!=lastChrom and countDict:
        print chrom
        for _wStart in xrange(1, chrSizes[lastChrom], tileSize):
            counts = countDict.get((lastChrom, _wStart), 0)
            print >> oFile, format % (lastChrom, _wStart, counts)
        oFile.flush()
        countDict = {}
    
    wStart = 1+tileSize*int((start-1)/tileSize)
    key = (chrom, wStart)
    try:
        countDict[key] += 1
    except KeyError:
        countDict[key] = 1
    lastChrom = copy.copy(chrom)

if countDict:
    for _wStart in xrange(1, chrSizes[chrom], 35):
        counts = countDict.get((chrom, _wStart), 0)
        print >> oFile, format % (chrom, _wStart, counts)
oFile.close()
