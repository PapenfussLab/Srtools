#!/usr/bin/env python

"""
mapview2profile.py <filtered mapview file> <output file>

Convert mapview format to a profile (suitable for gbrowse)

Author: Tony Papenfuss
Date: Mon Jun 16 15:17:20 EST 2008

"""

import os, sys, copy
from useful import progressMessage
from gff import Feature


if '-h' in sys.argv:
    sys.exit(__doc__)

iFilename = '/Users/papenfuss/databases/platypus/venom/solexa/mapview_filtered.txt'  # sys.argv[1]
oFilename = '/Users/papenfuss/platy/venom/gbrowse/solexa250.gff'  # sys.argv[2]
windowSize = 250

iFile = open(iFilename)
headers = iFile.readline().strip().split('\t')
oFile = open(oFilename, 'w')

chrom = None
lastChrom = None
countDict = {}
for i,line in enumerate(iFile):
    if (i % 1000)==0:
        progressMessage('# reads %s', i, 28000000)
    tokens = line.strip().split('\t')
    d = dict(zip(headers, tokens))
    chrom = "%s" % d['chrom']
    start = int(d['start'])
    # print chrom, start, lastChrom, (chrom!=lastChrom and len(countDict)!=0)
    
    if chrom!=lastChrom and len(countDict)!=0:
        countData = countDict.items()
        countData.sort(key=lambda x: x[0])
        for (_chrom,_wStart),_counts in countData:
            g = Feature(
                reference=_chrom,
                source='solexa250',
                type='tlevel',
                start=_wStart,
                end=_wStart+windowSize-1,
                score=_counts,
                group='Solexa normal'
            )
            print >> oFile, g
        oFile.flush()
        countDict = {}
        lastChrom = copy.copy(chrom)
    
    wStart = 1+windowSize*int((start-1)/windowSize)
    key = (chrom, wStart)
    try:
        countDict[key] += 1
    except KeyError:
        countDict[key] = 1

if len(countDict)!=0:
    countData = countDict.items()
    countData.sort(key=lambda x: x[0])
    for (chrom,wStart),counts in countData:
        g = Feature(
            reference=chrom,
            source='solexa250',
            type='tlevel',
            start=wStart,
            end=wStart+windowSize-1,
            score=counts,
            group='Solexa Venom_250'
        )
        print >> oFile, g
    oFile.flush()
oFile.close()
