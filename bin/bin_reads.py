#!/usr/bin/env python

"""
binReads.py <binSize> <bamFilename> <outputFilename>

Author: Tony Papenfuss
Date: Fri 26 Jun 2009 16:17:59 EST

Cumulate all reads into a total count for each bin
"""

import os, sys, math
from optparse import OptionParser
import numpy
import pysam
from srt.useful import progressMessage


def loadChromSizes(iFilename):
    sizes = {}
    for line in open(iFilename):
        tokens = line.strip().split("\t")
        sizes[tokens[0]] = int(tokens[1])
    return sizes


def createBins(L, binSize):
    nBins = int(math.ceil(float(L)/binSize))
    return numpy.zeros(nBins, dtype=numpy.int32)


def binReads(bamFilename, chromSizeFilename, oFilename, binSize=1000, igvFormat=True):
    chromSize = loadChromSizes(chromSizeFilename)
    samfile = pysam.Samfile(bamFilename, "rb")
    oFile = open(oFilename, "w")
    oFile.write("chrom\tstart\tend\tscore\n")
    for chrom in chromSize:
        if "chrUn" in chrom or "_gl" in chrom or "random" in chrom: continue
        L = chromSize[chrom]
        bins = createBins(L, binSize)

        for aln in samfile.fetch(reference=chrom):
            binIndex = int(aln.pos/binSize)
            if not aln.is_unmapped and aln.mapq!=0:
                bins[binIndex] += 1

        nBins = bins.shape[0]
        if igvFormat:
            for binIndex in xrange(nBins):
                start = binIndex*binSize
                end = (binIndex+1)*binSize
                oFile.write("%s\t%i\t%i\t.\t%i\n" % (chrom, start, end, bins[binIndex]))
        else:
            for binIndex in xrange(nBins):
                start = binIndex*binSize
                end = (binIndex+1)*binSize
                oFile.write("%s\t%i\t%i\t%i\n" % (chrom, start, end, bins[binIndex]))
        print chrom, L, sum(bins)
    oFile.close()


if __name__=="__main__":
    usage = "usage: %prog [options] <bam filename> <chrom size filename> <output filename>"
    parser = OptionParser(usage=usage)
    parser.add_option("-w", "--width", dest="winSize", help="Window size (default 1000nt)", 
        type="int", default=1000)
    parser.add_option("-i", "--igv", dest="igvFormat", help="Output in IGV format", 
        action="store_true")
    (options, args) = parser.parse_args()
    if len(args)==0:
        print usage
        sys.exit()
    bamFilename = args[0]
    chromSizeFilename = args[1]
    oFilename = args[2]
    binReads(bamFilename, chromSizeFilename, oFilename, options.winSize, options.igvFormat)
