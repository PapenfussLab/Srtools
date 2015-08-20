#!/usr/bin/env python

"""
readBiases.py <input fastq file> <output data file>

Author: Tony Papenfuss
Date: Sat Apr 17 17:22:21 EST 2010

"""

import os, sys
import numpy
from gx.fastq import FastqFile
from gx.useful import progressMessage


iFilename = sys.argv[1]
oFilename = sys.argv[2]

MAX_READ_LENGTH = 600
nSymbols = 5
link = dict(zip(["A","T","G","C","N", "a","t","g","c","n"], range(nSymbols)+range(nSymbols)))
content = numpy.zeros((MAX_READ_LENGTH, nSymbols))

L = 0
n = 0
for h,s,q in FastqFile(iFilename):
    if (n % 10000)==0: progressMessage("# seqs %s", n)
    n += 1
    L = max(L, len(s))
    for pos,nt in enumerate(list(s)):
        content[pos,link[nt]] += 1
progressMessage("# seqs %s\n", n)


oFile = open(oFilename, "w")
oFile.write("\t".join(["A","T","G","C","N"]) + "\n")
for pos in xrange(L):
    output = []
    for nt in xrange(nSymbols):
        output.append("%0.2f" % (100*content[pos, nt]/float(n)))
    oFile.write("\t".join(output) + "\n")
oFile.close()
