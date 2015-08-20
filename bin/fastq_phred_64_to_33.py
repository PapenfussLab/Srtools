#!/usr/bin/env python

"""
fastqPhred64to33.py <input filename> <output filename>

Author: Tony Papenfuss
Date: Fri Nov  6 15:44:08 EST 2009

-Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 (although in raw read data the Phred quality score rarely exceeds 60, higher scores are possible in assemblies or read maps).
-Solexa/Illumina 1.0 format can encode a Solexa/Illumina quality score from -5 to 62 using ASCII 59 to 126 (although in raw read data Solexa scores from -5 to 40 only are expected)
-Illumina 1.3+ format can encode a Phred quality score from 0 to 62 using ASCII 64 to 126 (although in raw read data Phred scores from 0 to 40 only are expected).

"""

import os, sys
from srt.fastq import FastqFile
from srt.useful import progressMessage


iFilename = sys.argv[1]
oFilename = sys.argv[2]
if oFilename=="-":
    oFilename = sys.stdout

iFile = FastqFile(iFilename)
oFile = FastqFile(oFilename, "w")

i = 0
for header,seq,qual in iFile:
    i += 1
    if (i % 1000)==0:
        progressMessage("# lines %s", i)
    qual2 = "".join([chr(ord(q)-64+33) for q in qual])
    oFile.write(header, seq, qual2)
progressMessage("# lines %s\n", i)
oFile.close()

