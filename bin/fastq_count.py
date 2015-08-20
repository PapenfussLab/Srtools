#!/usr/bin/env python

"""
fastq_count.py
"""

from srt.fastq import FastqFile
from optparse import OptionParser

parser = OptionParser()
(options, args) = parser.parse_args()

for filename in args:
    i = 0
    for h,s,q in FastqFile(filename):
        i += 1
    print "%s\t%i" % (filename, i)

