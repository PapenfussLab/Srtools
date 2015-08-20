#!/usr/bin/env python

"""
maq_remove_multi_mappers.py <input maqview file>

Author: Tony Papenfuss
Date: Wed 10 Jun 2009 14:48:16 EST

"""

import os, sys
from shorty import maq


iFilename = sys.argv[1]
oFilename = sys.argv[2]

maqFile = maq.MaqViewReader(iFilename)
oFile = open(oFilename, "w")
for aln in maqFile:
    if aln.mQ>0:
        oFile.write(str(aln) + "\n")
oFile.close()
