#!/usr/bin/env python

"""
loadSolexa.py

Author: Tony Papenfuss
Date: Tue Jun 24 14:27:34 EST 2008

"""

import os, sys
from maq import *
from useful import progressMessage


oFilename = 'tmp/PlatySolexa.txt'
if not os.path.exists(oFilename):
    oFile = open(oFilename, 'w')
    
    dataDir = '/Users/papenfuss/databases/platypus/venom/solexa/'
    for i,read in enumerate(MaqViewFile(os.path.join(dataDir, 'mapview.txt'), mQ_cutoff=40)):
        if (i % 1000)==0:
            progressMessage("# maq %s", i, 28395347)
        tokens = str(read).split('\t')
        tokens.append(i)
        print >> oFile, "|".join([str(x) for x in tokens])
    oFile.close()
    progressMessage("# maq %s\n", i, 28395347)


os.system("""sqlite3 alignedReads.db '.import "tmp/PlatySolexa.txt" PlatySolexa'""")
