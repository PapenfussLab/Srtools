#!/usr/bin/env python

"""
dbCreate.py <tab-delimited file>

Author: Tony Papenfuss
Date: Mon Jun 23 10:30:04 EST 2008

"""

import os, sys
from alchemy import *
from blast import HSP, BlastFile
from useful import progressMessage


homeDir = os.environ['HOME']
dsn = 'sqlite:///alignedReads.db'

case = 2
if case==1:
    tableName = 'Devil454'
    iFilename = os.path.join(homeDir, 
        'databases/devil/transcriptome/reads/blast_results/blastn_tumour_v_opossum.txt')
    n = 497921
elif case==2:
    tableName = 'Devil454_alt'
    iFilename = os.path.join(homeDir, 
        'databases/devil/transcriptome/reads/blast_results/blastn_tumour_v_opossum2.txt')
    n = 472244
elif case==3:
    tableName = 'Platy454'
    iFilename = os.path.join(homeDir, 
        'databases/platypus/venom/454/blat_platy_venom_reads.txt')
    n = 208523


metadata = MetaData()

# Create the table and define the mapping
h = HSP()
hsps_table = createTable(tableName, metadata, h.attributes, h.converters, 
    indexedAttributes=['subjectId', 'sStart', 'sEnd'])
mapper(HSP, hsps_table)

# Start a session & initialize database
session = createSession(dsn, metadata)

if case in [1,2]:
    # Devil 454 reads
    for i,line in enumerate(open(iFilename)):
        tokens = line.strip().split('\t')
        h = HSP(tokens[0:-2])
        h.convertBlockToGenomeCoords()
        session.save(h)
        if (i % 5000)==0:
            progressMessage("# HSPs %s", i, n)
            session.commit()
    progressMessage("# HSPs %s\n", i, n)
    session.commit()
elif case==3:
    # Platypus 454 reads
    for i,h in enumerate(BlastFile(iFilename)):
        h.subjectId = h.subjectId.split('|')[1]
        h.convertBlockToGenomeCoords()
        session.save(h)
        if (i % 5000)==0:
            progressMessage("# HSPs %s", i, n)
            session.commit()
    progressMessage("# HSPs %s\n", i, n)
    session.commit()
