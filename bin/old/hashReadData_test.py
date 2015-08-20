#!/usr/bin/env python

"""
hashReadData.py

Author: Tony Papenfuss
Date: Mon Aug  4 11:28:31 EST 2008

"""

import os, sys
import cPickle
from useful import progressMessage
from alchemy import *
from maq import MaqViewFile


# 1. Load annotation
class Gene(object):
    def __init__(self, _id, geneId, chrom, start, end, strand):
        self.id = _id
        self.geneId = geneId
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
    
    def __repr__(self):
        out = [self.geneId, self.chrom]
        out.append('%i' % self.start)
        out.append('%i' % self.end)
        out.append(self.strand)
        return '\t'.join(out)


def initDatabase(dbFilename, iFilename, clobber=False):
    db = Alchemy()
    gene_table = Table("Gene", db.metadata,
        Column('id', Integer, primary_key=True),
        Column('geneId', Text),
        Column('chrom', Text, index=True),
        Column('start', Integer, index=True),
        Column('end', Integer, index=True),
        Column('strand', String(1))
    )
    mapper(Gene, gene_table)
    
    populateDb = clobber or not os.path.exists(dbFilename)
    session = db.startSession('sqlite:///%s' % dbFilename)
    if populateDb:
        iFile = open(iFilename)
        headers = iFile.readline()
        strandTable = {'1': '+', '-1': '-'}
        for i,line in enumerate(iFile):
            if (i % 1000): progressMessage('# genes %s', i)
            tokens = line.strip().split('\t')
            g = Gene(None, tokens[0], tokens[5], int(tokens[6]), int(tokens[7]), 
                strandTable[tokens[8]])
            session.save(g)
        progressMessage('# genes %s', i)
        session.commit()
    
    return session


iDir = '/Users/papenfuss/databases/platypus/ensembl/Release50'
iFilename = os.path.join(iDir, 'mart_names_locations.txt')
dbFilename = os.path.join(iDir, 'mart_names_locations.sqlite')
# dbFilename = ':memory:'
session = initDatabase(dbFilename, iFilename)


# 2a. Parse read alignment file and has results
# 2b. Write out unannotated reads
flankSize = 1000
maqFilename = '/Users/papenfuss/databases/platypus/venom/solexa/mapview.txt'
data = {}
unannFile = open('test_unannotated.txt', 'w')
multFile = open('test_multiple.txt', 'w')
multFile2 = open('test_multiple_gene.txt', 'w')
for i,m in enumerate(MaqViewFile(maqFilename)):
    if i==10000: break
    if (i % 1000)==0: 
        progressMessage('# maq alns %s', i)
    
    q = session.query(Gene).filter(Gene.chrom==m.chrom) \
        .filter(Gene.start<m.start+flankSize) \
        .filter(Gene.end>m.start+32-flankSize).all()
    
    if len(q)==0:
        print >> unannFile, m
        continue
    elif len(q)>1:
        x = set([r.geneId for r in q])
        if len(x)>1:
            print >> multFile, m
            print >> multFile2, "%s\t%s" % (m.name, ','.join(x))
            continue
    
    try:
        data[q[0].geneId] += 1
    except KeyError:
        data[q[0].geneId] = 1

unannFile.close()
multFile.close()
multFile2.close()

cPickle.dump(data, open('test_counts.dat', 'w'))
