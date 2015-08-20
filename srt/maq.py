#!/usr/bin/env python

"""
srt.maq.py

Still has external dependencies and uses the mungoCore.AbstractFeature class

Author: Tony Papenfuss
Date: Tue Jun 24 14:32:36 EST 2008

"""

from srt.core import *


class Maq(AbstractFeature):
    attributes = ['name','chrom','start','strand','insertSize','isPaired',
        'mQ','mQ1','mQa','numMismatches','penalty','numZeroMismatchHits',
        'numOneMismatchHits','L','seq','quality']
    converters = [
        ('start', int),
        ('insertSize', int),
        ('mQ', int),
        ('mQ1', int),
        ('mQa', int),
        ('numMismatches', int),
        ('penalty', int),
        ('numZeroMismatchHits', int),
        ('numOneMismatchHits', int),
        ('L', int)
    ]
    # format = attributesToFormat(attributes)
    
    def __init__(self, *args):
        super(Maq, self).__init__(*args)
    
    def __repr__(self):
        return '\t'.join(['%%(%s)s' % str(a) for a in self.attributes]) % self.__dict__


def MaqViewFile(iFileHandle, **kw):
    return MaqViewReader(iFileHandle, **kw)


class MaqViewReader(AbstractDataReader):
    def __init__(self, iFileHandle):
        super(MaqViewReader, self).__init__(iFileHandle)
    
    def _generator(self):
        for line in self.iFile:
            line = line.strip()
            if line and line[0]!='#':
                tokens = line.split('\t')
                yield Maq(tokens)

