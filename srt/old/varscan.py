"""
srt.varscan
"""

from srt.core import *
import sys


class Varscan(AbstractFeature):
    attributes = ["chrom","position","ref","var","reads1","reads2",
        "varFreq","strands1","strands2","qual1","qual2","pvalue",
        "mapqual1","mapqual2","Reads1Plus","Reads1Minus",
        "Reads2Plus","Reads2Minus"]
    converters = [
        ("position", int),
        ("reads1", int),
        ("reads2", int),
        ("strands1", int),
        ("strands2", int),
        ("qual1", int),
        ("qual2", int),
        ("pvalue", float)
    ]
    
    def __init__(self, *args):
        super(Varscan, self).__init__(*args)
    
    def __repr__(self):
        return '\t'.join(['%%(%s)s' % str(a) for a in self.attributes]) % self.__dict__


def VarscanFile(iFileHandle, **kw):
    return VarscanFileReader(iFileHandle, **kw)


class VarscanFileReader(AbstractDataReader):
    def __init__(self, iFileHandle, skip=1):
        super(VarscanFileReader, self).__init__(iFileHandle)
        self.skip = skip
    
    def _generator(self):
        for i in xrange(self.skip):
            self.iFile.next()
        
        for line in self.iFile:
            tokens = line.strip().split('\t')
            try:
                yield Varscan(tokens)
            except KeyError, e:
                print >> sys.stderr, e
                print line
