"""
srt.fastaQuality module
"""

from srt.core import *

import os.path, sys, re, cPickle, sqlite3, tempfile
import srt.fasta as fasta
from srt.fasta import Interface, IndexMethod
from srt.useful import smartopen

QualityIndexFile = fasta.FastaIndexFile


def QualityFile(iFilename, mode='r', indexed=False, **kw):
    """Factory function for Reader and Writer classes
    
    @param iFilename: Fasta filename
    @keyword mode: read(r), append(a) or write(w)
    """
    if mode=='r' and not indexed:
        return QualityReader(iFilename, **kw)
    elif mode=='r' and indexed:
        return QualityReaderIndexed(iFilename, **kw)
    elif mode in ('w', 'a'):
        return QualityWriter(iFilename, mode=mode, **kw)


class Quality:
    """Fasta quality class"""
    
    def __init__(self, header='', seq=''):
        self.header = header
        self.seq = seq
    
    def __repr__(self):
        out = ">%s\n%s" % (self.header, pretty(self.seq))
        return out


class QualityWriter(fasta.FastaWriter):
    def write(self, header, seq):
        self.iFile.write('>%s\n' % header)
        self.iFile.write(pretty(seq))
        self.iFile.write('\n')


class QualityReader(fasta.FastaReader):
    def _generator(self):
        """Return an iterator to a multi-fasta file."""
        header = ''
        seq = []
        for line in self.iFile:
            line = line.strip()
            if line:
                if line[0]=='>':
                    if seq:
                        yield header, ' '.join(seq)
                    header = line[1:]
                    seq = []
                else:
                    seq.append(line)
        yield header, ' '.join(seq)


class QualityReaderIndexed(fasta.FastaReaderIndexed, QualityReader):
    """Class for accessing fasta files using an index"""
    
    def __init__(self, iFileHandle,  clobber=False, 
        interface=Interface.CONTAINER, method=IndexMethod.SQLITE, **kw):
        """Constructor
        
        @param iFileHandle: Name of input file
        """
        self.iFile = smartopen(iFileHandle)
        self.iFilename = self.iFile.name
        if method==IndexMethod.PICKLE:
            self.indexFile = fasta.FastaIndexPickleFile(self.iFilename)
        elif method==IndexMethod.TEXT:
            self.indexFile = fasta.FastaIndexTextFile(self.iFilename)
        else: # sqlite3 method is the default
            self.indexFile = fasta.FastaIndexFile(self.iFilename)
        
        self.indexFile.build(clobber=clobber)
        self.interface = interface
        self._iter = None
        self._initIter = True
        
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
#===============================================================================
#    def seek(self, i):
#        """Seek position in fasta file
#        
#        @param i: ith fasta entry
#        """
#        pos = self.indexFile.get(i)
#        self.iFile.seek(pos)
#        self._iter = self._generator()
#    
#    def get(self, i):
#        """Return slice of fasta file.
#        
#        @param i: Index or slice
#        """
#        if type(i) in (int,long):
#            self.seek(i)
#            return self.readOne()
#        elif type(i)==slice:
#            assert i.step in (None,1)
#            self.seek(i.start)
#            j = i.start
#            seqs = []
#            for h,s in self._iter:
#                seqs.append((h,s))
#                j += 1
#                if j==i.stop: break
#        return seqs
#    
#    def search(self, accession):
#        """Search for accession - private mapping interface
#        
#        @param accession: Accession key
#        """
#        pos = self.indexFile.search(accession)
#        self.iFile.seek(pos)
#        return self.readOne()
#    
#    def asContainer(self):
#        """Define a container-like interface to the fasta entries."""
#        self.__getitem__ = self.get
#    
#    def asMapping(self, convert=lambda x: x.split()[0]):
#        """Define a mapping interface to the fasta entries.
#        
#        @keyword convert: Anonymous function defining how to convert fasta header into accession
#        """
#        self.__getitem__ = self.search
#    
#    def __len__(self):
#        return len(self.indexFile)
#    
#    def _inMemoryMappingInterface(self, accession):
#        pass
#    
#    def _reprForInMemory(self):
#        pass
#===============================================================================


def pretty(qseq, width=30, joinChar='\n'):
    """
    Return a prettified version of seq. Default returns the string reformatted
    to 60 chars wide, e.g.
    
    pretty(seq, width=10, joinChar=' ') returns the string with a space every 10 chars.
    
    @param qseq: Sequence string
    @param width: Sequence width (default 60)
    @param joinChar: Character to join on (default \\n)
    @rtype: string
    @return: a pretty-version of seq.
    """
    qseq = qseq.split()
    output = []
    for i in xrange(0, len(qseq), width):
        output.append(' '.join(qseq[i:i+width]))
    return joinChar.join(output)
