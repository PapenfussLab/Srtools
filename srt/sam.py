"""
srt.sam module
"""

from srt.core import *


class Sam(AbstractFeature):
    attributes = ["qname","flag","rname","pos","mapq","cigar","mrnm","mpos","isize","seq","qual","tag","vtype","value"]
    converters = [
        ("flag", int),
        ("pos", int),
        ("mpos", int),
        ("isize", int),
    ]
    
    def __init__(self, *args):
        super(Sam, self).__init__(*args)
    
    def __repr__(self):
        return '\t'.join(['%%(%s)s' % str(a) for a in self.attributes]) % self.__dict__


def SamFile(iFileHandle, **kw):
    return SamFileReader(iFileHandle, **kw)


class SamFileReader(AbstractDataReader):
    def __init__(self, iFileHandle):
        super(SamFileReader, self).__init__(iFileHandle)
    
    def readHeader(self):
        self.iFile.seek(0)
        header = []
        for line in self.iFile:
            line = line.strip()
            if line[0]=="@":
                header.append(line)
            elif line:
                break
        return header
    
    def _generator(self):
        for line in self.iFile:
            line = line.strip()
            if line[0]=="@":
                continue
            elif line:
                tokens = line.split('\t')
                yield Sam(tokens)


class SamFileWriter(AbstractDataFile):
    def __init__(self, iFileHandle, **kw):
        """
        @param iFileHandle: Output file or name
        """
        self.iFile = smartopen(iFileHandle, "w")
        self.iFilename = self.iFile.name
    
    def writeHeader(self, header):
        for line in header:
            self.iFile.write(line + "\n")
    
    def write(self, samAln):
        self.iFile.write(str(samAln) + "\n")
    
    def __call__(self, samAln):
        self.write(samAln)
