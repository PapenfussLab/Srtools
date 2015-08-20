"""
srt.SNVs
"""

from srt.core import *


class Varscan:
    def __init__(self, tokens):
        # Chrom   Position        Ref     Var     Reads1  Reads2  VarFreq 
        # Strands1        Strands2        Qual1   Qual2   Pvalue  MapQual1        MapQual2
        # Reads1Plus      Reads1Minus     Reads2Plus      Reads2Minus   VarAllele
        attributes = ["chrom","position","ref","var","nref","nvar","varFreq",
            "ref_strands","var_strands","ref_qual","var_qual","pvalue",
            "ref_mapqual","var_mapqual","ref_plus","ref_minus","var_plus","var_minus","var_allele"]
        converters = [
            ("position", int),
            ("nref", int),
            ("nvar", int),
            ("ref_strands", int),
            ("var_strands", int),
            ("ref_qual", int),
            ("var_qual", int),
            ("pvalue", float),
            ("ref_plus", int),
            ("ref_minus", int),
            ("var_plus", int),
            ("var_minus", int),
        ]
        self.__dict__ = dict(zip(attributes, tokens))
        for key,converter in converters:
            self.__dict__[key] = converter(self.__dict__[key])
        self.line = "\t".join(tokens)
    
    def __repr__(self):
        return self.line


def VarscanFile(iFileHandle, **kw):
    return VarscanFileReader(iFileHandle, **kw)


class VarscanFileReader(AbstractDataReader):
    def __init__(self, iFileHandle, header=True):
        super(VarscanFileReader, self).__init__(iFileHandle)
        self.header = header
    
    def _generator(self):
        if self.header: self.iFile.next()
        for line in self.iFile:
            tokens = line.strip().split('\t')
            try:
                yield Varscan(tokens)
            except KeyError, e:
                print >> sys.stderr, e
                print line


class Snvmix2:
    def __init__(self, tokens):
        self.line = "\t".join(tokens)
        tokens2 = tokens[0].split(":")
        self.chrom = tokens2[0]
        self.position = int(tokens2[1])
        self.ref = tokens[1]
        self.var = tokens[2]
        details = tokens[3].split(",")
        self.nref = int(details[0].split(":")[1])
        self.nvar = int(details[1].split(":")[1])
    
    def __repr__(self):
        return self.line


def Snvmix2File(iFileHandle, **kw):
    return Snvmix2FileReader(iFileHandle, **kw)


class Snvmix2FileReader(AbstractDataReader):
    def __init__(self, iFileHandle, header=False):
        super(Snvmix2FileReader, self).__init__(iFileHandle)
        self.header = header
    
    def _generator(self):
        if self.header: self.iFile.next()
        for line in self.iFile:
            tokens = line.strip().split('\t')
            try:
                yield Snvmix2(tokens)
            except KeyError, e:
                print >> sys.stderr, e
                print line
