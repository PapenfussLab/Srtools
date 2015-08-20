"""
srt.annotation module

Classes to support reading in tab-delimited annotation files,
e.g. biomart, gff, ...

"""

from srt.core import *

import sys, re, warnings
from exceptions import NotImplementedError
from srt.intervals import Interval,Intersector
from srt.useful import smartopen



def loadAnnotationList(filename, columns, converters=None, header=True, separator="\t", **kw):
    """Basic annotation loader; return the annotation as a list
    
    If no converter is passed, an attempt is made to automatically build one.
    Columns named start and end are converted to ints. Column named strand is
    converted from 1/-1 to +/-.
    
    @param filename: File name
    @param columns: Names of columns/fields
    @keyword converters: List of field converting functions in same order as columns (default None)
    @keyword header: Skip first line (default True)
    @keyword separator: Column delimiter (default "\\t")
    
    """
    
    reader = AnnotationReader(filename,columns,converters=converters,header=header,separator=separator,**kw)
    annotation = []
    for gene in reader:
        annotation.append(gene)
    return annotation


def hashListOnAccession(annotation, geneIdAttribName="geneId", 
    startAttribName="start", endAttribName="end"):
    """Hash the annotation list on gene accession and sort annotations for each transcript by length
    
    @param annotation: Annotation list
    @param geneIdAttribName: Name of geneId field (default "geneId")
    @param startAttribName: Name of start field (default "start")
    @param endAttribName: Name of end field (default "end")
    
    """
    
    annotationDict = {}
    for gene in annotation:
        geneId = gene[geneIdAttribName]
        try:
            annotationDict[geneId].append(gene)
        except KeyError:
            annotationDict[geneId] = [gene]
    
    for geneId in annotationDict:
        annotationDict[geneId].sort(key=lambda x: -abs(x[endAttribName]-x[startAttribName]))
    return annotationDict


def loadAnnotationIntersectors(filename, columns, converters=None, 
    header=True, separator="\t", referenceColumn="chrom", startColumn="start", 
    endColumn="end", strandColumn="strand", pad=0, returnAnnotation=False, **kw):
    """Load annotation into Intersectors
    
    @param filename: File name
    @param columns: Names of columns in annotation file (list)
    @keyword converters: List of field converting functions in same order as columns (default None)
    @keyword header: Skip first line (default True)
    @keyword separator: Column delimiter (default "\\t")
    @keyword referenceColumn: Reference column/attribute name (default "chrom")
    @keyword startColumn: Start column/attribute name (default "start")
    @keyword endColumn: End column/attribute name (default "end")
    @keyword strandColumn: Strand column/attribute name (default "strand")
    @keyword pad: Pad the annotation with this much sequence
    @keyword returnAnnotation: Return annotation; function returns a tuple (default False)
    
    """
    annotation = loadAnnotationList(filename,columns,converters=converters,header=header,separator=separator,**kw)
    geneIntersectors = {}
    for gene in annotation:
        chrom = gene[referenceColumn]
        strand = gene[strandColumn]
        start = gene[startColumn]
        end = gene[endColumn]
        if start>end:
            start,end = end,start
        elif start==end:
            end += 1
        
        key = chrom
        interval = Interval(start-pad, end+pad, value=gene)
        try:
            geneIntersectors[key].add_interval(interval)
        except:
            geneIntersectors[key] = Intersector()
            geneIntersectors[key].add_interval(interval)
    
    if returnAnnotation:
        return geneIntersectors, annotation
    else:
        return geneIntersectors


class Feature(object):
    """Class for representing features on a sequence.
    Standard fields are name,chrom,start,end and fields (which name these).
    More fields can be added.
    
    """
    def __init__(self, columns=None):
        """Constructor"""
        self.name = ""
        self.chrom = ""
        self.start = -1
        self.end = -1
        if not columns:
            self.fields = ["name","chrom","start","end"]
        else:
            self.fields = columns
    
    def __getitem__(self, key):
        return self.__dict__[key]
    
    def __repr__(self):
        output = []
        for field in self.fields:
            output.append(str(self.__dict__[field]))
        return '\t'.join(output)
    
    def asInterval(self):
        """Represent the feature as an Interval object"""
        return Interval(self.start, self.end, value=self)


def AnnotationFile(iFileHandle, columns, converters=None, header=True, separator="\t", **kw):
    """Generic interface to annotation readers and writers.
    Additional keywords are passed to Reader or Writer object.
    
    @param iFileHandle: File name or object
    @param columns: Names of columns/fields
    @keyword converters: List of field converting functions in same order as columns (default None)
    @keyword header: Skip first line (default True)
    @keyword separator: Column delimiter (default "\\t")
    
    """
    return AnnotationReader(iFileHandle, columns, converters=converters, header=header, separator=separator, **kw)


class AnnotationReader(AbstractDataReader):
    """
    Class supporting the reading of annotation files.
    """
    
    def __init__(self, iFileHandle, columns, converters=None, header=True, separator='\t', **kw):
        """
        Constructor
        
        @param iFileHandle: File name or object
        @keyword columns: List of column names (default None)
        @keyword converters: List of types or conversion functions (default None)
        @keyword header: Skip header line (default True)
        @keyword separator: Column separator/delimiter (default '\\t')
        
        """
        self.iFile = smartopen(iFileHandle)
        self.columns = columns
        if converters==None:
            converters = []
            for name in columns:
                if "start" in name or "end" in name:
                    converters.append(int)
                elif "strand" in name:
                    converters.append(strandConverter)
                else:
                    converters.append(None)
        elif type(converters)==list and type(converters[0])==tuple: # list of tuples
            converters = dict(converters)
            converters = [converters.get(column) for column in columns]
        self.converters = converters
        
        self.separator = separator
        if header:
            headerline = self.iFile.readline()
        
        if 'fields' in kw:
            warnings.warn("Deprecated keyword fields. Use columns",
                category=DeprecationWarning)
            self.columns = kw['fields']
        
        if 'sep' in kw:
            warnings.warn("Deprecated keyword fields. Use columns",
                category=DeprecationWarning)
            self.separator = kw['sep']
    
    def initFromDict(self, d):
        """Set column names and converters using a dictionary"""
        self.columns,self.converters = d.keys(),d.values()
    
    def _generator(self):
        data = []
        for line in self.iFile:
            tokens = line.strip().split(self.separator)
            f = Feature(self.columns)
            for field,converter,token in zip(self.columns, self.converters, tokens):
                # print "_generator", field,converter,token
                try:
                    f.__dict__[field] = converter(token)
                except:
                    f.__dict__[field] = token
            yield f


def DelimitedDataFile(iFileHandle, mode="r", **kw):
    """Factory function for Delimited Data Reader and Writer classes
    
    @param iFileHandle: file name or object
    @keyword mode: file mode (r,w,a)
    
    """
    if "r" in mode:
        return DelimitedDataReader(iFileHandle, **kw)
    elif "w" in mode or "a" in mode:
        return DelimitedDataWriter(iFileHandle, mode=mode, **kw)


class DelimitedDataReader(object):
    def __init__(self, iFileHandle, header=False, separator="\t"):
        """
        Constructor
        
        @param iFileHandle: File name or object
        @keyword header: Skip header line (default False)
        @keyword separator: Column separator/delimiter (default '\\t')
        
        """
        self.iFile = smartopen(iFileHandle)
        self.separator = separator
        if header:
            headerline = self.iFile.readline()
    
    def __iter__(self):
        self._iter = self._generator()
        return self
    
    def next(self):
        for x in self._iter:
            return x
        raise StopIteration
    
    def _generator(self):
        for line in self.iFile:
            yield line.strip().split(self.separator)


class DelimitedDataWriter(AbstractDataFile):
    """Simple class for writing delimited files"""
    
    def __init__(self, fileHandle, header=None, separator="\t", mode='w', **kw):
        """
        @param iFileHandle: Output file or name
        @keyword mode: File mode - write(w) or append(a)
        """
        assert mode in ('w', 'a')
        self.iFile = smartopen(fileHandle, mode)
        self.iFilename = self.iFile.name
        self.separator = separator
        if header:
            self.iFile.write(header + "\n")
    
    def write(self, row):
        """Write a row
        
        @param row: list of row entries
        """
        self.iFile.write(self.separator.join([str(_) for _ in row]) + "\n")
    
    def __call__(self, row):
        """Call interface to write.
        
        @param row: list of row entries
        """
        self.write(row)


def strandConverter(istrand):
    """Convert strand from '1'/'-1' to '+'/'-'. Also handles integers.
    
    @param istrand: strand (integer)
    @return: '+' or '-'
    
    """
    if istrand in ['1', '+1', 1]:
        strand = '+'
    elif istrand in ['-1', -1]:
        strand = '-'
    else:
        raise Exception('Inappropriate strand value %s' % str(strand))
    return strand



def readGeneAnnotationTest():
    import time
    
    filename = "/Users/papenfuss/databases/platypus/annotation/ensembl/Release54/mart_gene_location.txt"
    columns = ["geneId", "chrom", "start", "end", "strand"]
    genes = loadAnnotationList(filename, columns)
    print genes[0:10]
    print
    time.sleep(3)
    
    print hashListOnAccession(genes)
    print
    time.sleep(3)
    
    
    geneIntersectors = loadAnnotationIntersectors(filename, columns)
    print geneIntersectors[("Ultra187", "+")].find(1, 1000000)


if __name__=="__main__":
    readGeneAnnotationTest()
