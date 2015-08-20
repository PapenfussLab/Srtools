"""
srt.useful module

Module of really useful things
"""

import sys


def smartopen(fileHandle, mode='r'):
    """Return an open file object, if passed a filename;
    otherwise just return the object.
    
    @param fileHandle: Input or output file or filename
    @param mode: IO mode ('r', 'w', ...)
    @return: the open file object
    """
    if type(fileHandle)==str:
        ioFile = open(fileHandle, mode)
        return ioFile
    else:
        return fileHandle


# def progressMessage(format, i, n=None, w=80):
#     if n is None:
#         s = ' %i' % i
#     else:
#         f = 100*float(i)/n
#         s = ' %i/%i (%0.1f%%)' % (i,n,f)
#     message = format % s
#     sys.stderr.write("\b"*w + message)
#     sys.stderr.flush()


def progressMessage(template, *args, **kw):
    """Output a progress message to stderr
    After the first (template) parameter, the args or keywords are treated as
    arguments for template.
    
    @param template:
    """
    w = kw.get('w', 80)
    if len(args)>0:
        message = template % args
    else:
        message = template % kw
    sys.stderr.write("\b"*w + message)
    sys.stderr.flush()
