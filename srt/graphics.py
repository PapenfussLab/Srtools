"""
srt.graphics
"""

import rpy2.robjects as robjects


def plot(x,y, filename="R.pdf", **kw):
    r = robjects.r
    r.pdf(file=filename)
    r.plot(x, y, **kw)
    r("dev.off()")


def test():
    import numpy
    plot(range(10), list(numpy.arange(10)**2))
