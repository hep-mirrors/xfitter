#!/usr/bin/env python

import sys
import commands
import numpy

# make sure  lhapdf is in path:
line = commands.getstatusoutput('lhapdf-config  --libdir')
sys.path.append(line[1]+"/python2.7/site-packages")
import lhapdf

def status():
    print ("Installed")


def loadPDF(outDir="output", lhapdfDir='hf_pdf'):
    lhapdf.pathsAppend(outDir)
    pdf = lhapdf.mkPDF(lhapdfDir)
    return pdf


def fxq( set, dist, x, q):
    return set.xfxQ(dist,x,q)

fq = numpy.vectorize(fxq, excluded = ['set','dist','x'])
fx = numpy.vectorize(fxq, excluded = ['set','dist','q'])
