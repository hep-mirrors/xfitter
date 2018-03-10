#!/usr/bin/env python

import sys

sys.setdlopenflags(256|1)   # MAGIC !!!

# at the moment xfitter.so file is  located in ./src directory:
sys.path.append("./lib")

import libxfitter_fit as xfitter

chi2l = []

class a(xfitter.Action):
    def Initialize(self):
        print "haha"
    def AtIteration(self):
        chi2 = xfitter.chi2()
        print "hoho ",chi2
        chi2l.append(chi2)

b = a()
xfitter.addAction(b)

xfitter.run()

from matplotlib.pyplot import plot,savefig
plot(chi2l)
savefig("chi2.pdf")

