#!/usr/bin/env python
import sys
sys.setdlopenflags(256|1)   # MAGIC !!!

# at the moment xfitter.so file is  located in ./lib directory:
sys.path.append("./lib")

import libxfitter_fit as xfitter

chi2l = []
theo1 = []

class a(xfitter.Action):
    def AtIteration(self):
        chi2l.append(xfitter.chi2())
        theo1.append(xfitter.theo(1))
        
b = a()
xfitter.addAction(b)

xfitter.run()

from matplotlib.pyplot import *
subplot(121)
plot(chi2l[1:])
subplot(122)
plot(theo1[1:])
savefig("chi2.pdf")

