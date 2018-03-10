#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
sys.setdlopenflags(256|1)   # MAGIC !!!

# at the moment xfitter.so file is  located in ./lib directory:
sys.path.append("./lib")

import libxfitter_fit as xfitter

class a(xfitter.Action):
    def Initialize(self):
        print ' \n\n Initialize python test action \n ' 
        self.chi2l = []
        self.alphaS = []        
    def AtIteration(self):
        self.chi2l.append(xfitter.chi2())
        self.alphaS.append(xfitter.getParamD("alphas"))
    def Finalize(self):
        print ' \n\n Finalize python test action \n ' 
        plt.subplot(121)
        plt.plot(self.chi2l[1:],label='$\\chi^2$')
        plt.legend()
        plt.subplot(122)
        plt.plot(self.alphaS[1:],label='$\\alpha_S$')
        plt.legend()
        plt.savefig("chi2.pdf")
        
b = a()
xfitter.addAction(b)
xfitter.run()



