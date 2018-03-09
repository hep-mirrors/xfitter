#!/usr/bin/env python

import sys

sys.setdlopenflags(256|1)   # MAGIC !!!

# at the moment xfitter.so file is  located in ./src directory:
sys.path.append("./lib")

import libxfitter_fit as xfitter

# execute using standard steering files:
xfitter.logo()
#xfitter.read_steer()
#xfitter.init_pars()
#xfitter.read_data()
#xfitter.init_theory()
#xfitter.fit()

class a(xfitter.Action):
    def Initialize(self):
        print "haha"
    def AtIteration(self):
        print "hoho"


b = a()


xfitter.addAction(a())


base = xfitter.Base()
class Derived(xfitter.Base):
    def f(self):
        return 42

derived = Derived()

print base.f()
print derived.f()


print xfitter.calls_f(base)
print xfitter.calls_f(derived)

print "Call python"
xfitter.actionIteration(a())
print "Call base c++"
xfitter.actionIteration(xfitter.Action())
