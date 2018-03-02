#!/usr/bin/env python

import sys

sys.setdlopenflags(256|1)   # MAGIC !!!

# at the moment xfitter.so file is  located in ./src directory:
sys.path.append("./lib")

import libxfitter_fit as xfitter

# execute using standard steering files:
xfitter.logo()
xfitter.read_steer()
xfitter.init_pars()
xfitter.read_data()
xfitter.init_theory()
xfitter.fit()

