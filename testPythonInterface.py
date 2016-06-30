#!/usr/bin/env python

import sys

# at the moment xfitter.so file is  located in ./src directory:
sys.path.append("./src")

import xfitter

# execute using standard steering files:
xfitter.read_steer_py()
xfitter.read_data_py()
xfitter.init_theory_py()
xfitter.fit()

