#!/usr/bin/env python

''' loops over all steering and minuit cards in input_steering/ directory and runs corresponding tests '''

from glob import glob
from os import system

steerings = glob('input_steering/steering.txt*')
minuits   = glob('input_steering/minuit.in.txt*')
ewpars    = glob('input_steering/ewparam.txt*')

print minuits

# find what exaclty we can test
tests = set()

for st in steerings:
    tests.add(st.replace("input_steering/steering.txt.",""))

for mn in minuits:
    tests.add(mn.replace("input_steering/minuit.in.txt.",""))

for ew in ewpars:
    tests.add(ew.replace("input_steering/ewparam.txt.",""))

# now there is a list of tests:
for t in tests:
    print ("Starting test "+t)
    # reset to defaults:
    system('''cp input_steering/steering.txt.def steering.txt; cp input_steering/minuit.in.txt.def minuit.in.txt; 
            cp input_steering/ewparam.txt.def ewparam.txt''')

    # now maybe find corresponding file:
    for st in steerings:
        if st.find(t)>0:
            print ("copy "+st)
            system('cp '+st+' steering.txt')
    for mn in minuits:
        if mn.find(t)>0:
            print ("copy "+mn)
            system('cp '+mn+' minuit.in.txt')

    for ew in ewpars:
        if ew.find(t)>0:
            print ("copy "+ew)
            system('cp '+ew+' ewparam.txt')

    system('bin/xfitter')
    exit(0)
