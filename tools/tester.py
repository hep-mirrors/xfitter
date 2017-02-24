#!/usr/bin/env python

''' loops over all steering and minuit cards in input_steering/ directory and runs corresponding tests '''

from glob import glob
from os import system
from sys import argv


if len(argv)<2:
    print (
'''Usage:
to run locally:
  tester.py local
to submit (only for nafhh-herafitter): 
  tester.py Submit  
to analyse (only for nafhh-herafitter):
  tester.py Ana   
'''
    )
    exit(0)

steerings = glob('input_steering/steering.txt*')
minuits   = glob('input_steering/minuit.in.txt*')
ewpars    = glob('input_steering/ewparam.txt*')


# print minuits
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
        if st.replace("input_steering/steering.txt.","") == t:
           # print ("copy "+st)
            system('cp '+st+' steering.txt')
    for mn in minuits:
        if mn.replace("input_steering/minuit.in.txt.","") == t:
           # print ("copy "+mn)
            system('cp '+mn+' minuit.in.txt')
    for ew in ewpars:
        if ew.replace("input_steering/ewparam.txt.","")== t:
           # print ("copy "+ew)
            system('cp '+ew+' ewparam.txt')

    if argv[1] == "Submit":
        system('./run_one.py UT_'+t)
    elif argv[1] == "Ana":
        system("ls -l  batch_out/UT_"+t+'/0/xfitter.log')
        system("grep 'After'  batch_out/UT_"+t+'/0/xfitter.log')
    elif argv[1] == "local":
        print ("Run xfitter...")
        system("bin/xfitter > UT_"+t+".out")
