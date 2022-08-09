#!/usr/bin/env python3

import argparse
from glob import glob
from os import system

parser = argparse.ArgumentParser()

parser.add_argument('--dir',default='',help='Directory with unpacked artifacts')

args = parser.parse_args()

if (args.dir==''):
    print ('provide directory name')
    exit(0)
    
ll = glob(args.dir+"/temp/*/test.log")

for f in ll:
    with open(f,"r+") as fo:
        for l in fo:
            if (l.find("... FAILED")>0) and (l.find("numdiff")==0):
                s = l.split(" ")
                command = f"cp {args.dir}/{s[5]} {s[6]}"
                system(command)
