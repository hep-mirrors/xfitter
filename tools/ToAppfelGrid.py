#!/usr/bin/env python

from sys import argv
from glob import glob
from os import path

def Fk_file(name):
    ''' For a given data file generate _fk file '''
    # Check if file esist:
    if not path.exists(name):
        return 0

    # first check if applgrid is present:
    present = 0
    with open(name, "r+") as r:
        for l in r:
            if l.find('applgrid')>0:
                present = 1

    if present:
        print ('Converting '+name)
    else:
        print ('Skiping '+name)
        return 0

    with open(name+'_fk',"w+") as w:
        with open(name, "r+") as r:
            for l in r:
                if l.lower().find('termtype')>=0:
                    l = l.replace('applgrid','apfelgrid')
                l = l.replace('.root','.fk')
                w.write(l)
    return name
def fk_steering():
    with open("steering.txt"+"_fk","w+") as w:
        with open("steering.txt","r+") as f:
            infiles = 0
            indata  = 0
            for l in f:
                res = 0
                if l.lower().find('&infiles')>=0:
                    infiles = 1
                if infiles>0:
                    if l.lower().find('&end')>=0:
                        infiles = 0
                if infiles>0:
#                    print (l,infiles)
                    if (l.lower().find('inputfilenames')>=0) & (l.find('=')>=0):
                        n = l.rstrip().split("=")
                        res = Fk_file(n[1].replace(' ','')[1:-1])
                    else:
                        res = Fk_file(l.replace(' ','').replace(',','').rstrip()[1:-1])
                        
                if res == 0:
                    w.write(l)
                else:
                    w.write(l.replace(res,res+'_fk'))
def fk_all():
    list = glob('datafiles/*/*/*/*/*.dat')
    for name in list:
        Fk_file(name)


# main
if len(argv)<2:
    print (
''' 
Usage: 
   file.dat      -- convert data file, generate file.dat_fk
   steering.txt  -- convert steering.txt and connected data files, generate steering.txt_fx (plus data files)
   ALL           -- convert all data files following standard datafiles/  data tree, generate _fk files
''')
    exit (1)

file = argv[1]
    
if file == "steering.txt":
    print ('steering')
    fk_steering()
elif file == "ALL":
    print ('all')
    fk_all()
else:
    Fk_file(file)
