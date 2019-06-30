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
    ''' Read steering.txt, produce steering.txt_fk '''
    list = []
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
                    list.append(res+'_fk')
    return list

def fk_all():    
    listAll = glob('datafiles/*/*/*/*/*.dat')
    list = []
    for name in listAll:
        res = Fk_file(name)
        if res != 0:
            list.append(res+'_fk')
    return list

def split_steer(list):
    ''' generate steering files to speed up generation of .fk files '''
    
    for i in range(len(list)):
        file = list[i]
        with open("steering.txt_fk"+str(i),"w+") as w:
            with open("steering.txt","r+") as f:
                for l in f:
                    infiles = 0
                    for l in f:
                        res = 0
                        if l.lower().find('&infiles')>=0:
                            infiles = 1
                        if infiles>0:
                            if l.lower().find('&end')>=0:
                                infiles = 0
                                w.write('''
&InFiles 
  NInputFiles = 1
  InputFileNames = ''' +"'"+file+"'"+'\n&End\n')
                        else:
                            if infiles == 0:
                                w.write(l)

# main
if len(argv)<2:
    print (
''' 
Usage: 
   file.dat      -- convert data file, generate file.dat_fk
   steering.txt  -- convert steering.txt and connected data files, generate steering.txt_fx (plus data files)
   ALL           -- convert all data files following standard datafiles/  data tree, generate _fk files
   steering.txt SPLIT -- prepare steering.txt files to run in parallel on the batch, to speedup fk grids generation
''')
    exit (1)

file = argv[1]
    
if file == "steering.txt":
    print ('steering')
    list = fk_steering()
    if len(argv)>2:
        if argv[2] == "SPLIT":
            split_steer(list)
elif file == "ALL":
    print ('all')
    fk_all()
else:
    Fk_file(file)
