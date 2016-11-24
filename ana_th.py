#!/bin/env python

import os
import sys

from numpy import *

cent,up,down,std = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

os.system('mkdir -p ' + std + '_2')
os.system('mkdir -p ' + std + '_th')

dirs = [cent,up, down]

preds = zeros((3,7,72))
preds.fill(1)

for i in range(3):
    dir = dirs[i]
    with open(dir+"/fittedresults.txt","r+") as f:
        for s in f:
            vals = s.split()
            if (len(vals) >= 13) and (vals[11] != "mod"):
                bin, th, dset = float(vals[0]), float(vals[6]), int(vals[11])
                ibin = int(bin+0.1) - 1
                iset = dset - 100
                preds[i,iset,ibin] = th


up = (preds[1]-preds[0])/(preds[0]+1.e-10)*100
dn = (preds[2]-preds[0])/(preds[0]+1.e-10)*100

# modify theo files

for i in range(7):
    data = 0
    with open(std+"/theo_0"+str(i+1)+".dat","r") as f:
        with open(std+"_2/theo_0"+str(i+1)+".dat","w+") as w:
            for s in f:
                if data == 0:
                    if s.find("NColumn")>0:
                        t = s.split("=")
                        w.write("  NColumn = "+str(int(t[1])+2)+"\n")
                    elif s.find("ColumnType")>0:
                        t = s.split(",")
                        v = t[2].split("*")
                        ss = s.replace(v[0],str(int(v[0])+2))
                        w.write(ss)
                    elif s.find("ColumnName")>0:
                        ss = s.strip() + ", \"thetaW+\" , \"thetaW-\" \n"
                        w.write(ss)
                    elif s.find("Percent")>0:
                        t = s.split("=")
                        v = t[1].split("*")
                        ss = s.replace(v[0],str(int(v[0])+2))
#                        ss = s.strip() + ", 2*F \n"
                        w.write(ss)
                    elif s.find("&End")>-1:
                        data = 1
                        w.write(s)
                    else:
                        w.write(s)
                else:
                    v = s.split()
                    ibin = int(float(v[0])+0.1)-1

                    upv = up[i,ibin]
                    dnv = dn[i,ibin]

                    ss = s.strip() + " " + str(upv) + " " + str(dnv) + "\n"
                    w.write(ss)

for i in range(7):
    data = 0
    with open(std+"/theo_0"+str(i+1)+".dat","r") as f:
        with open(std+"_th/theo_0"+str(i+1)+".dat","w+") as w:
            for s in f:
                if data == 0:
                    if s.find("NColumn")>0:
                        t = s.split("=")
                        w.write("  NColumn = 10 \n")
                    elif s.find("ColumnType")>0:
                        t = s.split(",")
                        v = t[2].split("*")
                        ss = s.replace(v[0],"2")
                        w.write(ss)
                    elif s.find("ColumnName")>0:
                        ss = "  ColumnName = \"Bin\",\"CosThLow\",\"CosThHigh\",\"RapLow\",\"RapHigh\",\"MassLow\",\"MassHigh\",\"theory\", \"thetaW+\" , \"thetaW-\" \n"
                        w.write(ss)
                    elif s.find("Percent")>0:
                        t = s.split("=")
                        v = t[1].split("*")
                        ss = s.replace(v[0],str(int(v[0])+2))
#                        ss = s.strip() + ", 2*F \n"
                        w.write(ss)
                    elif s.find("&End")>-1:
                        data = 1
                        w.write(s)
                    else:
                        w.write(s)
                else:
                    v = s.split()
                    ibin = int(float(v[0])+0.1)-1

                    upv = up[i,ibin]
                    dnv = dn[i,ibin]

                    ss = v[0] + " "  +  v[1] + " " + v[2] + " "  + v[3] + " " + v[4] +" " +v[5] +" "+ v[6] + " "+ v[7] + " " + str(upv) + " " + str(dnv) + "\n"
                    w.write(ss)



