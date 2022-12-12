#!/usr/bin/env python3

import f90nml
import json
import numpy as np

def read_theo_file(fname):
    '''Read theory prediction'''
    nmlf =  f90nml.read(fname)
    name = nmlf['Data']['name']
    ndata = int(nmlf['Data']['NData'])
    ncol  = int(nmlf['Data']['NColumn'])
    ctype = nmlf['Data']['ColumnType']
    cname = nmlf['Data']['ColumnName']
    perc  = nmlf['Data']['Percent']
    ctype = [c.upper() for c in ctype]
    cname = [c.upper() for c in cname]

    # determine error colums:
    ec = []
    for i,ct in enumerate(ctype):
        if ct == 'ERROR':
            ec.append(i)
    ec = np.array(ec)
    idxTheo = cname.index('THEORY')
    
    vals = np.zeros( (ndata,1+len(ec)) )
    # read from the back:
    irow = ndata
    for line in reversed(list(open(fname))):
        l = line.rstrip()
        if len(l.strip()) ==0:
            continue        
        if (l.find('*') == 0) or (l.upper().find('&END')>=0)  or (l.upper().find('/') ==0):
            break        
        vv = np.array([float(fl) for fl in l.split()])
        errs = vv[ec]
        # Scale to absolute, where needed
        errs = np.where(perc,errs*vv[idxTheo]/100,errs)
        irow -= 1
        vals[irow,1:] = errs
        vals[irow,0] = vv[idxTheo]
    return vals,list(np.array(cname)[ec])

def read_data_file(fname, cuts=None, valsTheo=None):
    '''Read namelist xFitter data file, output a tuple of lists'''
    nmlf =  f90nml.read(fname)
    name = nmlf['Data']['name']
    reaction = nmlf['Data']['reaction']
    ndata = int(nmlf['Data']['NData'])
    ncol  = int(nmlf['Data']['NColumn'])
    ctype = nmlf['Data']['ColumnType']
    cname = nmlf['Data']['ColumnName']
    perc  = nmlf['Data']['Percent']
    ctype = [c.upper() for c in ctype]
    cname = [c.upper() for c in cname]
    
    vals = np.zeros( (ndata,ncol) )

    # read from the end:
    irow = ndata        
    for line in reversed(list(open(fname))):
        l = line.rstrip()
        if len(l.strip()) ==0:
            continue
        if (l.find('*') == 0) or (l.upper().find('&END')>=0)  or (l.upper().find('/') ==0):
            break
        vv = np.array([float(fl) for fl in l.split()])
        irow -= 1
        vals[irow] = vv

    # apply cuts
    vals = apply_cuts(cuts,vals,reaction,cname,ctype)
        
    # transform uncertainties to absolute
    if valsTheo is not None:
        ref = valsTheo
    else:
        idxSig = cname.index('SIGMA')
        ref = vals[:,idxSig]
    ierrC = 0
    for i,cn in enumerate(cname):
        if ctype[i] == 'ERROR':
            if perc[ierrC]:
                vals[:,i] *=  ref / 100.
            ierrC += 1
    return name,reaction,vals,cname,ctype

def read_cuts(nml):
    '''Read cuts from the steering.txt nml, returns cuts dictionary'''
    cuts = dict()
    if 'cuts' in nml:
        processes = nml['cuts']['processname']
        var       = nml['cuts']['variable']
        cutmin    = nml['cuts']['cutvaluemin']
        cutmax    = nml['cuts']['cutvaluemax']
        for p,v,c1,c2 in zip(processes,var,cutmin,cutmax):
            cuts[(p,v.upper())] = (c1,c2)
    return cuts

def apply_cuts(cuts,vals, reaction, cname, ctype):
    '''Apply cuts. Returns vals with removed columns'''
    vout = []
    if reaction in list(cuts.keys())[0]:
        nms = [p[1] for p in list(cuts.keys()) if p[0] == reaction]
        idx = [cname.index(n) for n in nms]
        lims = [cuts[(reaction,n)] for n in nms]
    else:
        idx = None
    if 'FLAG' in ctype:
        idxFlag = ctype['FLAG']
    else:
        idxFlag = None
        
    for i in range(vals.shape[0]):
        add = True
        if idxFlag is not None:
            if vals[i,idxFlag] == 0:
                continue
        if idx is not None:
            for ii,l in zip(idx,lims):
#                print (i,ii,l,vals[i,ii])
                if ( vals[i,ii] < l[0] ) or ( vals[i,ii] >= l[1] ):
                    add = False
        if add:
            vout.append(list(vals[i]))
    return np.array(vout)

def to_pyhf(name,reaction, vals, cname, ctype, valsTheo = None, nameTh = None):
    '''Prepare samples and observables'''

    idxSig = cname.index('SIGMA')
    if valsTheo is not None:
        ref = valsTheo[:,0]
    else:
        ref = vals[:,idxSig]
    sample = dict()
    sample['name'] = 'xfitter'
    sample['data'] = list(ref)
    sample['modifiers'] = []
    statCnt = 0
    uncCnt = 0
    for i,cn in enumerate(cname):
        ct = ctype[i] # the list can be longer
        if ct == 'ERROR':
            md = dict()
            
            if cn == 'STAT':
                md['type'] = 'staterror'
                md['data'] = list(vals[:,i])
                md['name'] = f"{cn}_{statCnt}"
                statCnt += 1
            elif cn == 'UNCOR':
                md['type'] = 'shapesys'
                md['data'] = list(vals[:,i])
                md['name'] = f"{cn}_{uncCnt}"
                uncCnt += 1
            elif cn == 'IGNORE':
                continue
            else:
                md['type'] = 'histosys'
                md['data'] = {
                    "hi_data" : list(ref+vals[:,i]),
                    "lo_data" : list(ref-vals[:,i])
                }
                md['name'] = cn
            sample['modifiers'].append(md)

    # Also add theory errors:
    if nameTh is not None:
        for i,tn in enumerate(nameTh):
           md = dict()
           md['type'] = 'histosys'
           md['data'] = {
               "hi_data" : list(ref+valsTheo[:,1+i]),
               "lo_data" : list(ref-valsTheo[:,1+i])
               }
           md['name'] = tn
           sample['modifiers'].append(md)
           
    # also ovservation part
    observation = dict()
    observation['data'] = list(vals[:,idxSig])
    observation['name'] = name
    return name,sample,observation
