#!/usr/bin/env python

from pandas import *
from matplotlib.pyplot import *
from numpy import *



def plot_pom(a):
    figure(figsize=(10,10))
    subplot(221)
    plot(a.logQ2/log(10.),a.b0)
    legend()
    xlabel('$\log_{10} Q^2$')
    subplot(222)
    plot(a.logQ2/log(10.),a.b1)
    legend()
    xlabel('$\log Q_{10}^2$')
    subplot(223)
    plot(a.logQ2/log(10.),a.q2a0)
    legend()
    xlabel('$\log Q_{10}^2$')
    subplot(224)
    plot(a.logQ2/log(10.),a.q2a1)
    legend()
    xlabel('$\log Q_{10}^2$')
    #plot(a.logQ2,a.b0)
    #plot(a.logQ2,a.q2a0)
    #plot(a.logQ2,a.q2a1)
    #savefig("q2dep.pdf")


def symband(c,v,var):
    cc = array(c[var])
    l = []
    for i in range(len(v)):
        a = array(v[i][var])
        l.append(a)
    vv = array(l)
    o = vv-cc
    b = o.reshape(len(v)/2,2,len(c))
    esym = sqrt(sum(0.25*(b[:,0,:]-b[:,1,:])**2,axis=0))
    return esym

def asymband(c,v,var):
    cc = array(c[var])
    l = []
    for i in range(len(v)):
        a = array(v[i][var])
        l.append(a)
    vv = array(l)
    o = vv-cc
    b = o.reshape(len(v)/2,2,len(c))

    eup = sqrt( sum ( maximum(b[:,0,:],maximum(b[:,1,:],0.0))**2,axis=0))
    edn = sqrt( sum ( maximum(-b[:,0,:],maximum(-b[:,1,:],0.0))**2,axis=0))
    return edn,eup

def plotsub(x,v,err1,err2):
    fill_between(x,v-err1,v+err2, alpha=0.8)
    plot(x,v,'-',color='r')

def add_cols(a):
    a['q2'] = exp(a['logQ2'])
    a['a0'] = a['q2a0']/a['q2']
    a['a1'] = a['q2a1']/a['q2']
    a['r0'] = a['q2a0']/(a['b0']-a['q2a0'])
    a['r1'] = a['q2a1']/(a['b1']-a['q2a1'])
    return a

# main

a= read_csv("pom/pomeron.csv",sep="\s+")
a = add_cols(a)
d = []
for i in range(42):
    t = read_csv("pom/pom_"+str(i+1)+".csv","\s+")
    t = add_cols(t)
    d.append(t)


vars = ["b0", "b1", "a0", "a1","q2a0","q2a1",'r0','r1']
# vars = ["b0", "b1", "q2a0", "q2a1"]
ylab = ["$b_0$", "$b_1$", "a0", "a1", "$Q^2 a_0$", "$Q^2 a_1$",'r0','r1']

x = array(a.logQ2/log(10.))

figure(figsize=(12,12))

for i in range(len(vars)):
#    esym = symband(a,d,vars[i])
    edn,eup   = asymband(a,d,vars[i])
    subplot(420+i+1)

#    plotsub(x,array(a[vars[i]]),esym,esym)
    plotsub(x,array(a[vars[i]]),edn,eup)
    xlabel("$\log_{10} Q^2$",size=16)

    xx = gca()
    xx.yaxis.set_label_coords(-0.13,0.95)

    ylabel(ylab[i],size=16)

savefig("pom.pdf")



