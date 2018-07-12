#!/usr/bin/env python 

''' Script to generate templates for a new evolution module '''

import sys
import os
import datetime

if len(sys.argv)<2:
    print ''' 
 Usage: AddEvolution.py NAME
    '''
    exit(0)

name = sys.argv[1]

# First check if the name is already used

with open("Evolutions.txt","r+") as f:
    for l in f:
        a = l.split()
        if a[0] == name:
            print "Interface for evolution "+name+" already exists, exit"
            exit(0)

# Not present, add new line to the Evolutions.txt file

with  open("Evolutions.txt","a") as f:
    f.write(name+" "+"lib"+name.lower()+"_xfitter.so\n")

# Create directory structures:

print "Creating directories in evolutions/"+name

os.system("mkdir -p evolutions/"+name+"/include")
os.system("mkdir -p evolutions/"+name+"/src")
os.system("mkdir -p evolutions/"+name+"/yaml")
os.system("touch evolutions/"+name+"/yaml/parameters.yaml")


print "Creating header file  evolutions/"+name+"/include/Evolution"+name+".h"

with open("evolutions/"+name+"/include/Evolution"+name+".h","w+") as f:
    f.write(
'''
#pragma once

/**
  @class' Evolution'''+name+'''

  @brief A wrapper class for '''+name+''' evolution 

  @version 0.1
  @date ''' + datetime.date.today().isoformat() + '''
  */

class Evolution'''+name+'''
{
  public:
    Evolution'''+name+'''(){};

//    ~Evolution'''+name+'''(){};
//    ~Evolution'''+name+'''(const Evolution'''+name+''' &){};
//    Evolution'''+name+''' & operator =(const Evolution'''+name+''' &r){return *(new Evolution'''+name+'''(r));};

  public:
    virtual string getEvolutionName() const { return  "'''+name+ '''" ;};
    int initAtStart(const string &); 
    virtual int compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err);
  protected:
    virtual int parseOptions(){ return 0;};
};

''')

print "Creating source file  evolutions/"+name+"/src/Evolution"+name+".cc"
with open("evolutions/"+name+"/src/Evolution"+name+".cc","w+") as f:
    f.write(''' 
/*
   @file Evolution'''+name+'''.cc
   @date ''' + datetime.date.today().isoformat() + '''
   @author  AddEvolution.py
   Created by  AddEvolution.py on ''' + datetime.date.today().isoformat() + '''
*/

#include "Evolution'''+name+'''.h"

// the class factories
extern "C" Evolution'''+name+'''* create() {
  return new Evolution'''+name+'''();
}


// Initialize at the start of the computation
int Evolution'''+name+'''::initAtStart(const string &s)
{
  return 0;
}

// Main function to compute results at an iteration
int Evolution'''+name+'''::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  return 0;
}

''')


print "Creating autoconf file  evolutions/"+name+"/src/Makefile.am"
with open("evolutions/"+name+"/src/Makefile.am","w+") as f:
    f.write('''
# Created by AddEvolution.py on ''' + datetime.date.today().isoformat() + '''

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = lib'''+ name.lower() + '''_xfitter.la
lib'''+ name.lower()+'''_xfitter_la_SOURCES = Evolution'''+name+'''.cc

# lib'''+ name.lower()+'''_xfitter_la_LDFLAGS = place_if_needed  

datadir = ${prefix}/yaml/evolutions/'''+name+'''
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
 ''')


print "Update configure.ac file"
os.system("sed 's|xfitter-config|xfitter-config\\n		 evolutions/" +name +"/src/Makefile|' configure.ac  >/tmp/configure.ac")
os.system("cp /tmp/configure.ac configure.ac")

print "Update Makefile.am"
os.system("sed 's|tools/process|tools/process evolutions/" +name +"/src|' Makefile.am > /tmp/Makefile.am")
os.system("cp /tmp/Makefile.am Makefile.am")

print "Update doxygen.cfg"
os.system("sed 's|reactions/APPLgrid/include|reactions/APPLgrid/include evolutions/" +name +"/src evolutions/" +name +"/include|' doxygen.cfg > /tmp/doxygen.cfg   ")
os.system("cp /tmp/doxygen.cfg  doxygen.cfg")
