#!/usr/bin/env python 

''' Script to generate templates for new PDF parameterisation '''

import sys
import os
import datetime

if len(sys.argv)<2:
    print ''' 
 Usage: AddPdfParam.py NAME
    '''
    exit(0)

name = sys.argv[1]

# First check if the name is already used

print "Creating directories in pdfparams/"+name

os.system("mkdir -p pdfparams/"+name+"PdfParam/include")
os.system("mkdir -p pdfparams/"+name+"PdfParam/src")
os.system("mkdir -p pdfparams/"+name+"PdfParam/yaml")
os.system("touch pdfparams/"+name+"PdfParam/yaml/parameters.yaml")

hFile = "pdfparams/{:s}PdfParam/include/{:s}PdfParam.h".format(name,name)

print "Creating header file  "+hFile

with open(hFile,"w+") as f:
    f.write(
        '''
#pragma once

#include "BasePdfParam.h"

/**
  @class {name:s}PdfParam 

  @brief A class for {name:s} pdf parameterisation

  @version 0.1
  @date {date:s}
  */

class {name:s}PdfParam:public BasePdfParam{{
  public:
    {name:s}PdfParam(const std::string&inName):BasePdfParam(inName){{}}
    //Evaluate xf(x) at given x with current parameters
    virtual double operator()(double x)const override final;
    // (Optional) compute moments:
    // virtual double  moment(int nMoment=-1)const override final;
    // (Optional) set moments:
    // virtual void setMoment(int nMoment,double value)override final;
    // (Optional)
    //Initialize from a yaml node. Uses node[getName] as the basis
    // virtual void initFromYaml(YAML::Node value)override final;
}};
'''.format(name=name,date=datetime.date.today().isoformat())
)


sFile = "pdfparams/{:s}PdfParam/src/{:s}PdfParam.cc".format(name,name)

print "Creating source file "+sFile

with open(sFile,"w+") as f:
    f.write(''' 
/*
   @file {name:s}PdfParam.cc
   @date {date:s}
   @author AddPdfParam.py
   Created by AddPdfParam.py on {date:s}
*/

#include "{name:s}PdfParam.h"

double {name:s}PdfParam::operator()(double x){{
  //Your code here
}}
'''.format(name=name,date=datetime.date.today().isoformat())
)

aFile = "pdfparams/{:s}PdfParam/src/Makefile.am".format(name)

print "Creating autoconf file " + aFile


with open(aFile,"w+") as f:
    f.write('''
# Created by AddPdfParam.py on {:s}

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../BasePdfParam/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = lib{:s}PdfParam_xfitter.la
lib{:s}PdfParam_xfitter_la_SOURCES = {:s}PdfParam.cc

datadir = ${{prefix}}/yaml/pdfparams/{:s}
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
'''.format(datetime.date.today().isoformat(),name,name,name,name))


print "Update configure.ac file"
os.system("sed 's|xfitter-config|xfitter-config\\n		 pdfparams/{:s}PdfParam/src/Makefile|' configure.ac  >/tmp/configure.ac".format(name))
os.system("cp /tmp/configure.ac configure.ac")

print "Update Makefile.am"
os.system("sed 's|pdfdecompositions/BasePdfDecomposition/src|pdfdecompositions/BasePdfDecomposition/src pdfparams/{:s}PdfParam/src|' Makefile.am > /tmp/Makefile.am".format(name))
os.system("cp /tmp/Makefile.am Makefile.am")

print "Update doxygen.cfg"
os.system("sed 's|pdfparams/BasePdfParam/include|pdfparams/BasePdfParam/include  pdfparams/{:s}PdfParam/include|' doxygen.cfg > /tmp/doxygen.cfg".format(name))
os.system("cp /tmp/doxygen.cfg  doxygen.cfg")

exit(0)
