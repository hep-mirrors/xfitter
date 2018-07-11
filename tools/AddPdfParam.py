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
  @class {:s}PdfParam 

  @brief A class for {:s} pdf parameterisation

  @version 0.1
  @date {:s}
  */

class {:s}PdfParam : public BasePdfParam
{{
  public:
     /// Default constructor. Name is the PDF name
    {:s}PdfParam (const std::string& inName) : BasePdfParam(inName) {{}}
     /// Compute xf(x,pars). Pure virtual method
     virtual double compute( double const x, double const* pars) override final;
    
     /// (optionally) compute moments:
     // virtual double moment( double const* pars, int const iMoment = 1) override final;
}};
'''.format( name, name, datetime.date.today().isoformat(),name,name)
)


sFile = "pdfparams/{:s}PdfParam/src/{:s}PdfParam.cc".format(name,name)

print "Creating source file "+sFile

with open(sFile,"w+") as f:
    f.write(''' 
/*
   @file {:s}PdfParam.cc
   @date {:s}
   @author  AddPdfParam.py
   Created by  AddPdfParam.py on {:s}
*/

#include "{:s}PdfParam.h"


// Main function to compute PDF
double {:s}PdfParam::compute( double const x, double const* pars)
{{
  return 0;
}}

// Optional function to compute integrals:

// double {:s}PdfParam::moment( double const* pars, int const iMoment)
// {{
//   return 0
// }}

'''.format(name,datetime.date.today().isoformat(),datetime.date.today().isoformat(),name,name,name)
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
