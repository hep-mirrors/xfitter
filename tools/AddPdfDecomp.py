#!/usr/bin/env python

''' Script to generate templates for new PDF decomposition '''

import sys
import os
import datetime

if len(sys.argv)<2:
    print '''
 Usage: AddPdfDecomp.py NAME
    '''
    exit(0)

name = sys.argv[1]

# First check if the name is already used

with open("Reactions.txt","r+") as f:
    for l in f:
        a = l.split()
        if a[0] == name:
            print "Interface for reaction "+name+" already exists, exit"
            exit(0)

# Not present, add new line to the Reactions.txt file

with  open("Reactions.txt","a") as f:
    f.write(name+" "+"lib"+name.lower()+"_xfitter.so\n")


print "Creating directories in pdfdecompositions/"+name+"PdfDecomposition"

os.system("mkdir -p pdfdecompositions/"+name+"PdfDecomposition/include")
os.system("mkdir -p pdfdecompositions/"+name+"PdfDecomposition/src")
os.system("mkdir -p pdfdecompositions/"+name+"PdfDecomposition/yaml")
os.system("touch pdfdecompositions/"+name+"PdfDecomposition/yaml/parameters.yaml")

hFile = "pdfdecompositions/{:s}PdfDecomposition/include/{:s}PdfDecomposition.h".format(name,name)

print "Creating header file  "+hFile


with open(hFile,"w+") as f:
    f.write(
        '''
#pragma once

#include "BasePdfDecomposition.h"

/**
  @class {:s}PdfDecomposition

  @brief A class for {:s} pdf decomposition

  @version 0.1
  @date {:s}
  */

namespace xfitter {{

class {:s}PdfDecomposition : public BasePdfDecomposition
{{
  public:
     /// Default constructor.
    {:s}PdfDecomposition ();

     /// Default constructor. Name is the PDF name
    {:s}PdfDecomposition (const std::string& inName);

    /// Optional initialization at the first call
    virtual void initAtStart(const std::string & pars) override final;

    /// Compute PDF in a physical base in LHAPDF format for given x and Q
    virtual std::function<std::map<int,double>(const double& x)> f0() const  override final;

}};
}}
'''.format( name, name, datetime.date.today().isoformat(),name,name,name)
)


sFile = "pdfdecompositions/{:s}PdfDecomposition/src/{:s}PdfDecomposition.cc".format(name,name)

print "Creating source file "+sFile

with open(sFile,"w+") as f:
    f.write('''
/*
   @file {:s}PdfDecomposition.cc
   @date {:s}
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on {:s}
*/

#include "{:s}PdfDecomposition.h"

namespace xfitter {{

/// the class factories, for dynamic loading
extern "C" {:s}PdfDecomposition* create() {{
    return new {:s}PdfDecomposition();
}}


// Constructor
    {:s}PdfDecomposition::{:s}PdfDecomposition() : BasePdfDecomposition("{:s}") {{
}}

// Constructor
{:s}PdfDecomposition::{:s}PdfDecomposition(const std::string& inName) : BasePdfDecomposition(inName) {{
}}

// Init at start:
void {:s}PdfDecomposition::initAtStart(const std::string & pars) {{
  return;
}}

// Returns a LHAPDF-style function, that returns PDFs in a physical basis for given x
std::function<std::map<int,double>(const double& x)>  {:s}PdfDecomposition::f0() const
{{
  const auto f_ = [=](double const& x)->std::map<int, double> {{
      std::map<int, double> res_  = {{
	{{-6,0}},
	{{-5,0}},
	{{-4,0}},
	{{-3,0}},
	{{-2,0}},
	{{-1,0}},
	{{ 1,0}},
	{{ 2,0}},
	{{ 3,0}},
	{{ 4,0}},
	{{ 5,0}},
	{{ 6,0}},
        {{22,0}}
      }};
      return res_;
  }};
  return f_;
}}

}}
'''.format(name,datetime.date.today().isoformat(),datetime.date.today().isoformat()
           ,name,name,name,name,name,name,name,name,name,name)
)


aFile = "pdfdecompositions/{:s}PdfDecomposition/src/Makefile.am".format(name)

print "Creating autoconf file " + aFile


with open(aFile,"w+") as f:
    f.write('''
# Created by AddPdfDecomposition.py on {:s}

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../pdfparams/BasePdfParam/include/   -I$(srcdir)/../../BasePdfDecomposition/include -Wall -fPIC -Wno-deprecated

lib_LTLIBRARIES = lib{:s}PdfDecomposition_xfitter.la
lib{:s}PdfDecomposition_xfitter_la_SOURCES = {:s}PdfDecomposition.cc

datadir = ${{prefix}}/yaml/pdfdecompositions/{:s}
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml
'''.format(datetime.date.today().isoformat(),name,name,name,name))



print "Update configure.ac file"
os.system("sed 's|xfitter-config|xfitter-config\\n		 pdfdecompositions/{:s}PdfDecomposition/src/Makefile|' configure.ac  >/tmp/configure.ac".format(name))
os.system("cp /tmp/configure.ac configure.ac")

print "Update Makefile.am"
os.system("sed 's|pdfdecompositions/BasePdfDecomposition/src|pdfdecompositions/BasePdfDecomposition/src pdfdecompositions/{:s}PdfDecomposition/src|' Makefile.am > /tmp/Makefile.am".format(name))
os.system("cp /tmp/Makefile.am Makefile.am")

print "Update doxygen.cfg"
os.system("sed 's|pdfdecompositions/BasePdfDecomposition/include|pdfdecompositions/BasePdfDecomposition/include  pdfdecompositions/{:s}PdfDecomposition/include|' doxygen.cfg > /tmp/doxygen.cfg".format(name))
os.system("cp /tmp/doxygen.cfg  doxygen.cfg")

