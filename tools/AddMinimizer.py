#!/usr/bin/env python 
print "This code-generating script was written for an older version and is currently broken and needs to be rewritten.\nSorry!"
exit(-1)

''' Script to generate templates for new PDF decomposition '''

import sys
import os
import datetime

if len(sys.argv)<2:
    print ''' 
 Usage: AddMinimizer.py NAME
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


print "Creating directories in minimizers/"+name+"Minimizer"

os.system("mkdir -p minimizers/"+name+"Minimizer/include")
os.system("mkdir -p minimizers/"+name+"Minimizer/src")
os.system("mkdir -p minimizers/"+name+"Minimizer/yaml")
os.system("touch minimizers/"+name+"Minimizer/yaml/parameters.yaml")

hFile = "minimizers/{:s}Minimizer/include/{:s}Minimizer.h".format(name,name)

print "Creating header file  "+hFile


with open(hFile,"w+") as f:
    f.write(
        '''
#pragma once

#include "BaseMinimizer.h"

/**
  @class {NAME:s}Minimizer 

  @brief A class for {NAME:s} 

  @version 0.1
  @date {DATE:s}
  */

namespace xfitter {{

class {NAME:s}Minimizer : public BaseMinimizer
{{
  public:
     /// Default constructor. 
    {NAME:s}Minimizer () ;

     /// Default constructor. Name is the PDF name
    {NAME:s}Minimizer (const std::string& inName);

    /// Optional initialization at the first call
    virtual void initAtStart() override final;

    /// Miniminzation loop
    virtual void doMimimization() override final;

    /// Action at last iteration 
    virtual void actionAtFCN3() override final;

    /// Error analysis
    virtual void errorAnalysis() override final;

    /// Parameter transfer
    virtual void addParameter(double par, std::string const &name, double step = 0.01, double const* bounds = nullptr , double  const* priors  = nullptr ) override final;
}};
}}
'''.format( NAME=name, DATE= datetime.date.today().isoformat() )
)


sFile = "minimizers/{:s}Minimizer/src/{:s}Minimizer.cc".format(name,name)

print "Creating source file "+sFile

with open(sFile,"w+") as f:
    f.write(''' 
/*
   @file {NAME:s}Minimizer.cc
   @date {DATE:s}
   @author  AddMinimizer.py
   Created by  AddMinimizer.py on {DATE:s}
*/

#include "{NAME:s}Minimizer.h"

namespace xfitter {{
  
/// the class factories, for dynamic loading
extern "C" {NAME:s}Minimizer* create() {{
    return new {NAME:s}Minimizer();
}}


// Constructor
{NAME:s}Minimizer::{NAME:s}Minimizer() : BaseMinimizer("{NAME:s}") 
{{  
}}

// Constructor
{NAME:s}Minimizer::{NAME:s}Minimizer(const std::string& inName) : BaseMinimizer(inName) 
{{  
}}

// Init at start:
void {NAME:s}Minimizer::initAtStart() {{
  return;
}}

/// Miniminzation loop
void {NAME:s}Minimizer::doMimimization() 
{{
    return;
}}

/// Action at last iteration 
void {NAME:s}Minimizer::actionAtFCN3() 
{{
    return;
}}

/// Error analysis
void {NAME:s}Minimizer::errorAnalysis() 
{{
    return;
}}

/// Parameter transfer
void {NAME:s}Minimizer::addParameter(double par, std::string const &name, double step, double const* bounds , double  const* priors  )
{{
  BaseMinimizer::addParameter(par,name,step,bounds,priors);

  return; 
}}

}}

'''.format(NAME=name,DATE=datetime.date.today().isoformat())
)

    
aFile = "minimizers/{:s}Minimizer/src/Makefile.am".format(name)

print "Creating autoconf file " + aFile


with open(aFile,"w+") as f:
    f.write('''
# Created by AddMinimizer.py on {DATE:s}

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include   -I$(srcdir)/../../BaseMinimizer/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = lib{NAME:s}Minimizer_xfitter.la
lib{NAME:s}Minimizer_xfitter_la_SOURCES = {NAME:s}Minimizer.cc

datadir = ${{prefix}}/yaml/minimizers/{NAME:s}
data_DATA = ../yaml/parameters.yaml

dist_noinst_HEADERS = ../include ../yaml

lib{NAME:s}Minimizer_xfitter_la_LDFLAGS = -lbaseminimizer_xfitter -L$(libdir)

'''.format(DATE=datetime.date.today().isoformat(),NAME=name)
)


    
print "Update configure.ac file"
os.system("sed 's|xfitter-config|xfitter-config\\n		 minimizers/{:s}Minimizer/src/Makefile|' configure.ac  >/tmp/configure.ac".format(name))
os.system("cp /tmp/configure.ac configure.ac")

print "Update Makefile.am"
os.system("sed 's|minimizers/BaseMinimizer/src|minimizers/BaseMinimizer/src minimizers/{:s}Minimizer/src|' Makefile.am > /tmp/Makefile.am".format(name))
os.system("cp /tmp/Makefile.am Makefile.am")

print "Update doxygen.cfg"
os.system("sed 's|minimizers/BaseMinimizer/include|minimizers/BaseMinimizer/include  minimizers/{:s}Minimizer/include|' doxygen.cfg > /tmp/doxygen.cfg".format(name))
os.system("cp /tmp/doxygen.cfg  doxygen.cfg")

