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

with open("Reactions.txt","r+") as f:
    for l in f:
        a = l.split()
        if a[0] == name:
            print "Interface for evolution "+name+" already exists, exit"
            exit(0)

# Not present, add new line to the Evolutions.txt file

with  open("Reactions.txt","a") as f:
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
  @class' Evolution{NAME:s}

  @brief A wrapper class for {NAME:s} evolution 

  @version 0.1
  @date {DATE:s}
  */

#include "BaseEvolution.h"

namespace xfitter 
{{

class Evolution{NAME:s} : BaseEvolution 
{{
  public:
    /// Empty constructor (needed for the dynamic loading)
    Evolution{NAME:s}():  BaseEvolution("{NAME:s}",nullptr) {{}};

  public:
  /// Constructor setting the name
    virtual std::string getEvolutionName() const {{ return  "{NAME:s}" ;}};
  /// Global initialization
    virtual void initAtStart() override final;
  /// Init at each iteration
    virtual void initAtIteration() override final;  

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)   
    virtual std::function<std::map<int,double>(double const& x, double const& Q)> xfxQMap() override final;

  /// Returns PDFs as a function of i, x, Q
    virtual std::function<double(int const& i, double const& x, double const& Q)> xfxQDouble() override final;
    
  /// Returns PDFs as double pdfs* --> double[13] from -6 to 6.  
    virtual std::function<void(double const& x, double const& Q, double* pdfs)> xfxQArray() override final;

  /// Returns alphaS
    virtual std::function<double(double const& Q)> AlphaQCD() override final;
}};

}};  // namespace xfitter

'''.format(NAME=name, DATE=datetime.date.today().isoformat())
)

print "Creating source file  evolutions/"+name+"/src/Evolution"+name+".cc"
with open("evolutions/"+name+"/src/Evolution"+name+".cc","w+") as f:
    f.write(''' 
/*
   @file Evolution{NAME:s}.cc
   @date  {DATE:s}
   @author  AddEvolution.py
   Created by  AddEvolution.py on {DATE:s}
*/

#include "Evolution{NAME:s}.h"

namespace xfitter
{{

// the class factories
extern "C" Evolution{NAME:s}* create() {{
  return new Evolution{NAME:s}();
}}

    
/// Global initialization
  void Evolution{NAME:s}::initAtStart() {{
    return ;
 }};
    
  /// Init at each iteration
  void Evolution{NAME:s}::initAtIteration() {{
    return ;
  }};

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)   
  std::function<std::map<int,double>(double const& x, double const& Q)> Evolution{NAME:s}::xfxQMap() {{
   }};

  /// Returns PDFs as a function of i, x, Q
  std::function<double(int const& i, double const& x, double const& Q)> Evolution{NAME:s}::xfxQDouble() {{
  }};
    
  /// Returns PDFs as double pdfs* --> double[13] from -6 to 6.  
  std::function<void(double const& x, double const& Q, double* pdfs)> Evolution{NAME:s}::xfxQArray() {{
  }};

  /// Returns alphaS
  std::function<double(double const& Q)> Evolution{NAME:s}::AlphaQCD() {{
  }};
}}
 
'''.format(NAME=name,DATE=datetime.date.today().isoformat())
) 


print "Creating autoconf file  evolutions/"+name+"/src/Makefile.am"
with open("evolutions/"+name+"/src/Makefile.am","w+") as f:
    f.write('''
# Created by AddEvolution.py on ''' + datetime.date.today().isoformat() + '''

AM_CXXFLAGS = -I$(srcdir)/../include -I$(srcdir)/../../BaseEvolution/include -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

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
