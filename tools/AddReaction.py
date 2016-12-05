#!/usr/bin/env python 

''' Script to generate templates for a new theory module '''

import sys
import os
import datetime

if len(sys.argv)<2:
    print ''' 
 Usage: AddReaction.py NAME
    '''
    exit(0)

name = sys.argv[1]

# First check if the name is already used

with open("Reactions.txt","r+") as f:
    for l in f:
        a = l.split()
        if a[0] == name:
            print "Interface for reaction "+name+" already exist, exit"
            exit(0)

# Not present, add new line to the Reactions.txt file

with  open("Reactions.txt","a") as f:
    f.write(name+" "+"lib"+name.lower()+"_xfitter.so\n")

# Create directory structures:

print "Creating directories in reactions/"+name

os.system("mkdir -p reactions/"+name+"/include")
os.system("mkdir -p reactions/"+name+"/src")

print "Creating header file  reactions/"+name+"/include/Reaction"+name+".h"

with open("reactions/"+name+"/include/Reaction"+name+".h","w+") as f:
    f.write(
'''
#pragma once

#include "ReactionTheory.h"

/**
  @class' Reaction'''+name+'''

  @brief A wrapper class for '''+name+''' reaction 

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date ''' + datetime.date.today().isoformat() + '''
  */

class Reaction'''+name+''' : public ReactionTheory
{
  public:
    Reaction'''+name+'''(){};

//    ~Reaction'''+name+'''(){};
//    ~Reaction'''+name+'''(const Reaction'''+name+''' &){};
//    Reaction'''+name+''' & operator =(const ReactionA'''+name+''' &r){return *(new Reaction'''+name+'''(r));};

  public:
    virtual string getReactionName() const { return  "'''+name+ '''" ;};
    int compute();
  private:

    void initAtStart(const string &); 
    int parseOptions(){};
};

''')

print "Creating source file  reactions/"+name+"/src/Reaction"+name+".cc"
with open("reactions/"+name+"/src/Reaction"+name+".cc","w+") as f:
    f.write(''' 
/*
   @file Reaction'''+name+'''.cc
   @date ''' + datetime.date.today().isoformat() + '''
   @author  AddReaction.py
   Created by  AddReaction.py on ''' + datetime.date.today().isoformat() + '''
*/

#include "Reaction'''+name+'''.h"

// the class factories
extern "C" Reaction'''+name+'''* create() {
  return new Reaction'''+name+'''();
}


// Initialize at the start of the computation
void Reaction'''+name+'''::initAtStart(const string &s)
{
}

// Main function to compute results at an iteration
int Reaction'''+name+'''::compute()
{
}

''')


print "Creating source file  reactions/"+name+"/src/Makefile.am"
with open("reactions/"+name+"/src/Makefile.am","w+") as f:
    f.write('''
# Created by AddReaction.py on ''' + datetime.date.today().isoformat() + '''

AM_CXXFLAGS = -I$(srcdir)/../include  -I$(srcdir)/../../../include  -I$(srcdir)/../../../interfaces/include -Wall -fPIC -Wno-deprecated 

lib_LTLIBRARIES = lib'''+ name.lower() + '''xfitter.la
lib'''+ name.lower()+'''xfitter_la_SOURCES = Reaction'''+name+'''.cc

# lib'''+ name.lower()+'''xfitter_la_LDFLAGS = place_if_needed  

 ''')


print "Update configure.ac file"
os.system("sed 's|xfitter-config|xfitter-config\\n		 reactions/" +name +"/src/Makefile|' configure.ac  >/tmp/configure.ac")
# os.system("cp /tmp/configure.ac configure.ac")

print "Update Makefile.am"
os.system("sed 's|tools/process|tools/process reactions/" +name +"/src|' Makefile.am > /tmp/Makefile.am")
# os.system("cp /tmp/Makefile.am Makefile.am")
