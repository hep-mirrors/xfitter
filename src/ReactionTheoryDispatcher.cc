/*!
 @file ReactionTheoryDispatcher.cc
 @date Thu Jan 21 2016
 @author Andrey Sapronov <sapronov@ifh.de>

 Contains implementations of ReactionTheoryDispatcher class member functions.
 */

#include <list>
#include <string>
#include <cstring>

#include "ReactionTheoryDispatcher.h"
#include "ReactionTheory.h"
#include "ReactionDIS.h"
#include "xfitter_cpp.h"

using std::list;
using std::string;

ReactionTheoryDispatcher::~ReactionTheoryDispatcher()
{
}

ReactionTheory *ReactionTheoryDispatcher::getReactionTheory(const string &reaction_type)
{
  ReactionTheory *rt(NULL);
  if ( reaction_type == string("NC e+-p") ) rt = new ReactionDIS(string("NCDIS"));
  else if ( reaction_type == string("CC e+-p") ) rt = new ReactionDIS(string("CCDIS"));
  else if ( reaction_type == string("NC e+-p charm") ) rt = new ReactionDIS(string("CHARMDIS"));
  else if ( reaction_type == string("NC e+-p beauty") ) rt = new ReactionDIS(string("BEAUTYDIS"));
  else if ( reaction_type == string("NC e+-p FL") ) rt = new ReactionDIS(string("FL"));
  else if ( reaction_type == string("NC e+-p F2") ) rt = new ReactionDIS(string("F2"));
  else {
    int id = 16012101;
    char text[] = "S: Unknown reaction type in theory expression.";
    int textlen = strlen(text);
    hf_errlog_(id, text, textlen);

    return NULL;
  }

  _rt_list.push_back(rt);

  return rt;
}


void ReactionTheoryDispatcher::releaseReactionTheories()
{
  while (!_rt_list.size() ) {
    delete _rt_list.back();
    _rt_list.pop_back();
  }
}

