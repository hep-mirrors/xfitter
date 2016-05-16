/*!
 @file ReactionTheoryDispatcher.cc
 @date Thu Jan 21 2016
 @author Andrey Sapronov <sapronov@ifh.de>

 Contains implementations of ReactionTheoryDispatcher class member functions.
 */

#include <list>
#include <string>

#include "ReactionTheoryDispatcher.h"
#include "ReactionTheory.h"
#include "Reaction_FTDY_NC.h"

using std::list;
using std::string;

ReactionTheoryDispatcher::~ReactionTheoryDispatcher()
{
}

ReactionTheory &ReactionTheoryDispatcher::getReactionTheory(const string &reaction_type)
{
  ReactionTheory *rt(NULL);
  if ( reaction_type == string("FTDY_NC_pp") ) rt = new Reaction_FTDY_NC(string("pp"));
  else if ( reaction_type == string("FTDY_NC_pn") ) rt = new Reaction_FTDY_NC(string("pn"));
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


ReactionTheory &ReactionTheoryDispatcher::releaseReactionTheories()
{
  while (!_rt_list.size() ) do {
    delete _rt_list.back();
    _rt_list.pop_back();
  }
}

