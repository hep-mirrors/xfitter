/*!
 @file ReactionTheory.cc
 @date Thu Jan 21 2016
 @author Andrey Sapronov <sapronov@ifh.de>

 Contains implementations of ReactionTheory class member functions.
 */

#include <list>
#include <iostream>
#include <string>

#include "ReactionTheory.h"

using std::list;
using std::string;

ReactionTheory::ReactionTheory(const ReactionTheory &rt)
{
  _val = new valarray<double>(rt._val->size());
}

ReactionTheory &
ReactionTheory::operator=(const ReactionTheory &rt)
{
  /*
  _subtype = rt._subtype;
  _ro = rt._ro;
  _binFlags = rt._binFlags;
  _dsBins = rt._dsBins;

  _val = new valarray<double>(*(rt._val));
  */

  return *this;
}

bool ReactionTheory::notMasked(int DSID, int Bin) {
  auto bins = _dsBins[DSID];
  auto flag = bins->find("binFlag");
  if ( flag == bins->end()) {  // DS has no bin "binFlag"
    return true;
  }
  else {
    return flag->second[Bin];
  }
}

string ReactionTheory::GetParamY(const string& name, int dsID ) const {
  if ( _xfitter_pars_node.find(name) != _xfitter_pars_node.end() ) {
  }
  else {
    return "";
  }
}
