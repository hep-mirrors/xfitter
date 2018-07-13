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

// Global variable to hold current function

std::function<void(double const& x, double const& Q, double* pdfs)> gProtonPdf;

void protonPDF(double const& x, double const& Q, double* pdfs) {
  gProtonPdf(x,Q,pdfs);
}


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

void ReactionTheory::initAtIteration() {
  // do some magic
  gProtonPdf = XFITTER_PARS::retrieveXfxQArray("APFELxx:p");
}

const pXFXlike  ReactionTheory::getXFX(const string& type) {
  //  gProtonPdf = XFITTER_PARS::retrieveXfxQArray("APFELxx:p");

  // return _xfx[type];
  std::cout << " here "   <<std::endl;
  double dd[13];
  gProtonPdf(0.01,1.,dd);
  std::cout << " ho " << dd[6] << std::endl;


  
  return &protonPDF;
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

void ReactionTheory::resetParameters(const YAML::Node& node) {
  XFITTER_PARS::parse_node( node, _xfitter_pars, _xfitter_pars_i, _xfitter_pars_s, _xfitter_pars_vec, _xfitter_pars_node);
}

/// Dump local parameters for the reaction. Update to the current state for the fitted parameters.
std::string ReactionTheory::emitReactionLocalPars() const {

  if ( _xfitter_pars_node.find(getReactionName()) ==  _xfitter_pars_node.end()) {
    return "";
  }


  YAML::Emitter out;
  YAML::Node pars = _xfitter_pars_node.at(getReactionName());

  // update pars:
  for ( YAML::iterator it = pars.begin(); it != pars.end(); ++it) {
    auto value = it->second;
    string name  = (it->first).as<string>();
    if ( value.IsMap() ) {
      if (value["step"]) {
	// Update:
	value["value"] = GetParam(name);
      }
    }
  }
  // Dump with prefix:
  YAML::Node n;
  n[getReactionName()] = pars;

  out.SetIndent(2);
  out.SetSeqFormat(YAML::Flow);

  out << n;
  return out.c_str();
}


