#include <iostream>
#include "TermData.h"
#include "ReactionTheory.h"
#include "TheorEval.h"
#include "xfitter_cpp_base.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include "BaseEvolution.h"
using namespace std;
using xfitter::BaseEvolution;
//Wrappers
BaseEvolution *wrappedPDFs[2];
extern "C"
{
  void PDF_xfxQ_wrapper(const double &x, const double &Q, double *r)
  {
    wrappedPDFs[0]->xfxQarray(x, Q, r);
  }
  void PDF_xfxQ_wrapper1(const double &x, const double &Q, double *r)
  {
    wrappedPDFs[1]->xfxQarray(x, Q, r);
  }
  double AlphaS_wrapper(const double &Q)
  {
    return wrappedPDFs[0]->getAlphaS(Q);
  }

  // Add fortran versions:
  void pdf_xfxq_wrapper_(const double &x, const double &Q, double *r) {
    PDF_xfxQ_wrapper(x,Q,r);
  }
  void pdf_xfxq_wrapper1_(const double &x, const double &Q, double *r) {
    PDF_xfxQ_wrapper1(x,Q,r);
  }

  double alphas_wrapper_(const double &Q) {
    return AlphaS_wrapper(Q);
  }
}

//Note: there was some TermInfo-parsing logic in TheorEval and ftheor_eval code, but I did not like it and wrote my own. Now it can handle spaces!
void splitTermInfo(map<string, string> &out, const char *s)
{
  //s is string from TermInfo
  //out is the map to fill with key=value pairs
  //Both ':' or ';' can be used to separate pairs
  //keys and values can have spaces in them, for example
  //"  key with  spaces  =  value with spaces   "
  //parses into
  //"key with spaces"->"value with spaces"
  const char *p = s; //current position
  while (true)
  {
    const char *k0, *k1; //beginning and end of current key
    const char *v0, *v1; //beginning and end of current value
    while (*p == ' ')
      ++p; //skip whitespace leading to key
    if (*p == 0)
      return;
    if (*p == ':' || *p == ';')
    {
      cerr << "[ERROR] No key before separator in TermInfo=\"" << s << "\" (position " << p - s << ")\n"
                                                                                                   "Expected format: TermInfo='key1=value1;key2=value2;...'"
           << endl;
      hf_errlog(18041201, "F: Bad TermInfo format, see stderr");
    }
    if (*p == '=')
    {
      cerr << "[ERROR] Empty key in TermInfo=\"" << s << "\" (position " << p - s << ")\n"
                                                                                     "Expected format: TermInfo='key1=value1;key2=value2;...'"
           << endl;
      hf_errlog(18041201, "F: Bad TermInfo format, see stderr");
    }
    k0 = p;
    while (*p != '=')
    {
      if (*p == 0 || *p == ':' || *p == ';')
      {
        cerr << "[ERROR] No '=' sign found for a key in TermInfo=\"" << s << "\" (key starts at position " << k0 - s << ")\n"
                                                                                                                        "Expected format: TermInfo='key1=value1;key2=value2;...'"
             << endl;
        hf_errlog(18041201, "F: Bad TermInfo format, see stderr");
      }
      ++p;
    }
    k1 = p - 1;
    ++p;
    while (*k1 == ' ')
      --k1; //trim trailing whitespace in key
    ++k1;
    while (*p == ' ')
      ++p; //skip whitespace leading to value
    if (*p == 0 || *p == ':' || *p == ';' || *p == '=')
    {
      cerr << "[ERROR] No value found after key \"" << string(k0, k1 - k0) << "\" in TermInfo=\"" << s << "\"\n"
                                                                                                          "Expected format: TermInfo='key1=value1;key2=value2;...'"
           << endl;
      hf_errlog(18041201, "F: Bad TermInfo format, see stderr");
    }
    v0 = p;
    ++p;
    while (*p != ':' && *p != ';' && *p)
    {
      if (*p == '=')
      {
        cerr << "[ERROR] Unexpected '=' sign after value corresponding to key \"" << string(k0, k1 - k0) << "\" in TermInfo=\"" << s << "\"\n"
                                                                                                                                        "Expected format: TermInfo='key1=value1;key2=value2;...'"
             << endl;
        hf_errlog(18041201, "F: Bad TermInfo format, see stderr");
      }
      ++p;
    }
    v1 = p - 1;
    if (*p)
      ++p;
    while (*v1 == ' ')
      --v1; //trim trailing whitespace in value
    ++v1;
    string key = string(k0, k1 - k0);
    if (out.count(key) != 0)
    {
      cerr << "[ERROR] Duplicate key \"" << key << "\" in TermInfo=\"" << s << "\"" << endl;
      hf_errlog(18041200, "F: Duplicate key in a TermInfo, see stderr");
    }
    out[key] = string(v0, v1 - v0);
  }
  //NOTE: only spaces are treated as whitespace, tabs and other stuff is not handled
}
const char *const BY_REACTION = "byReaction";
using XFITTER_PARS::createConstantParameter;
using XFITTER_PARS::gParameters;
using XFITTER_PARS::rootNode;
YAML::Node getFromByReaction(const string &parName, const ReactionTheory *const reaction)
{
  //Return a node corresponding to reation-specific parameter
  //Or, if such a node is not found, return an invalid(undefined) node
  YAML::Node byReactionNode = rootNode[BY_REACTION];
  if (byReactionNode.IsMap())
  {
    YAML::Node reactionNode = byReactionNode[reaction->getReactionName()];
    if (reactionNode.IsMap())
      return reactionNode[parName];
    else if (reactionNode.IsDefined())
    {
      cerr << "[ERROR] Reaction-specific parameters node " << BY_REACTION << '/' << reaction->getReactionName() << " is not a YAML map" << endl;
      hf_errlog(19041301, "F: Reaction parameter node is not a map, see stderr");
    }
  }
  else
  { //Check in case somebody still uses the old syntax
    //When everyone has migrated, you can remove this whole block
    static bool once = true;
    if (once)
    {
      YAML::Node oldstyleNode = XFITTER_PARS::rootNode[reaction->getReactionName()];
      if (oldstyleNode.IsDefined())
      {
        once = false;
        cerr << "[WARN] When searching for parameter \"" << parName << "\" requested by reaction \"" << reaction->getReactionName() << "\" I noticed that there is no \"" << BY_REACTION << "\" node, but there is a node named \"" << reaction->getReactionName() << "\" at YAML root. "
                                                                                                                                                                                                                                                                      "It seems you are using the old syntax to pass reaction-specific parameters. By the new syntax, which was introduced during the \"TermData rewrite\", reaction node should be under \""
             << BY_REACTION << "\", like this:\n"
             << BY_REACTION << ":\n  " << reaction->getReactionName() << ":\n    #YOUR PARAMETERS HERE\n  some_other_reaction:\n    #other parameters\n  ...\n"
                                                                         "I am NOT reading old-style reaction-specific parameters this time, go change your YAML input. I will try to get this parameter from YAML root now.\n"
                                                                         "This type of warning will not be issued again."
             << endl;
        hf_errlog(19041302, "W: Old-style reaction parameter syntax detected, see stderr");
      }
    }
  }
  return YAML::Node(YAML::NodeType::Undefined);
}
enum class Type
{
  None = 0,
  DoublePtr = 1,
  String = 2,
  Int = 3,
  Array = 4
}; //None is unused
const char *to_cstring(Type t)
{
  static const char *a[] = {"None", "double*", "string", "int", "vector<double*>"};
  return a[int(t)];
}
std::ostream &operator<<(std::ostream &os, Type t) { return os << to_cstring(t); }
enum class ParamScope
{
  Term = 0,
  Reaction = 1,
  Global = 2
};
void reportUndefinedParameter[[noreturn]](const string &parName, const Type type, const ReactionTheory *reaction)
{
  cerr << "[ERROR] Undefined parameter \"" << parName << "\" requested as " << type << " by reaction \"" << reaction->getReactionName() << "\"" << endl;
  hf_errlog(19041303, "F: Reaction requested undefined parameter, see stderr");
  abort(); //unreachable
}
void reportFailedConversion[[noreturn]](const string &parName, const Type type, const ParamScope scope, const ReactionTheory *reaction, const YAML::Node *node = nullptr, const string *definition = nullptr)
{
  cerr << "[ERROR] Failed to convert parameter \"" << parName << "\" to " << type << ".\nThis parameter was requested as " << type << " by reaction \"" << reaction->getReactionName() << "\", and was found in ";
  switch (scope)
  {
  case ParamScope::Term:
    cerr << "corresponding TermInfo";
    break;
  case ParamScope::Reaction:
    cerr << "YAML steering as " << BY_REACTION << '/' << reaction->getReactionName() << '/' << parName;
    break;
  case ParamScope::Global:
    cerr << "root of YAML steering";
    break;
  }
  if (node)
  {
    cerr << ". Trying to print the node:" << endl;
    cerr << *node << endl;
  }
  if (definition)
  { //for TermInfo
    cerr << ". Parameter definition: \"" << parName << '=' << *definition << "\"" << endl;
  }
  cerr << "[/ERROR]" << endl;
  hf_errlog(19041300, "F: Failed to convert parameter to requested type, see stderr");
  abort(); //unreachable
}
const char *const PAR_TAG = "!parameter";
void reportUndefinedReference[[noreturn]](const string &parName, const string &refName, const ParamScope scope, const ReactionTheory *reaction)
{
  cerr << "[ERROR] Reference to undefined parameter \"" << refName << "\" in ";
  switch (scope)
  {
  case ParamScope::Term:
    cerr << "TermInfo parameter \"" << parName << "\" for reaction \"" << reaction->getReactionName() << '\"';
    break;
  case ParamScope::Reaction:
    cerr << "YAML steering at node " << BY_REACTION << '/' << reaction->getReactionName() << '/' << parName;
    break;
  case ParamScope::Global:
    cerr << "parameter \"" << parName << "\" at root of YAML steering";
    break;
  }
  cerr << endl;
  hf_errlog(19041710, "F: Reference to undefined parameter, see stderr");
  abort(); //unreachable
}
//TODO: Implement TheorEval managing instances of TermData
TermData::TermData(unsigned _id, ReactionTheory *_reaction, TheorEval *_parent, const char *term_info_string) : id{_id}, reaction{_reaction}, parent{_parent}
{
  splitTermInfo(term_info, term_info_string); //fills term_info
}
int TermData::getStringFromTermOrReaction(const string &parName, string *out)
{
  //On sucess return 0 and set out
  //On failure return 1
  //1. TermInfo from datafile
  {
    const auto it = term_info.find(parName);
    if (it != term_info.end())
    {
      *out = it->second;
      return 0;
    }
  }
  //2. Reaction-specific from yaml steering
  {
    YAML::Node parNode = getFromByReaction(parName, reaction);
    if (parNode.IsDefined())
    {
      try
      {
        *out = parNode.as<string>();
      }
      catch (const YAML::TypedBadConversion<string> &ex)
      {
        reportFailedConversion(parName, Type::String, ParamScope::Reaction, reaction, &parNode);
      }
      return 0;
    }
  }
  return 1;
}
string TermData::getParamS(const string &parName)
{
  //Try to get the parameter from the following places, in the following order/priority:
  //1. Term-specific from datafile
  //2. Reaction-specfic from yaml steering
  //3. Global from yaml steering

  //1, 2. Term and Reaction
  {
    string out;
    if (getStringFromTermOrReaction(parName, &out) == 0)
      return out;
  }
  //3. Global from YAML
  //NOTE: I do not use gParametersI, gParametersS, gParametersV, gParametersY because I do not like them --Ivan
  {
    YAML::Node parNode = rootNode[parName];
    if (parNode.IsDefined())
    {
      try
      {
        return parNode.as<string>();
      }
      catch (const YAML::TypedBadConversion<string> &ex)
      {
        reportFailedConversion(parName, Type::String, ParamScope::Global, reaction, &parNode);
      }
    }
  }
  reportUndefinedParameter(parName, Type::String, reaction);
}
const double *TermData::getParamD(const string &parName)
{
  //1. TermInfo from datafile
  {
    const auto it = term_info.find(parName);
    if (it != term_info.end())
    {
      //Try reference to parameter defined in YAML: "!parameter Bv"
      string definition = it->second;
#define PAR_PREFIX "!parameter "
      if (beginsWith(definition, PAR_PREFIX))
      {
        size_t p = definition.find_first_not_of(' ', sizeof(PAR_PREFIX) - 1); //because definition cannot end in ' ', this will always succeed
        string refName = definition.substr(p);
        auto itp = gParameters.find(refName);
        if (itp == gParameters.end())
        {
          reportUndefinedReference(parName, refName, ParamScope::Term, reaction);
        }
        return itp->second;
      }
#undef PAR_PREFIX
      //Try constant: "10" or "2.0", then create a constant parameter specifically for this value
      //Name it "termID/NAME", with some ID and NAME
      const string parFullname = "term" + to_string(id) + '/' + parName;
      //First check if this constant parameter has already been created
      const auto itp = gParameters.find(parFullname);
      if (itp != gParameters.end())
        return itp->second;
      //Else try converting to double:
      char *endp;
      double value = strtod(definition.c_str(), &endp);
      if (endp == definition.c_str())
      {
        reportFailedConversion(parName, Type::DoublePtr, ParamScope::Term, reaction, nullptr, &definition);
      }
      return createConstantParameter(parFullname, value);
    }
  }
  //2. Reaction-specfic from yaml steering
  {
    YAML::Node parNode = getFromByReaction(parName, reaction);
    if (parNode.IsDefined())
    {
      //2.1 reference to another parameter, like
      //alphas: !parameter fitted_alphas
      if (parNode.Tag() == PAR_TAG)
      {
        try
        {
          string refName = parNode.as<string>();
          const auto it = gParameters.find(refName);
          if (it == gParameters.end())
          {
            reportUndefinedReference(parName, refName, ParamScope::Reaction, reaction);
          }
          return it->second;
        }
        catch (YAML::TypedBadConversion<string>)
        {
          cerr << "[ERROR] Failed to convert name of reference parameter to string in reference\n"
               << BY_REACTION << ":\n  " << reaction->getReactionName() << ":\n    " << parName << ": " << PAR_TAG << " ???\nTrying to print it:" << endl;
          cerr << parNode << endl;
          cerr << "[/ERROR]";
          hf_errlog(19041700, "F: Failed to convert parameter name to string, see stderr");
        }
      }
      //2.2 constant
      string parFullname = reaction->getReactionName() + "/" + parName;
      //first check if it exists already
      const auto it = gParameters.find(parFullname);
      if (it != gParameters.end())
        return it->second;
      //else create it
      try
      {
        return createConstantParameter(parFullname, parNode.as<double>());
      }
      catch (YAML::TypedBadConversion<double>)
      {
        reportFailedConversion(parName, Type::DoublePtr, ParamScope::Reaction, reaction, &parNode);
      }
    }
  }
  //3. Global
  {
    //First check if it exists already
    const auto it = gParameters.find(parName);
    if (it != gParameters.end())
      return it->second;
    //Else try to create a new constant parameter
    const YAML::Node parNode = rootNode[parName];
    if (parNode.IsDefined())
    {
      try
      {
        return createConstantParameter(parName, parNode.as<double>());
      }
      catch (YAML::TypedBadConversion<double>)
      {
        reportFailedConversion(parName, Type::DoublePtr, ParamScope::Global, reaction, &parNode);
      }
    }
  }
  reportUndefinedParameter(parName, Type::DoublePtr, reaction);
}
int TermData::getParamI(const string &parName)
{
  //1. TermInfo from datafile
  {
    const auto it = term_info.find(parName);
    if (it != term_info.end())
    {
      try
      {
        return stoi(it->second);
      }
      catch (...)
      {
        reportFailedConversion(parName, Type::Int, ParamScope::Term, reaction, nullptr, &it->second);
      }
    }
  }
  //2. Reaction-specific from yaml steering
  {
    YAML::Node parNode = getFromByReaction(parName, reaction);
    if (parNode.IsDefined())
    {
      try
      {
        return parNode.as<int>();
      }
      catch (const YAML::TypedBadConversion<int> &ex)
      {
        reportFailedConversion(parName, Type::Int, ParamScope::Reaction, reaction, &parNode);
      }
    }
  }
  //3. Global from YAML
  {
    YAML::Node parNode = rootNode[parName];
    if (parNode.IsDefined())
    {
      try
      {
        return parNode.as<int>();
      }
      catch (const YAML::TypedBadConversion<int> &ex)
      {
        reportFailedConversion(parName, Type::Int, ParamScope::Global, reaction, &parNode);
      }
    }
  }
  reportUndefinedParameter(parName, Type::Int, reaction);
}
bool TermData::hasParam(const string &parName)
{
  //1. TermInfo from datafile
  if (term_info.count(parName) > 0)
    return true;
  { //2. Reaction-specific from yaml steering
    YAML::Node parNode = getFromByReaction(parName, reaction);
    if (parNode.IsDefined())
      return true;
  }
  { //3. Global from YAML
    YAML::Node parNode = rootNode[parName];
    if (parNode.IsDefined())
      return true;
  }
  return false;
}
BaseEvolution *TermData::getPDF(int ind)
{
  using xfitter::get_evolution;
  string pdfName;
  if (ind == 0)
  {
    if (getStringFromTermOrReaction("evolution", &pdfName) == 0)
    {
      return get_evolution(pdfName);
    }
  }
  string parName = "evolution" + to_string(ind + 1);
  if (getStringFromTermOrReaction(parName, &pdfName) == 0)
  {
    return get_evolution(pdfName);
  }
  return get_evolution(); //returns default evolution, which always exists
}
void TermData::actualizeWrappers()
{
  //TODO: possible to speed up by caching calls to getPDF
  wrappedPDFs[0] = getPDF(0);
  wrappedPDFs[1] = getPDF(1);
}
const valarray<double> &TermData::getBinColumn(const string &n)
{
  const valarray<double> *ret = parent->getBinColumn(n);
  if (ret)
    return *ret;
  cerr << "[ERROR] Reaction \"" << reaction->getReactionName() << "\" requested bins of nonexistant column \"" << n << "\"" << endl;
  hf_errlog(19042700, "F: Reaction requested nonexistant bin column, see stderr");
  abort(); //unreachable
}
