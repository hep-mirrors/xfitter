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
  void pdf_xfxq_wrapper_(const double &x, const double &Q, double *r) {
    wrappedPDFs[0]->xfxQarray(x, Q, r);
  }
  void pdf_xfxq_wrapper1_(const double &x, const double &Q, double *r) {
    wrappedPDFs[1]->xfxQarray(x, Q, r);
  }

  double alphas_wrapper_(const double &Q) {
    return wrappedPDFs[0]->getAlphaS(Q);
  }

  // also for apfel:
  void externalsetapfel1_(double &x, double &Q, double *xf)
  {
    pdf_xfxq_wrapper_(x,Q,xf);
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
const char *const BY_DATASET = "byDataset";
using XFITTER_PARS::createConstantParameter;
using XFITTER_PARS::gParameters;
using XFITTER_PARS::rootNode;
YAML::Node getFromByDataset(const string& parName, const string& datasetName) {
  //Return a node corresponding to dataset-specific parameter in the YAML steering
  //Or, if such a node is not found, return an invalid(undefined) node
  YAML::Node byDatasetNode = rootNode[BY_DATASET];
  if (byDatasetNode.IsMap())
  {
    YAML::Node reactionNode = byDatasetNode[datasetName];
    if (reactionNode.IsMap()) return reactionNode[parName];
    else if (reactionNode.IsDefined())
    {
      cerr<< "[ERROR] Dataset-specific parameters node " << BY_DATASET << '/' << datasetName << " is not a YAML map" << endl;
      hf_errlog(20032800, "F: Dataset parameter node is not a map, see stderr");
    }
  }
  return YAML::Node(YAML::NodeType::Undefined);
}
YAML::Node getFromByReaction(const string& parName, const string& reaction_name) {
  //Return a node corresponding to reaction-specific parameter
  //Or, if such a node is not found, return an invalid(undefined) node
  YAML::Node byReactionNode = rootNode[BY_REACTION];
  if (byReactionNode.IsMap()) {
    YAML::Node reactionNode = byReactionNode[reaction_name];
    if (reactionNode.IsMap()) return reactionNode[parName];
    else if (reactionNode.IsDefined()) {
      cerr << "[ERROR] Reaction-specific parameters node " << BY_REACTION << '/' << reactionNode << " is not a YAML map" << endl;
      hf_errlog(19041301, "F: Reaction parameter node is not a map, see stderr");
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
  static const char* a[] = {"None", "double*", "string", "int", "vector<double*>"};
  return a[int(t)];
}
std::ostream &operator<<(std::ostream &os, Type t) { return os << to_cstring(t); }
enum class ParamScope
{
  Dataset = 0, //Dataset-specific override from yaml steering
  Term = 1, //Term-specific from TermInfo in datafile
  Reaction = 2, //Reaction-specfic from yaml steering
  Global = 3, //Global from yaml steering
  Undefined = 4
  //Smaller number means higher priority
  //Scope with smaller number overrides scope with larger number
};
ParamScope getParameterScope(const string& parName, const map<string,string>& term_info, const string& dataset_name, const string& reaction_name){
  if (getFromByDataset(parName, dataset_name).IsDefined()) return ParamScope::Dataset;
  if (term_info.count(parName) > 0) return ParamScope::Term;
  if (getFromByReaction(parName, reaction_name).IsDefined()) return ParamScope::Reaction;
  if (gParameters.count(parName) > 0 or rootNode[parName].IsDefined()) return ParamScope::Global;
  return ParamScope::Undefined;
}
/// \brief Print a debug message to cerr describing where the parameter was found
void printScopeDescription(
  const string& parameter_name,
  const ParamScope scope,
  const Type type,
  const string& dataset_name,
  const string& reaction_name
){
  cerr<<
    "Parameter \""<<parameter_name<<"\" "
    "was requested as type \""<<type<<"\" "
    "by reaction \""<<reaction_name<<"\" "
    "for dataset \""<<dataset_name<<"\". "
    "The parameter is ";
  if (scope == ParamScope::Undefined) {
    cerr<<"undefined";
    return;
  }
  cerr<<"defined in the ";
  if (scope == ParamScope::Term) {
    cerr<<"corresponding TermInfo";
    return;
  }
  cerr<<"YAML steering at ";
  switch (scope) {
  case ParamScope::Dataset:
    cerr<<BY_DATASET<<"/\""<<dataset_name<<"\"/"<<parameter_name;
    break;
  case ParamScope::Reaction:
    cerr<<BY_REACTION<<'/'<<reaction_name<<'/'<<parameter_name;
    break;
  case ParamScope::Global:
    cerr<<"global scope";
    break;
  default:break;
  }
}
void reportUndefinedParameter[[noreturn]](
  const string& parameter_name,
  const Type type,
  const string& dataset_name,
  const string& reaction_name
){
  cerr<<"[ERROR] ";
  printScopeDescription(parameter_name, ParamScope::Undefined, type, dataset_name, reaction_name);
  cerr<<endl<<"[/ERROR]"<<endl;
  hf_errlog(19041303, "F: Reaction requested undefined parameter, see stderr"); //does not return
  abort(); //unreachable
}
void reportFailedConversion[[noreturn]](
  const string& parameter_name,
  const Type type,
  const ParamScope scope,
  const string& dataset_name,
  const string& reaction_name,
  const YAML::Node* node = nullptr,
  const string* definition = nullptr
){
  cerr<<"[ERROR] Type conversion failed"<<endl;
  printScopeDescription(parameter_name, scope, type, dataset_name, reaction_name);
  cerr<<endl;
  if (node) {
    cerr << "Trying to print the node:" << endl;
    cerr << *node << endl;
  }
  if (definition) { //for TermInfo
    cerr << "Parameter definition: \"" << parameter_name << '=' << *definition << "\"" << endl;
  }
  cerr << "[/ERROR]" << endl;
  hf_errlog(19041300, "F: Failed to convert parameter to requested type, see stderr");
  abort(); //unreachable
}
static const char* const PAR_TAG = "!parameter";
void reportUndefinedReference[[noreturn]](
  const string& parameter_name,
  const string& reference_name,
  const ParamScope scope,
  const string& dataset_name,
  const string& reaction_name
){
  cerr<<"[ERROR] Parameter \""<<parameter_name<<"\" is a reference to undefined parameter \""<<reference_name<<'\"'<<endl;
  printScopeDescription(parameter_name, scope, Type::DoublePtr, dataset_name, reaction_name);
  cerr<<"\n[/ERROR]"<<endl;
  hf_errlog(19041710, "F: Reference to undefined parameter, see stderr"); //does not return
  abort(); //unreachable
}
TermData::TermData(unsigned _id, ReactionTheory *_reaction, TheorEval *_parent, const char *term_info_string) : id{_id}, reaction{_reaction}, parent{_parent} {
  splitTermInfo(term_info, term_info_string); //fills term_info
}
string TermData::getParamS(const string& parameter_name) {
  const Type type = Type::String;
  const string& dataset_name = parent->_ds_name;
  const string& reaction_name = reaction->getReactionName();
  ParamScope scope = getParameterScope(parameter_name, term_info, dataset_name, reaction_name);
  if (scope == ParamScope::Undefined) reportUndefinedParameter(parameter_name, type, dataset_name, reaction_name);
  if (scope == ParamScope::Term) return term_info.at(parameter_name);
  YAML::Node node;
  if (scope == ParamScope::Dataset) node = getFromByDataset(parameter_name, dataset_name);
  else if (scope == ParamScope::Reaction) node = getFromByReaction(parameter_name, reaction_name);
  else if (scope == ParamScope::Global) node = rootNode[parameter_name];
  else {//this branch is unreachable
    hf_errlog(20040400, "F: Programming error: unhandled parameter scope");
    abort();
  }
  try{
    return node.as<string>();
  } catch (const YAML::TypedBadConversion<string> &ex) {
    reportFailedConversion(parameter_name, type, scope, dataset_name, reaction_name, &node);
  }
}
const double* TermData::getParamD(const string& parameter_name) {
  const Type type = Type::DoublePtr;
  const string& dataset_name = parent->_ds_name;
  const string& reaction_name = reaction->getReactionName();
  ParamScope scope = getParameterScope(parameter_name, term_info, dataset_name, reaction_name);
  if (scope == ParamScope::Undefined) reportUndefinedParameter(parameter_name, type, dataset_name, reaction_name);
  if (scope == ParamScope::Term) {
    const string& definition = term_info.at(parameter_name);
    //Try reference to parameter defined in YAML: "!parameter Bv"
    static const char REFERENCE_TAG[] = "!parameter ";
    if (beginsWith(definition, REFERENCE_TAG)) {
      const size_t p = definition.find_first_not_of(' ', sizeof(REFERENCE_TAG) - 1); //because definition cannot end in ' ', this will always succeed
      const string& reference_name = definition.substr(p);
      auto itp = gParameters.find(reference_name);
      if (itp == gParameters.end()) reportUndefinedReference(parameter_name, reference_name, scope, dataset_name, reaction_name);
      return itp->second;
    }
    //Try constant: "10" or "2.0", then create a constant parameter specifically for this value
    //Name it "termID/NAME", with some ID and NAME
    const string full_name = "term" + to_string(id) + '/' + parameter_name;

    /*
    //First check if this constant parameter has already been created
    const auto it = gParameters.find(full_name);
    if (it != gParameters.end()) return it->second;

    //Else try converting to double:
    char *endp;
    double value = strtod(definition.c_str(), &endp);
    if (endp == definition.c_str()) reportFailedConversion(parameter_name, type, scope, dataset_name, reaction_name, nullptr, &definition);
    */

    //First try converting to double:
    char *endp;
    double value = strtod(definition.c_str(), &endp);
    if (endp == definition.c_str()) reportFailedConversion(parameter_name, type, scope, dataset_name, reaction_name, nullptr, &definition);

    //Then check if this constant parameter has already been created
    const auto it = gParameters.find(full_name);
    if (it != gParameters.end())
      {
	//update the value
	*(it->second) = value;
	return it->second;
      }
    
    //Else create a new constant parameter
    return createConstantParameter(full_name, value);
  }
  YAML::Node node;
  if (scope == ParamScope::Dataset) node = getFromByDataset(parameter_name, dataset_name);
  else if (scope == ParamScope::Reaction) node = getFromByReaction(parameter_name, reaction_name);
  else if (scope == ParamScope::Global) node = rootNode[parameter_name];
  else {//this branch is unreachable
    hf_errlog(20040400, "F: Programming error: unhandled parameter scope");
    abort();
  }
  try{
    //Case 1: reference to another double-typed parameter
    if (scope != ParamScope::Global && node.Tag() == PAR_TAG) {
      const string reference_name = node.as<string>();
      try{
        return gParameters.at(reference_name);
      } catch (std::out_of_range& ex){
        reportUndefinedReference(parameter_name, reference_name, scope, dataset_name, reaction_name);
      }
    }
    //Case 2: constant double
    string full_name;
    if (scope == ParamScope::Dataset) full_name = "dataset" + dataset_name + '/' + parameter_name;
    else if (scope == ParamScope::Reaction) full_name = reaction_name + '/' + parameter_name;
    else if (scope == ParamScope::Global) full_name = parameter_name;
    else {//this branch is unreachable
      hf_errlog(20040401, "F: Programming error: unhandled parameter scope");
      abort();
    }
    //Check if already exists
    const auto it = gParameters.find(full_name);
    if (it != gParameters.end()) return it->second;
    //else create and return a new constant parameter
    return createConstantParameter(full_name, node.as<double>());
  } catch (const YAML::TypedBadConversion<double> &ex) {
    reportFailedConversion(parameter_name, Type::DoublePtr, scope, dataset_name, reaction_name, &node);
  } catch (const YAML::TypedBadConversion<string> &ex) {
    reportFailedConversion(parameter_name, Type::String, scope, dataset_name, reaction_name, &node);
  }
}
int TermData::getParamI(const string& parameter_name) {
  //Copy-pasted from getParamS with minimal changes
  const Type type = Type::Int;
  const string& dataset_name = parent->_ds_name;
  const string& reaction_name = reaction->getReactionName();
  ParamScope scope = getParameterScope(parameter_name, term_info, dataset_name, reaction_name);
  if (scope == ParamScope::Undefined) reportUndefinedParameter(parameter_name, type, dataset_name, reaction_name);
  if (scope == ParamScope::Term) {
    string definition = term_info.at(parameter_name);
    try{
      return std::stoi(definition);
    } catch (...) {
      reportFailedConversion(parameter_name, type, scope, dataset_name, reaction_name, nullptr, &definition);
    }
  }
  YAML::Node node;
  if (scope == ParamScope::Dataset) node = getFromByDataset(parameter_name, dataset_name);
  else if (scope == ParamScope::Reaction) node = getFromByReaction(parameter_name, reaction_name);
  else if (scope == ParamScope::Global) node = rootNode[parameter_name];
  else {//this branch is unreachable
    hf_errlog(20040400, "F: Programming error: unhandled parameter scope");
    abort();
  }
  try{
    return node.as<int>();
  } catch (const YAML::TypedBadConversion<int> &ex) {
    reportFailedConversion(parameter_name, type, scope, dataset_name, reaction_name, &node);
  }
}
bool TermData::hasParam(const string& parName) {
  return getParameterScope(parName, term_info, parent->_ds_name, reaction->getReactionName()) != ParamScope::Undefined;
}
BaseEvolution* TermData::getPDF(int ind) {
  using xfitter::get_evolution;
  if (ind == 0 and hasParam("evolution")) return get_evolution(getParamS("evolution"));
  string parName = "evolution" + to_string(ind + 1);
  if (hasParam(parName)) return get_evolution(getParamS(parName));
  return get_evolution(); //returns default evolution
}
void TermData::actualizeWrappers() {
  wrappedPDFs[0] = getPDF(0);
  wrappedPDFs[1] = getPDF(1);
}
const valarray<double>& TermData::getBinColumn(const string& n)
{
  const valarray<double>* ret = parent->getBinColumn(n);
  if (ret) return *ret;
  cerr << "[ERROR] Reaction \"" << reaction->getReactionName() << "\" requested bins of nonexistant column \"" << n << "\"" << endl;
  hf_errlog(19042700, "F: Reaction requested nonexistant bin column, see stderr");
  abort(); //unreachable
}
