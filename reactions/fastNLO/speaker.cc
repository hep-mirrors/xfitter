// Author: Daniel Britzger
// DESY, 06/08/2012

#include <iostream>
#include <map>
#include <string>
#include "fastnlotk/speaker.h"

using namespace std;

std::map<unsigned long,speaker*>* speaker::list = NULL;
std::ostream* speaker::weg = NULL;
// Set default logging level to start with; should be INFO.
// Use DEBUG only for debugging code executed before
// SetGlobalVerbosity(say::Verbosity volume)
// can be called from within the program.
say::Verbosity speaker::fverb = say::INFO;
unsigned long speaker::ct = 0;
bool speaker::fe2cerr = true;

speaker::speaker(std::string prefix,say::Verbosity volume,bool err,bool quiet) {
   if (list==NULL) list = new map<unsigned long,speaker*>();
   if (weg==NULL) {
      weg = new ostream(0);
      weg->clear(std::ios::badbit);
   }
   pref=prefix;
   fii=ct;
   (*list)[ct++] = this;
   fvol=volume;
   errs=err;
   fquiet= (quiet || fvol<fverb);
}


speaker::speaker(const speaker& other) :
   fquiet(other.fquiet), pref(other.pref),
   errs(other.errs), fvol(other.fvol),
   cn(other.cn)
{
   (*list)[ct++] = this;
}


speaker::~speaker() {
   list->erase(fii);
   if (list->empty()) {
      delete list;
      list=NULL;
      delete weg;
      weg=NULL;
   }
}

const speaker& speaker::operator=(const speaker& other)
{
   list->erase(fii);
   fii=ct;
   (*list)[ct++] = this;

   fquiet=other.fquiet;
   pref=other.pref;
   errs=other.errs;
   fvol=other.fvol;
   cn=other.cn;
   return *this;
}

std::ostream& speaker::operator()(const std::string& fct) const {
   if (fquiet) return *weg;
   //       *this<<"In "<<fct<<". ";
   if (errs && fe2cerr) return std::cerr<<fct;
   else return std::cout<<fct;
}

std::ostream& speaker::operator>> (const std::string& arg) const {
   return print(arg);
}

std::ostream& speaker::print(const std::string& mes) const {
   if (fquiet) return *weg;
   else {
      if (errs&&fe2cerr) return std::cerr<<mes;
      else  return std::cout<<mes;
   }
}

std::ostream& speaker::operator[](const std::string& fct) const {
   if (fquiet) return *weg;
   if (!cn.empty()) return *this<<"["<<cn<<"::"<<fct<<"] ";
   else return *this<<"["<<fct<<"] ";
}

const speaker& speaker::prefix(const std::string& fct) const {
   if (!fquiet) {
      if (errs&&fe2cerr) std::cerr<<fct;
      else std::cout<<fct;
   }
   return *this;
}

int speaker::SetGlobalVerbosity(say::Verbosity volume) {
   fverb=volume;
   int c=0;
   for (map<unsigned long, speaker*>::const_iterator ii=(*list).begin(); ii!=(*list).end(); ++ii) {
      (*ii).second->DoSpeak((*ii).second->GetVolume()>=volume);
      c++;
   }
   return c;
}

say::Verbosity speaker::GetGlobalVerbosity() {
   return fverb;
}

PrimalScream::PrimalScream(std::string classname) {
   debug = speaker(" # DEBUG.   ",say::DEBUG);
   man   = speaker(" # USAGE.   ",say::MANUAL);
   info  = speaker(" # INFO.    ",say::INFO);
   warn  = speaker(" # WARNING! ",say::WARNING);
   error = speaker(" # ERROR!   ",say::ERROR,true);
   debugsep = speaker("",say::DEBUG);
   mansep   = speaker("",say::MANUAL);
   infosep  = speaker("",say::INFO);
   warnsep  = speaker("",say::WARNING);
   errorsep = speaker("",say::ERROR,true);
   shout = speaker(" # ",say::ERROR,false);
   yell = speaker("",say::ERROR,false);
   SetClassName(classname);
}

void PrimalScream::SetClassName(const std::string classname){
   ___cn=classname;
   debug.SetClassName(___cn);
   man.SetClassName(___cn);
   info.SetClassName(___cn);
   warn.SetClassName(___cn);
   error.SetClassName(___cn);
   debugsep.SetClassName(___cn);
   mansep.SetClassName(___cn);
   infosep.SetClassName(___cn);
   warnsep.SetClassName(___cn);
   errorsep.SetClassName(___cn);
   shout.SetClassName(___cn);
   yell.SetClassName(___cn);
}

void PrimalScream::SetVerbosity(say::Verbosity volume) {
   debug.DoSpeak(debug.GetVolume() >= volume);
   man.DoSpeak(man.GetVolume() >= volume);
   info.DoSpeak(info.GetVolume() >= volume);
   warn.DoSpeak(warn.GetVolume() >= volume);
   error.DoSpeak(error.GetVolume() >= volume);
   debugsep.DoSpeak(debug.GetVolume() >= volume);
   mansep.DoSpeak(man.GetVolume() >= volume);
   infosep.DoSpeak(info.GetVolume() >= volume);
   warnsep.DoSpeak(warn.GetVolume() >= volume);
   errorsep.DoSpeak(error.GetVolume() >= volume);
   shout.DoSpeak(shout.GetVolume() >= volume);
   yell.DoSpeak(yell.GetVolume() >= volume);
}

namespace say {
   speaker debug(" # DEBUG.   ",say::DEBUG);
   speaker man  (" # USAGE.   ",say::MANUAL);
   speaker info (" # INFO.    ",say::INFO);
   speaker warn (" # WARNING! ",say::WARNING);
   speaker error(" # ERROR!   ",say::ERROR,true);
   speaker debugsep("",say::DEBUG);
   speaker mansep  ("",say::MANUAL);
   speaker infosep ("",say::INFO);
   speaker warnsep ("",say::WARNING);
   speaker errorsep("",say::ERROR,true);
   speaker shout(" # ",say::ERROR,false);
   speaker yell("",say::ERROR,false);
   int SetGlobalVerbosity(Verbosity verbosity) {
      return speaker::SetGlobalVerbosity(verbosity);
   };
   Verbosity GetGlobalVerbosity() {
      return speaker::fverb;
   }
}
