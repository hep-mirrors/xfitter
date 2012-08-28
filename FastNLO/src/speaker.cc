// Author: Daniel Britzger
// DESY, 06/08/2012

#include "speaker.h"
#include <string>
#include <iostream>
#include <map>


using namespace std;


std::map<unsigned long,speaker*>* speaker::list = NULL;
std::ostream* speaker::weg = NULL;
say::Verbosity speaker::fverb = say::INFO;
unsigned long speaker::ct = 0;
bool speaker::fe2cerr = true;

speaker::speaker(std::string prefix,say::Verbosity volume,bool err,bool quiet){
   if ( list==NULL ) list = new map<unsigned long,speaker*>();
   if ( weg==NULL ) {
      weg = new ostream(0);
      weg->clear(std::ios::badbit);
   }
   pref=prefix;
   fii=ct;
   (*list)[ct++] = this;
   fvol=volume;
   errs=err;
   fquiet= ( quiet || fvol<fverb );
}

speaker::~speaker(){
   list->erase(fii);
   if(list->empty()){
      delete list; list=NULL;
      delete weg; weg=NULL;
   }
}

std::ostream& speaker::operator() (std::string fct) const {
   if (fquiet) return *weg;
   //       *this<<"In "<<fct<<". ";
   if(errs && fe2cerr) return std::cerr<<fct;
   else return std::cout<<fct;
}

std::ostream& speaker::operator>> (std::string arg) const {
   return print(arg);
}

std::ostream& speaker::print(std::string mes) const {
   if (fquiet) return *weg;
   else {	
      if (errs&&fe2cerr) return std::cerr<<mes;
      else  return std::cout<<mes;
   }
}

std::ostream& speaker::operator[] (std::string fct) const {
   if ( fquiet ) return *weg;
   if ( !cn.empty()) return *this<<"["<<cn<<"::"<<fct<<"] ";
   else return *this<<"["<<fct<<"] ";
}

const speaker& speaker::prefix(std::string fct) const {
   if ( !fquiet ) {
      if (errs&&fe2cerr) std::cerr<<fct;
      else std::cout<<fct;
   }
   return *this; 
}


int speaker::SetGlobalVerbosity(say::Verbosity volume){
   fverb=volume;
   int c=0;
   for( map<unsigned long, speaker*>::const_iterator ii=(*list).begin(); ii!=(*list).end(); ++ii){
      (*ii).second->DoSpeak( (*ii).second->GetVolume()>=volume );
      c++;
   }
   return c;
}


PrimalScream::PrimalScream(std::string classname){//,std::string prefix=""){
   cn=classname;
   debug = speaker("Debug. ",say::DEBUG);
   debug.SetClassName(cn);
   man   = speaker("",say::MANUAL);
   man.SetClassName(cn);
   info  = speaker("Info. ",say::INFO);
   info.SetClassName(cn);
   warn  = speaker("Warning. ",say::WARNING);
   warn.SetClassName(cn);
   error = speaker("Error! ",say::ERROR,true);
   error.SetClassName(cn);
   debug["PrimalScream"]<<"Primal Scream initialized."<<std::endl;
}

void PrimalScream::SetVerbosity(say::Verbosity volume){
   debug.DoSpeak( debug.GetVolume() >= volume );
   man.DoSpeak( man.GetVolume() >= volume );
   info.DoSpeak( info.GetVolume() >= volume );
   warn.DoSpeak( warn.GetVolume() >= volume );
   error.DoSpeak( error.GetVolume() >= volume );
}

namespace say {
   speaker debug("Debug. ",say::DEBUG);
   speaker man("",say::MANUAL);
   speaker info("Info. ",say::INFO);
   speaker warn("Warning. ",say::WARNING);
   speaker error("Error! ",say::ERROR,true);
   //debug["namespace say"]<<"speakers initialized."<<std::endl;
   int SetGlobalVerbosity(Verbosity verbosity){ return speaker::SetGlobalVerbosity(verbosity);};
}

