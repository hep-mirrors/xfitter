// Author: Daniel Britzger
// DESY, 06/08/2012

#ifndef SPEAKER_H_
#define SPEAKER_H_

#include <string>
#include <iostream>
#include <map>

namespace say {
   enum Verbosity {DEBUG=-1000, MANUAL=2, INFO=0, WARNING=1, ERROR=2, SILENT=1000};

   // String to enum
   struct toVerbosity : public std::map<std::string, Verbosity>
      {
         toVerbosity()
            {
               this->operator[]("DEBUG")   = DEBUG;
               this->operator[]("MANUAL")  = MANUAL;
               this->operator[]("INFO")    = INFO;
               this->operator[]("WARNING") = WARNING;
               this->operator[]("ERROR")   = ERROR;
               this->operator[]("SILENT")  = SILENT;
            };
         ~toVerbosity(){}
      };

}

class speaker {
public:
   speaker(std::string prefix="",say::Verbosity volume=say::INFO,bool err=false,bool quiet=false);
   speaker(const speaker& spk);
   ~speaker();
   const speaker& operator= (const speaker& other);
   std::ostream& operator[](const std::string& fct) const ;
   const speaker& operator+ (const std::string& fct) const {
      return this->prefix(fct);
   }
   const speaker& prefix(const std::string& fct) const ;
   std::ostream& operator()(const std::string& fct) const ;

   template<typename T> std::ostream& operator<< (const T& arg) const {
      if (fquiet) return *weg;
      else {
         if (errs && fe2cerr) return std::cerr<<pref<<arg;
         else return std::cout<<pref<<arg;
      }
   }
#ifndef SWIG
   std::ostream& operator>> (const std::string& arg) const ;
#endif
   std::ostream& print(const std::string& mes) const ;
   void DoSpeak(bool loud) {
      fquiet=!loud;
   };
   bool GetSpeak() const {
      return !fquiet;
   };
   void SetPrefix(std::string prefix) {
      pref=prefix;
   };
   std::string GetPrefix() const {
      return pref;
   };
   void SetClassName(std::string classname) {
      cn=classname;
   };
   std::string GetClassName(void) const {
      return cn;
   };
   say::Verbosity GetVolume(void) const {
      return fvol;
   };
   void SetVolume(say::Verbosity volume) {
      fvol=volume;
   };
   static int SetGlobalVerbosity(say::Verbosity volume);
   say::Verbosity GetGlobalVerbosity();
   static say::Verbosity fverb;
   static void ErrorToErrStream(bool ToCerr) {
      fe2cerr=ToCerr;
   };

protected:
   static std::ostream* weg;
   bool fquiet;
   std::string pref;
   bool errs;
   say::Verbosity fvol;
   unsigned long fii;
   static unsigned long ct;
   static bool fe2cerr;
   static std::map<unsigned long,speaker*>* list;
   std::string cn;
};

namespace say {
   extern speaker debug;
   extern speaker man;
   extern speaker info;
   extern speaker warn;
   extern speaker error;
   extern speaker debugsep;
   extern speaker mansep;
   extern speaker infosep;
   extern speaker warnsep;
   extern speaker errorsep;
   extern speaker shout; // same as error but streamed to cout
   extern speaker yell;  // same as error but streamed to cout without prefix
   extern int SetGlobalVerbosity(Verbosity verbosity);
   extern Verbosity GetGlobalVerbosity(void);
}


class PrimalScream {
public:
   PrimalScream(std::string classname);
   void SetClassName(const std::string classname );
   void SetVerbosity(say::Verbosity volume);
   speaker debug;
   speaker man;
   speaker info;
   speaker warn;
   speaker error;
   speaker debugsep;
   speaker mansep;
   speaker infosep;
   speaker warnsep;
   speaker errorsep;
   speaker shout;
   speaker yell;
private:
   std::string ___cn;
};

#endif //SPEAKER_H_
