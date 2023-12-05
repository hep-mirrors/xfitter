/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _FIASCO_H_
#define _FIASCO_H_

// #include <iostream>
#include <cstdlib>
#include <cstdio>
// #include <cstring>
#include <cstdarg>
#include <stdexcept>

using namespace std;

/// Exception class
// ooooooooooooooooooooooooooooooo
class Fiasco {
  char msg[1024];
public:
  Fiasco(const char* fmt,...) {
	  va_list apt;
		va_start(apt,fmt);
		vsnprintf(msg, sizeof(msg), fmt, apt);
		va_end(apt); 
	}
	
	string what() const { return string(msg);}
	// const char* Msg() { return msg;}
};

#define START_XCODE try {
#define END_XCODE } \
  catch(const logic_error& x) {cerr <<"logic_error: " << x.what() << endl; } \
  catch(const runtime_error& x) {cerr <<"runtime_error: " << x.what() << endl; } \
  catch(const Fiasco& f) {cerr << f.what() << endl;}\
  catch(const char* msg) {cerr << msg << endl;} \
  cerr << "\nFATAL error - program terminated." << endl; \
  exit(1);

#endif
