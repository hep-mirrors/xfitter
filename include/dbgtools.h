/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef DBG_TOOLS_H_
#define DBG_TOOLS_H_

// #undef DBG_SHOW
// #undef DBG_HERE
#ifdef DBGT_ON
  #include <iostream>
  using namespace std;
  #define DBG_SHOW(a) cout <<"*** "<< __FILE__ <<"("<< __LINE__<<"): " << #a" = '" << (a) <<"'"<< endl;
  #define DBG_HERE(msg) cout <<"*** "<< __FILE__ <<"("<< __LINE__<<"): " << msg << endl;
#else
  #define DBG_SHOW(a)
  #define DBG_HERE(msg)
#endif

#endif
