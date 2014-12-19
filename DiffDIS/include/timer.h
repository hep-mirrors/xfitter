/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef TIMER_H_
#define TIMER_H_

#if __cplusplus < 201103L
#define constexpr 
#endif



#include <ctime>
// #include <iostream>
#include <sstream>
using namespace std;

// ooooooooooooooooooooooooooooooooo
class Timer_t {
  constexpr static const double sec = 1./CLOCKS_PER_SEC;
  clock_t M_start;
  
public:
  // ===========
  Timer_t() {
    M_start = clock();
  }
  
  // ===========
  void Reset() {
    M_start = clock();
  }
  
  // =========================
  double Seconds() const {
    return (clock() - M_start)*sec;
  }

  // =========================
  string Value() const {
    stringstream vs(ios_base::out);
    int cs = 100*Seconds();
    int h = cs / 360000;
    int m = cs % 360000;
    cs = m % 6000;
    m /= 6000;
    if(h) vs << h << "h ";
    if(h || m) vs << m << "m ";
    vs << cs*0.01 << "s";
    return vs.str();
  }

};

#endif
