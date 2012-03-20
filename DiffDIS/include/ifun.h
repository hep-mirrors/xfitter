/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _IFUN_H_
#define _IFUN_H_

#include <cstdio>
#include <cmath>

// oooooooooooooooooooooooooooooooooooooooooooooooooo
class ifun_t {
	virtual double f(double x) {return x*x;}
	double Gauss16(double a, double b, double eps);
	int gauss_error;
protected:
	double acc;
	ifun_t(double acc_=1e-5) {acc = acc_;}
public:
	void SetAcc(double acc_) {acc = acc_;}
	double GetAcc() {return acc;}
  double Integral(double a=0, double b=1, double eps=0) {
		return Gauss16(a,b,(eps <= 0 ? acc : eps));
	}
	bool isOK() {return !gauss_error;}
	bool isErr() {return gauss_error;}
};

#define NEW_IFUN(name,v)\
class name : public ifun_t {\
public: name(double acc_=1e-5) : ifun_t(acc_) {}\
private: double f(double v)

#endif
