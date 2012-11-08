/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string.h>
using namespace std;

class Vector;

//ooooooooooooooooooooooooooooooooooooooooo
class SqMatrix_t {
public:
  enum {RW_FULL, RW_LOLEFT, RW_UPRIGTH} RWmode;
  int first,last,len;
  double** mp;

  //===========================
  void M_aloc() {
    mp = new double*[len] - first;
    for(int row=first; row <= last; row++) {
      mp[row] = new double[len] - first;
      memset(mp[row]+first,0,len*sizeof(double));
    }
  }

  //===========================
  void M_free() {
    if(!len) return;
    for(int row=first; row <= last; row++) delete[] (mp[row]+first);
    delete[] (mp+first);
    len = 0;
  }

  //===========================
  SqMatrix_t(int lo=0, int n=0) {
    first = lo;
    last = lo+n-1;
    len = n;
    if(len) M_aloc();
    RWmode = RW_FULL;
  }
  //SqMatrix_t(int hi) { SqMatrix_t(1,hi);}

  //========================
  SqMatrix_t& Zero() {
    for(int row=first; row <= last; row++)
      memset(mp[row]+first,0,len*sizeof(double));
    return *this;
  }

  //========================
  SqMatrix_t& operator= (const SqMatrix_t& A) {
    if(this == &A) return *this;
    if(len) M_free();
    first = A.first;
    last = A.last;
    len = A.len;
    if(len) {
      M_aloc();
      Set(A);
    }
    return *this;
  }

  //=============================
  SqMatrix_t(const SqMatrix_t& A) {
    //if(len) M_free();
    first = A.first;
    last = A.last;
    len = A.len;
    if(!len) return;
    M_aloc();
    Set(A);
    /*
    int k;
    for(int row=first; row <= last; row++) {
      for(k=first; k <= last; k++) mp[row][k]=A[row][k];
    }
    */
    RWmode = RW_FULL;
  }

  //=============================
  ~SqMatrix_t(void) {
    M_free();
  }
  
  //========================
  void StartIndex(int f) {
    if(f == first) return;
    int df = first -f;
    for(int row=first; row <= last; row++) mp[row] += df;
    mp += df;
    first = f;
    last = first+len-1;
  }

  //===========================
  void Set(const SqMatrix_t& A) {
    if(A.len != len) {
      cerr << "SqMatrix_t::Set: diff. lengths"<<endl;
      exit(1);
    }
    for(int row=first; row <= last; row++) {
      //for(int k=first; k <= last; k++) mp[row][k]=A[row][k];
      memcpy(mp[row]+first, A[row-first + A.first]+A.first, len*sizeof(double));
    }
  }

  //==========================
  SqMatrix_t& Transpose() {
    double t;
    for(int row=first; row <= last; row++) {
      for(int k=first; k < row; k++) {
        t = mp[row][k];
        mp[row][k] = mp[k][row];
        mp[k][row] = t;
      }
    }
    return *this;
  }

  //==========================
  void Set(double* v) {
    for(int row=first; row <= last; row++) {
      for(int k=first; k <= last; k++) mp[row][k] = *v++;
    }
  }

  //==========================
  void Get(double* v) {
    for(int row=first; row <= last; row++) {
      for(int k=first; k <= last; k++) *v++ = mp[row][k];
    }
  }

  //==========================
  SqMatrix_t& Inverse(); //--- gaussj
  int Diagonalize(Vector& eigenvals);
  SqMatrix_t& Rmul(const SqMatrix_t& B);
  void GetCol(int k, Vector& c);
  void SetCol(int k, Vector& c);
  Vector& GetCol(int k);
  double MxElem(Vector& v);

  //===========================
  double* operator[](int r) {return mp[r];}
  double* operator[](int r) const {return mp[r];}

  //====================================
  void Write(ostream &ostr=cout) const {
    int prec = 6;
    ios_base::fmtflags ff = ostr.flags();
    streamsize w = ostr.width();
    streamsize p = ostr.precision();
    ostr.precision(prec);
    for(int j = first; j <= last; j++) {
      // ostr << j <<":";
      int k0 = (RWmode == RW_UPRIGTH)? j : first;
      int k1 = (RWmode == RW_LOLEFT)? j : last;
      for(int k = k0; k <= k1; k++)
        // ostr << " "<< setw(13) << mp[j][k];
        ostr << scientific << setw(prec+9) << right << mp[j][k];
        // ostr << fixed << setw(prec+9) << right << mp[j][k];
      ostr << endl;
    }
    ostr.flags(ff);
    ostr.width(w);
    // ostr.unsetf(ios::floatfield);
    ostr.precision(p);
  }

  // //====================================
  // friend ostream& operator<<(ostream &ostr, const SqMatrix_t& a) {
    // a.Write(ostr);
    // return ostr << endl;
  // }

};

//====================================
ostream& operator<<(ostream &ostr, const SqMatrix_t& a);

//ooooooooooooooooooooooooooooooooooooooooo
class Vector {
public:
  int first,last,len;
  double* vp;

  //========================
  Vector(int lo=0, int n=0) {
    first = lo;
    last = lo+n-1;
    len = n;
    if(len) {
      vp = new double[len];
      memset(vp, 0, len*sizeof(double));
      vp -= first;
    }
  }

  //========================
  Vector& Set(double *x) {
    memcpy(vp+first, x, len*sizeof(double));
    return *this;
  }

  //========================
  Vector& Zero() {
    memset(vp+first, 0, len*sizeof(double));
    return *this;
  }

  //===========================
  void Set(Vector& v) {
    if(v.len != len) {
      cerr << "Vector::Set: diff. lengths"<<endl;
      exit(1);
    }
    memcpy(vp+first, v.vp+v.first, len*sizeof(double));
  }

  //========================
  Vector(Vector& v) {
    first = v.first;
    last = v.last;
    len = v.len;
    vp = new double[len];
    vp -= first;
    Set(v.vp+v.first);
  }

  //========================
  ~Vector(void) {
    if(len) delete[] (vp+first);
    len = 0;
  }

  //========================
  void StartIndex(int f) {
    if(f == first) return;
    vp += first -f;
    first = f;
    last = first+len-1;
  }

  //========================
  double& operator[](int index) {return vp[index];}
  double operator[](int index) const {return vp[index];}

  /*
  //========================
  Vector& operator- () {
    Vector V(*this);
    return V.Scale(-1.);
  }
  */

  //========================
  Vector& operator= (const Vector& v) {
    if(this == &v) return *this;
    if(len) delete[] (vp+first);
    first = v.first;
    last = v.last;
    len = v.len;
    vp = new double[len];
    vp -= first;
    Set(v.vp+v.first);
    return *this;
  }

  //========================
  void Lmul(SqMatrix_t& MM);
  //========================
  double ScalarProduct(Vector& v);

  //========================
  Vector& Scale(double s) {
    for(int r=first; r <= last; r++) vp[r] *= s;
    return *this;
  }

  //========================
  double Length() {
    double s=0;
    for(int r=first; r <= last; r++) s += vp[r]*vp[r];
    return sqrt(s);
  }

  //========================
  void print() {
    int j;
    for(j=first; j <= last; j++) {
      cout << j <<": "<< vp[j] << endl;
    }
  }

  //========================
  friend ostream &operator<<(ostream &ostr, const Vector& a) {
    for(int j=a.first; j <= a.last; j++)
    ostr << j <<": "<< a.vp[j] << endl;
    return ostr << endl;
  }

};

  SqMatrix_t operator+(const SqMatrix_t& A, const SqMatrix_t& B);

#endif
