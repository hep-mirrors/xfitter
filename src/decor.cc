/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "decor.h"

struct cmp_vec_elems_t {
  Vector v;
  bool Ascend;
  bool operator() (int i,int j) {
    if(Ascend) return v[i] < v[j];
    else return v[i] > v[j];
  }
};
static cmp_vec_elems_t cmp_vec_elems;

//===========================
int Decor_t::SetCov(const char* fn) {
  std::ifstream dat(fn);
  int nn;
  
  #ifdef DEBUG_OFFSET
    cout << "SetCov opening: " << fn << endl;
  #endif
  if(!dat) return 1;
  dat >> N;
  DBG_SHOW(N)
  double v[nn=N*N], *vp = v;
  while(nn--) dat >> *vp++;
  if(!dat) {dat.close(); return 2;}
  dat.close();
  SetCov(N, v);
  return 0;
}

//===========================
void Decor_t::Sort(bool asc) {
  vector<int> Index(N);
  int k;
  for(k=0; k < N; k++) Index[k] = k;
  cmp_vec_elems.v = EigVals;
  cmp_vec_elems.Ascend = asc;
  stable_sort(Index.begin(), Index.end(), cmp_vec_elems);
  Vector u(EigVals);
  SqMatrix_t A(EigVecs);
  // for(k=0; k < N; k++) cout << k <<": "<< Index[k] <<"  "<< u[Index[k]] << endl;
  for(k=0; k < N; k++) {
    EigVals[k] = u[Index[k]];
    EigVecs.SetCol(k, A.GetCol(Index[k]));
    EigShifts[k] = sqrt(DeltaChiSqr*EigVals[k]);
  }
  // for(k=0; k < N; k++) cout << k <<": "<< EigVals[k] << endl;
}
