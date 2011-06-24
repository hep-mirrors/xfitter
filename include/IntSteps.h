#ifndef IntSteps_h
#define IntSteps_h 1

#include <string>

class IntSteps
{
 public:
  IntSteps(){};
  ~IntSteps();

  IntSteps(const IntSteps&);
  IntSteps(const std::string, const double *, const double *, const double);

 protected:
  std::string _boz;

  double *_mr;
  double *_yr;

  int _nms;
  int _nys;
  double *_msteps;
  double *_ysteps;

  double _leptPtCut;

 private:
  void _makeMstepsZ();
  void _makeMstepsW();
  void _makeYsteps();

 public:
  std::string getBozName(){ return _boz;}

  void getDims(int &nms, int &nys){
    nms = _nms;
    nys = _nys;
  }

  void getMsteps(int &nms, double *m_steps){
    nms=_nms;
    m_steps=_msteps;
  }

  void getYsteps(int &nys, double *y_steps){
    nys=_nys;
    y_steps=_ysteps;
  }

  // dirty lepton pt cut storage
  void setLeptPtCut(const double leptPtCut){
    _leptPtCut = leptPtCut;
  }

  double getLeptPtCut() { return _leptPtCut; }
};


#endif
