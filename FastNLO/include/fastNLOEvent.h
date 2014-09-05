#ifndef __fnloevent__
#define __fnloevent__

#include <string>
#include <vector>
#include <map>

class fnloScenario {
   //! useful class to keep all scenario specific quantities.
   //! e.g. observables and scales
   friend class fastNLOCreate;

 public:
 fnloScenario() : _iOB(-1){;}
   ~fnloScenario() {;};

   void SetObservableDimI(double o, int iDim) {_o[iDim]=o;}                                     //!< Set observable of dimension iDim (e.g. in case of multidimensional measurements)
   void SetObservable0(double o) {SetObservableDimI(o,0);}                                      //!< Set observable for '0th' dimension for single-differential calculation
   void SetObservable1(double o) {SetObservableDimI(o,1);}                                      //!< Set observable for '1st' dimension for single and double-differential calculations
   void SetObservable2(double o) {SetObservableDimI(o,2);}                                      //!< Set observable for '2nd' dimension for single/double/triple differential calculations
   void SetObsBin(int iBin) {_iOB = iBin; }                                                     //!< [optional] Set ObsBin (e.g. if binning is performed by generator, no other observables are then needed.)
   //! flexible scale table:
   void SetObsScale1(double mu) {_m1=mu;}                                                       //!< For flexible-scale tables. Set scale 1 (should be in 'GeV').

   void SetObsScale2(double mu) {_m2=mu;}                                                       //!< For flexible-scale tables. Set scale 2
   //! if not a flexible-scale table
   void SetScale(double mu) {_m1=mu;}                                                           //!< the ren. and fact. scale (not mu^2)
private:
   std::map<int,double> _o;
   double _m1, _m2;
   int _iOB;
};


class fnloEvent {
   //! useful class to keep all process related variables.
   //! e.g x-values, weights, process identifiers, etc...
   friend class fastNLOCreate;

public:
   fnloEvent(){Reset();}
   ~fnloEvent(){;}
   void ResetButX(){
      _w=0,_wf=0,_wr=0,_wrr=0,_wff=0,_wrf=0;
      _p=-1;
      _n=-1;
   }
   void Reset(){
      ResetButX();
      _x1=0,_x2=0;
   }
   // event specific quantites, which are required for every 'Fill()' step.
   void SetX(double x) {_x1=x;}                                                                 //!< set x-value of first hadron (if e.g. DIS)
   void SetX1(double x) {_x1=x;}                                                                //!< setx-value of first hadron
   void SetX2(double x) {_x2=x;}                                                                //!< set x-value of second hadron
   void SetProcessId(int n){_p=n;}                                                              //!< set identifier of specific subprocess (0<n<NSubproc), according to the corresponding PDF linear combination
   void SetEventCounter(long long int n){_n=n;}                                                 //!< Set event counter
   //! if not a flexible-scale table
   void SetWeight(double w) {_w=w;}                                                             //!< weights must be mutliplied with dummypdf (1/x)
   //! flexible scale table:
   void SetWeight_MuIndependent(double w) {_w=w;}                                               //!< weights must be mutliplied with dummypdf (1/x)
   void SetWeight_log_mur(double w) {_wr=w;}                                                    //!< set weight w, which will contribute with log_e(mur^2)*w
   void SetWeight_log_muf(double w) {_wf=w;}                                                    //!< set weight w, which will contribute with log_e(muf^2)*w
   void SetWeight_log_murr(double w) {_wrr=w;}                                                  //!< set weight w, which will contribute with log^2_e(mur^2)*w
   void SetWeight_log_muff(double w) {_wff=w;}                                                  //!< set weight w, which will contribute with log^2_e(muf^2)*w
   void SetWeight_log_murf(double w) {_wrf=w;}                                                  //!< set weight w, which will contribute with log_e(mur^2)*log_e(muf^2)*w
   void AddWeight_MuIndependent(double w) {_w+=w;}                                              //!< weights must be mutliplied with dummypdf (1/x)
   void AddWeight_log_mur(double w) {_wr+=w;}                                                   //!< set weight w, which will contribute with log_e(mur^2)*w
   void AddWeight_log_muf(double w) {_wf+=w;}                                                   //!< set weight w, which will contribute with log_e(muf^2)*w
   void AddWeight_log_murr(double w) {_wrr+=w;}                                                 //!< set weight w, which will contribute with log^2_e(mur^2)*w
   void AddWeight_log_muff(double w) {_wff+=w;}                                                 //!< set weight w, which will contribute with log^2_e(muf^2)*w
   void AddWeight_log_murf(double w) {_wrf+=w;}                                                 //!< set weight w, which will contribute with log_e(mur^2)*log_e(muf^2)*w
private:
   double _x1, _x2;                                                                             //!< an event has always identical x1 and x2;
   double _w, _wf, _wr, _wrr, _wff, _wrf;                                                       //!< weights
   int _p;                                                                                      //!< processId/channel. Must be consistent with PDF linear combination
   long long int _n;                                                                            //!< event count
};

#endif
