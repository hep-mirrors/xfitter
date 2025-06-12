#ifndef __fnloevent__
#define __fnloevent__

#include <string>
#include <vector>
#include <map>

class fnloScenario {
   //! Useful class to keep all scenario specific quantities,
   //! e.g. observables and scales.
   //! KR: Test including additional weight factor per observable entry
   //!     This could be helpful e.g. for jet cross section differences!
   friend class fastNLOCreate;

 public:
 fnloScenario() : _m1(0), _m2(0), _iOB(-1), _wo(1) {;}
   ~fnloScenario() {;};
   //! Store observable
   //!< Set observable of dimension iDim (e.g. in case of multidimensional measurements)
   inline void SetObservableDimI(double o, int iDim) {_o[iDim]=o;}
   //!< Set observable for '0th' dimension for single-differential calculation
   inline void SetObservable0(double o) {SetObservableDimI(o,0);}
   //!< Set observable for '1st' dimension for single and double-differential calculations
   inline void SetObservable1(double o) {SetObservableDimI(o,1);}
   //!< Set observable for '2nd' dimension for single/double/triple differential calculations
   inline void SetObservable2(double o) {SetObservableDimI(o,2);}
   //!< [ALTERNATIVELY] Directly set ObsBin of serialised one-dimensional representation,
   //!< e.g. when serialisation is performed already in generator.
   //!< The observables per dimension then are not needed, but take care with GetBin!
   inline void SetObsBin(int iBin) {_iOB = iBin; }

   //! Store central scale choice for mur and muf (not mu^2!)
   //!< Fixed-scale table (scale must be given in units of GeV)
   inline void SetScale(double mu) {_m1=mu;}
   //!< Flexible-scale table
   //!< Set scale 1 (scale must be given in units of GeV)
   inline void SetObsScale1(double mu) {_m1=mu;}
   //!< Set scale 2
   inline void SetObsScale2(double mu) {_m2=mu;}

   //! Store extra weight factor
   inline void SetObsWeight(double w) {_wo=w;}

 private:
   std::map<int,double> _o;
   double _m1, _m2;
   int _iOB;
   double _wo;
};



class fnloEvent {
   //! Useful class to keep all process related variables.
   //! e.g x-values, weights, process identifiers, etc.
   friend class fastNLOCreate;
   friend class fastNLOCoeffAddFix;
   friend class fastNLOCoeffAddFlex;

public:
   fnloEvent(){Reset();}
   ~fnloEvent(){;}

   inline void ResetButX(){
      _w=0, _wf=0, _wr=0, _wrr=0, _wff=0, _wrf=0;
      _sig = 0; // sigma
      _p = -1;
      _n = -1;
   }
   inline void Reset(){
      ResetButX();
      _x1 = 0, _x2 = 0;
   }

   // Event specific quantites, which are required for every 'Fill()' step.
   //!< Get x-value of first hadron
   inline const double& X1() const {return _x1;}
   //!< Get x-value of second hadron
   inline const double& X2() const {return _x2;}
   //!< Get subprocess ID
   inline const int& p() const {return _p;}
   //!< Set x-value of first hadron (if e.g. DIS)
   inline void SetX(double x) {_x1=x;}
   //!< Set x-value of first hadron
   inline void SetX1(double x) {_x1=x;}
   //!< Set x-value of second hadron
   inline void SetX2(double x) {_x2=x;}
   //!< Set subprocess ID
   inline void SetProcessId(int n){_p=n;}
   //!< Set event counter
   inline void SetEventCounter(long long int n){_n=n;}

   //! Fixed-scale table:
   //!< Weights must be multiplied with dummypdf (1/x)
   inline void SetWeight(double w) {_w=w;}
   //!< Set weight to calculate cross section (i.e. already multiplied by PDF,alpha_s).
   inline void SetSigma(double s) {_sig=s;}
   //!< Add weight
   inline void AddSigma(double s) {_sig+=s;}

   //! Flexible-scale table:
   //!< Weights must be multiplied with dummypdf (1/x)
   inline void SetWeight_MuIndependent(double w) {_w=w;}
   //!< Set weight w, which will contribute with log_e(mur^2)*w
   inline void SetWeight_log_mur(double w) {_wr=w;}
   //!< Set weight w, which will contribute with log_e(muf^2)*w
   inline void SetWeight_log_muf(double w) {_wf=w;}
   //!< Set weight w, which will contribute with log^2_e(mur^2)*w
   inline void SetWeight_log_murr(double w) {_wrr=w;}
   //!< Set weight w, which will contribute with log^2_e(muf^2)*w
   inline void SetWeight_log_muff(double w) {_wff=w;}
   //!< Set weight w, which will contribute with log_e(mur^2)*log_e(muf^2)*w
   inline void SetWeight_log_murf(double w) {_wrf=w;}
   //!< Weights must be multiplied with dummypdf (1/x)
   inline void AddWeight_MuIndependent(double w) {_w+=w;}
   //!< Add weight w, which will contribute with log_e(mur^2)*w
   inline void AddWeight_log_mur(double w) {_wr+=w;}
   //!< Add weight w, which will contribute with log_e(muf^2)*w
   inline void AddWeight_log_muf(double w) {_wf+=w;}
   //!< Add weight w, which will contribute with log^2_e(mur^2)*w
   inline void AddWeight_log_murr(double w) {_wrr+=w;}
   //!< Add weight w, which will contribute with log^2_e(muf^2)*w
   inline void AddWeight_log_muff(double w) {_wff+=w;}
   //!< Add weight w, which will contribute with log_e(mur^2)*log_e(muf^2)*w
   inline void AddWeight_log_murf(double w) {_wrf+=w;}

 private:
   //!< An event has always identical x1 and x2
   double _x1, _x2;
   //!< Sigma, i.e. weight including PDF & alpha_s
   double _sig;
   //!< Weights
   double _w, _wf, _wr, _wrr, _wff, _wrf;
   //!< SubprocessId/channel. Must be consistent with PDF linear combination
   int _p;
   //!< Event count
   long long int _n;
};

#endif
