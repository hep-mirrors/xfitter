/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 4.1.3
_____________________________________________________________*/

#ifndef CLS_THEMODEL_BASE_H_
#define CLS_THEMODEL_BASE_H_

// #include "genut.h"
// #include "dbgtools.h"
#include "Coca.h"
#include "pdf_base.h"

// ooooooooooooooooooooooooooooooooo
class TheModel_base_t {
protected:
  bool isInitialized;
  pdf_base_t* Pdf;
  PhysParams_t* Phys;
  string M_Name;
  
public:
  Coca_t Params;
  
  // ===========================
  TheModel_base_t(const string& name="") {
    isInitialized = false;
    Phys = NULL;
    Pdf = NULL;
    M_Name = name;
  }
  
  // ===========================
  string GetName() const {
    if(M_Name.empty()) return Params.GetName();
    return M_Name;
  }
  
  void SetPDF(pdf_base_t* pdf) {Pdf = pdf;}
  void SetPhys(PhysParams_t* phys) {Phys = phys;}
  
  /**
    \brief Read params from \c ExpData. 
    
    This method is called upon connecting to a \c ExpData_t object, i.e.
    by \c ExpData_t::SetModel() which supplies \c CC of the calling \c ExpData_t object.
  */
  virtual void GetExpParams(const Coca_t& EDcoca) {}
  // virtual void GetExpParams(const Coca_t* ExpCC) {}
  
  virtual double Value(double* x) = 0;
  
  virtual void SetParams(const vector<double>& ) {}
  virtual void SetParams(const double* ) {}
  /**
    Before calling \c Initialize or \c Configure
    \c Pdf and \c Phys must be set 
  */
  virtual void Initialize() {
    // cout << "==> TheModel_base_t::Initialize" << endl;
    isInitialized = true;
    Configure();
  }
  
  void Initialize(const Coca_t& cc) {
    Params.Set(cc);
    Initialize();
  }
  
  void Initialize(const string& ccfn) {
    Params.Read(ccfn.c_str()); 
    Initialize();
  }
  
  virtual void Configure() {}
  
  void Configure(const Coca_t& cc) {
    Params.Set(cc);
    Configure();
  }
  
  void Configure(const string& ccfn) {
    Params.Read(ccfn.c_str()); 
    Configure();
  }
  
  virtual void Configure(const vector<double>& ) {}
  virtual void Configure(const double* params, int n=-1) {}
  
  // virtual void Update(const Coca_t& coca) {Configure(coca);}
  // virtual void Update() {Configure();}
};

#endif
