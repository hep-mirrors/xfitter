#include "sock_DDIS.h"
#include "qcdnum_pdf.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FORTRAN export
*/

extern "C" {
  double ddisvalue_(double *kinvars);
  void ddisinit_();
  void setecmsq_(const double *s);
  void ddissetparams_(const double*);
}

static sock_DDIS_t* g_sock_DDIS;
static qcdnum_pdf_t* g_QCDNUM_pdf;
static PhysParams_t g_PhysPars(0,1);  //--- needed omly to make Initialize quiet (???)

// =========================================
// static void g_DDISInit() {
void ddisinit_() {
  delete g_sock_DDIS;
  delete g_QCDNUM_pdf;
  g_sock_DDIS = new sock_DDIS_t;
  g_QCDNUM_pdf = new qcdnum_pdf_t;
  string ccf("DDIS.coca");
  if(FileReadable(ccf)) g_sock_DDIS->Params.Read(ccf.c_str());
  g_sock_DDIS->Params.Show();
  g_sock_DDIS->SetPhys(&g_PhysPars);
  g_sock_DDIS->SetPDF(g_QCDNUM_pdf); 
  g_sock_DDIS->Initialize();
  // g_sock_DDIS->SetEcmSquared(4*27.5*920);
}

// =========================================
void setecmsq_(const double *s) {
  g_sock_DDIS->SetEcmSquared(*s);
}

// =========================================
double ddisvalue_(double *kinvars) {
  // --- needs ECM fixed or y specified for each data record
  //cout << "ddisvalue_: " << kinvars[0] << " "  << kinvars[1] << " "  << kinvars[2]
  //     <<"   "<< g_sock_DDIS->Value(kinvars) << endl;
  return g_sock_DDIS->Value(kinvars);
}
  
/**
  fitpars[3]:
  \arg Pomeron_a0 — Pomeron intercept
  \arg Reggeon_factor
  \arg Reggeon_a0 — Reggeon intercept
*/
// =========================================
void ddissetparams_(const double* fitpars) {
  g_sock_DDIS->SetParams(fitpars);
}
