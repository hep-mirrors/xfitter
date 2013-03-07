#include "sock_DDIS.h"
#include "qcdnum_pdf.h"
#include <map>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FORTRAN export
*/

extern "C" {
  double ddisvalue_(const int* IDataSet, double *kinvars);
  void ddisinit_(const int* IDataSet);
  void setecmsq_(const int* IDataSet, const double *s);
  void ddissetparams_(const double*);
}

// static sock_DDIS_t* g_sock_DDIS;
// static map<int,sock_DDIS_t*> g_sock_DDIS;
typedef map<int,sock_DDIS_t*> dict_t;
static dict_t g_sock_DDIS;
static qcdnum_pdf_t* g_QCDNUM_pdf;
static PhysParams_t g_PhysPars(0,1);  //--- needed omly to make Initialize quiet (???)

// =========================================
void ddisinit_(const int* IDataSet) {
  int DataSetIndex = *IDataSet;
  std::cout << "\n--- DDISinit for DatSet " << DataSetIndex << endl;
  // delete g_sock_DDIS;
  delete g_QCDNUM_pdf;
  g_sock_DDIS[DataSetIndex] = new sock_DDIS_t;
  g_QCDNUM_pdf = new qcdnum_pdf_t;
  string ccf("DDIS.coca");
  if(FileReadable(ccf)) g_sock_DDIS[DataSetIndex]->Params.Read(ccf.c_str());
  g_sock_DDIS[DataSetIndex]->Params.Show();
  g_sock_DDIS[DataSetIndex]->SetPhys(&g_PhysPars);
  g_sock_DDIS[DataSetIndex]->SetPDF(g_QCDNUM_pdf); 
  g_sock_DDIS[DataSetIndex]->Initialize();
  // g_sock_DDIS[DataSetIndex]->SetEcmSquared(4*27.5*920);
}

// =========================================
void setecmsq_(const int *DataSetIndex, const double *s) {
  g_sock_DDIS[*DataSetIndex]->SetEcmSquared(*s);
}

// =========================================
double ddisvalue_(const int *IDataSet, double *kinvars) {
  int DataSetIndex = *IDataSet;
  // --- needs ECM fixed or y specified for each data record
  //cout << "ddisvalue_: " << kinvars[0] << " "  << kinvars[1] << " "  << kinvars[2]
  //     <<"   "<< g_sock_DDIS[DataSetIndex]->Value(kinvars) << endl;
  if(!g_sock_DDIS.count(DataSetIndex)) {
    cout << "No DDIS model for DataSetIndex " << DataSetIndex << endl;
    exit(1);
  }
  return g_sock_DDIS[DataSetIndex]->Value(kinvars);
}
  
/**
  fitpars[3]:
  \arg Pomeron_a0 — Pomeron intercept
  \arg Reggeon_factor
  \arg Reggeon_a0 — Reggeon intercept
*/
// =========================================
void ddissetparams_(const double* fitpars) {
  // g_sock_DDIS[*DataSetIndex]->SetParams(fitpars);
  for (dict_t::iterator it=g_sock_DDIS.begin(); it!=g_sock_DDIS.end(); ++it) (it->second)->SetParams(fitpars);

}
