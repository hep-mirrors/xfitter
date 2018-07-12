 
/*
   @file LHAPDF.cc
   @date 2018-07-12
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-12
*/

#include "LHAPDFDecomposition.h"
#include "xfitter_pars.h"

namespace xfitter
{
  //_________________________________________________________________________________
  LHAPDFDecomposition::LHAPDFDecomposition(const std::string& PDFset, const int& mem):
    BasePdfDecomposition{"LHAPDF"},
    _mem(mem)
    {
      // Upload PDF set from LHAPDF
      _dist = LHAPDF::mkPDFs(PDFset);
    }

  //_________________________________________________________________________________
  void LHAPDFDecomposition::initAtStart(const std::string & pars) const
  {
  }

  //_________________________________________________________________________________
  std::function<std::map<int,double>(const double& x)> LHAPDFDecomposition::f0() const
  {
    return [=] (const double& x)->std::map<int,double>
      {
	return _dist[_mem]->xfxQ(x, *(XFITTER_PARS::gParameters.at("Q0")));
      };
  }
}



