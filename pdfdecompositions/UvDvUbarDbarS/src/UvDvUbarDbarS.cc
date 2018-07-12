 
/*
   @file UvDvUbarDbarS_PdfDecomposition.cc
   @date 2018-07-11
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "UvDvUbarDbarS.h"
#include "HERAPDF_PdfParam.h"

namespace xfitter
{
  //_________________________________________________________________________________
  UvDvUbarDbarS::UvDvUbarDbarS(): BasePdfDecomposition{"UvDvUbarDbarS"} { }

  //_________________________________________________________________________________
  void UvDvUbarDbarS::initAtStart(const std::string & pars) const
  {
  }

  //_________________________________________________________________________________
  std::function<std::map<int,double>(const double& x)> UvDvUbarDbarS::f0() const
  {
    return nullptr;
  }
}


