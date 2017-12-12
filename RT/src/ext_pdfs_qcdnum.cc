/**
  @file ext_pdfs_qcdnum.cc

  @brief Direct link to QCDNUM PDFs and alphaS calls for RT codes

  Uses ext_pdfs.cc calls

  @version 0.1
  @date 2017/04/16
 */

#include "get_pdfs.h"
#include "ext_pdfs.h"

extern "C" {
  void set_qcdnum_pdfs_rt_();
}

void set_qcdnum_pdfs_rt_() {
  rt_set_pdfs_alphaS( &HF_GET_PDFSQ_WRAP, &HF_GET_ALPHASQ_WRAP );
}
