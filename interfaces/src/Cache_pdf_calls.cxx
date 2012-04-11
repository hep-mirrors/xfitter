//
// Created 2/04/2012. Cache PDF calls
//
#include <map>
#include <stdio.h>
#include <cstring>
#include "get_pdfs.h"

// Fortran interface
extern "C" {
  void getcachedpdfs_(int *iteration, double *x, double *Q2, double *Pdfs);
  int  getcachesize_();
}

//
// x,Q2 grid point
// 
struct xQ2grid {
  double x;
  double Q2;
};

//
// PDFs for an iteration
// 
struct PDFsContainer {
  double PDFs[13];
  int    iteration;
};


//
// Compare operation
//
struct compare_xQ2 {
  bool operator()( const xQ2grid xq1, const xQ2grid xq2 ) const
  {
    return memcmp(&xq1,&xq2,sizeof(xQ2grid) ) < 0; // Do fast bit-wise comparison for sorting.
  }
};


static std::map<xQ2grid,PDFsContainer,compare_xQ2> pdf_cache;

// Cache PDF calls using associative map
void getcachedpdfs_(int *iteration, double *x, double *Q2, double *Pdfs) {
  //

  xQ2grid xQ2point;
  xQ2point.x  = *x;
  xQ2point.Q2 = *Q2;

  std::map<xQ2grid,PDFsContainer>::iterator it; 

  // Check if already present:
  it = pdf_cache.find(xQ2point);
  if ( it == pdf_cache.end() ) {
    // Not found. Add new point
    PDFsContainer pdf_store;
    // Get PDFs
    HF_GET_PDFS_UNCACHED_WRAP(x,Q2,Pdfs);

    // Store them
    for (int i=0; i<13; i++) {
      pdf_store.PDFs[i] = Pdfs[i];
    }
    // Add to the cache 
    pdf_store.iteration = *iteration;
    pdf_cache[xQ2point] = pdf_store;

    //    printf("Map size= %10i\n",pdf_cache.size());
  }
  else {
    // Got it. Check if up-to-date.
    if ( it->second.iteration == *iteration) {
      // Up-to-date, return:
      for (int i=0; i<13; i++) {
	Pdfs[i] = it->second.PDFs[i];
      }
    }
    else {
      // Out-of-date, update:
      HF_GET_PDFS_UNCACHED_WRAP(x,Q2,Pdfs);

      // Store them
      for (int i=0; i<13; i++) {
	it->second.PDFs[i] = Pdfs[i];
      }
      it->second.iteration = *iteration;
    }
  }
};

int getcachesize_(){
  return pdf_cache.size();
};


