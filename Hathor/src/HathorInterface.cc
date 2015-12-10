#include <iostream>
#include <map>
#include <stdio.h>
#include <cstring>
#include "Hathor.h"
#include "../interface/xFitterPdf.h"

extern "C" {
  int hathorinit_(const int* idataset, const double& sqrtS, const bool& ppbar, const double& mt,
		  const unsigned int& pertubOrder, const unsigned int& precisionLevel);
  int hathorcalc_(const int *idataset, double *xsec);
}

extern "C" {
  int rlxd_size(void);
  void rlxd_get(int state[]);
  void rlxd_reset(int state[]);
  void rlxd_init(int level,int seed);
}

extern "C" {
  int hf_errlog_(const int* ID, const char* TEXT, long length);
}

// FIXME: delete pointers at the end! (in some hathordestroy_ or so)
std::map<int, Hathor*> hathor_array;   
xFitterPdf* pdf;
int* rndStore;
double mtop;

int hathorinit_(const int* idataset, const double& sqrtS, const bool& ppbar, const double& mt,
		const unsigned int& pertubOrder, const unsigned int& precisionLevel) {

  if(hathor_array.size()==0) {
    pdf = new xFitterPdf();

    mtop = mt;
    std::cout << " Top mass and renorm./fact. scale used for Hathor [GeV]: " << mtop << std::endl;

    rlxd_init(1,1);
    int nRnd = rlxd_size();
    //std::cout << " Size of random number array = " << nRnd << "\n";
    rndStore = new int [nRnd];
    rlxd_get(rndStore);
  }

  Hathor* hathor = new Hathor(*pdf);

  if(ppbar)
    hathor->setColliderType(Hathor::PPBAR);
  else
    hathor->setColliderType(Hathor::PP);

  hathor->setSqrtShad(sqrtS);

  unsigned int scheme = Hathor::LO;
  if(pertubOrder>1)
    scheme = scheme | Hathor::NLO;
  if(pertubOrder>2)
    scheme = scheme | Hathor::NNLO;
  hathor->setScheme(scheme);

  hathor->PrintOptions();

  hathor->setPrecision((int)pow(10,2+precisionLevel));

  hathor_array.insert(std::pair<int, Hathor*>(*idataset, hathor));

  return 0;
}

int hathorcalc_(const int *idataset, double *xsec) {
  rlxd_reset(rndStore);

  std::map<int, Hathor*>::const_iterator hathorIter = hathor_array.find(*idataset);

  if(hathorIter==hathor_array.end()) {
    const int id = 12040501;
    char text[256];
    sprintf(text, "S: Can not find HathorInterface for DataSet: %d", *idataset);
    hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
  }

  hathorIter->second->getXsection(mtop, mtop, mtop);

  double val,err;
  hathorIter->second->getResult(0,val,err);

  xsec[0] = val;
  
  return 0;
}
