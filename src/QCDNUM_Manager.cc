#include"QCDNUM/QCDNUM.h"
#include"QCDNUM_Manager.h"
#include<iostream>
using namespace std;
namespace xfitter{
bool isQCDNUMinitialized = false;
void initQCDNUM(){
  if (isQCDNUMinitialized) return;
  QCDNUM::qcinit(-6," ");
  isQCDNUMinitialized = true;
}

bool isZMSTFinitialized = false;

void requireZMSTF(){
  if (isZMSTFinitialized) return;
  initQCDNUM();

  int nwords;
  QCDNUM::zmfillw(nwords);

  isZMSTFinitialized = true;
}
}
