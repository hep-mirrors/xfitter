#include <iostream>
#include <iomanip>
#include <string>

#include "DYcalc.h"

#include "PhysPar.h"
#include "PDFconv.h"
#include "BinMatrix.h"


using namespace std;


DYcalc::DYcalc(BinMatrix *bm, PDFconv* pc, const IntSteps* ist) 
      : IntSteps(*ist)
{
  _bm = bm;
  _pc = pc;

  _nbins = _bm->getNbins();
  // even - Wm, odd - Wp
  if ( string("W") == _boz ) _nbins*=2; 

  _bin_int = new double[_nbins];
  for ( int ib = 0; ib<_nbins; ib++ ) _bin_int[ib] = 0.;

  intY=NULL;
  if ( string("W") == _boz ){
    intY=&DYcalc::intY_W;
  } else if ( string("Z") == _boz ){
    if ( string("y") == _var_name ) {
      intY=&DYcalc::intYbins_Z;
    } else {
      intY=&DYcalc::intY_Z;
    }
  }

}

DYcalc::~DYcalc()
{
  delete[] _bin_int;
}

int DYcalc::Integrate()
{
  for (int ib = 0; ib<_nbins; ib++){
    _bin_int[ib]=0.;
  }

  double qya[_nbins], qyb[_nbins], qym[_nbins];
  for (int ib=0;ib<_nbins;ib++){
    qya[ib]=0.; qyb[ib]=0.; qym[ib]=0.;
  }

  (*this.*intY)(0, qya);
  for (int ims=0; ims<_nms; ims++){
    double a = _msteps[ims], b = _msteps[ims+1];

    (*this.*intY)(2*(ims+1), qyb);
    (*this.*intY)(2*ims+1, qym);
    for (int ib = 0; ib<_nbins; ib++){
      _bin_int[ib]+=(b-a)/6.*(qya[ib]+4.*qym[ib]+qyb[ib]);
      qya[ib]=qyb[ib];
    }
      //cout << _bin_int[7] << endl;
   // exit(1);
  }

  return 1;

}

int DYcalc::intY_W(const int imp, double *qy){
  const double scale(PhysPar::mw);

  double qca[_nbins], qcb[_nbins], qcm[_nbins];
  for (int ib=0;ib<_nbins;ib++){
    qca[ib]=0.; qcb[ib]=0.; qcm[ib]=0.;
  }

  double xfxc[4] = {0.}; // pdf convolutions for dir=-1,+1 and Wm, Wp
  //  dir \ chrg  |   - |   +
  //  1           |   0 |   1
  //  -1          |   2 |   3
  (_pc->*(_pc->getPDFconv))(imp, 0, 1., xfxc[0], xfxc[1]);
  (_pc->*(_pc->getPDFconv))(imp, 0, -1., xfxc[2], xfxc[3]);
  for (int ib=0;ib<_nbins;ib++){
    for (int icd = 0; icd<4; icd ++){
    // even - Wm, odd - Wp
      qca[ib] += (ib%2==icd%2?1.:0.)*xfxc[icd]*_bm->BM[imp][0][ib/2][icd/2];
    }
    qy[ib] = 0.;
  }

  for (int iys=0;iys<_nys;iys++){
    double ya = _ysteps[iys];
    double yb = _ysteps[iys+1];
    double xfxcb[4]={0.}, xfxcm[4]={0.};
    (_pc->*(_pc->getPDFconv))(imp, 2*(iys+1)  ,  1., xfxcb[0], xfxcb[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*(iys+1)  , -1., xfxcb[2], xfxcb[3]);
    (_pc->*(_pc->getPDFconv))(imp, 2*iys+1,  1., xfxcm[0], xfxcm[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*iys+1, -1., xfxcm[2], xfxcm[3]);
    for(int ib=0;ib<_nbins;ib++){
      qcb[ib] = 0.; qcm[ib] = 0.;
      for (int icd = 0; icd<4; icd ++){
          qcb[ib] += (ib%2==icd%2?1.:0.)*xfxcb[icd]*_bm->BM[imp][2*(iys+1)][ib/2][icd/2];
          qcm[ib] += (ib%2==icd%2?1.:0.)*xfxcm[icd]*_bm->BM[imp][2*iys+1][ib/2][icd/2];
      }
      qy[ib] += (yb-ya)/6.*(qca[ib]+4*qcm[ib]+qcb[ib]);
      qca[ib] = qcb[ib];
    }

  }
    
      
  return 1;
}

int DYcalc::intY_Z(const int imp, double *qy){

  const double scale(PhysPar::mz);
  double qca[_nbins], qcb[_nbins], qcm[_nbins];
  for (int ib=0;ib<_nbins;ib++){
    qca[ib]=0.; qcb[ib]=0.; qcm[ib]=0.;
  }
  double xfxc[4] = {0.}; // pdf convolutions for dir=-1,+1 and fl=d,u
  //  dir \ flav  |   d |   u
  //  1           |   0 |   1
  //  -1          |   2 |   3
  (_pc->*(_pc->getPDFconv))(imp, 0, 1., xfxc[0], xfxc[1]);
  (_pc->*(_pc->getPDFconv))(imp, 0, -1., xfxc[2], xfxc[3]);
  for (int ib=0;ib<_nbins;ib++){
    for (int icdf = 0; icdf<12; icdf ++){
      qca[ib] += xfxc[icdf%4]*_bm->BM[imp][0][ib][icdf];
    }
    qy[ib] = 0.;
  }

  for (int iys=0;iys<_nys;iys++){
    double ya = _ysteps[iys];
    double yb = _ysteps[iys+1];
    double pcb[4]={0.}, pcm[4]={0.};
    (_pc->*(_pc->getPDFconv))(imp, 2*(iys+1)  ,  1.,  pcb[0], pcb[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*(iys+1)  , -1.,  pcb[2], pcb[3]);
    (_pc->*(_pc->getPDFconv))(imp, 2*iys+1,  1.,  pcm[0], pcm[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*iys+1, -1.,  pcm[2], pcm[3]);
    for(int ieb=0;ieb<_nbins;ieb++){
      qcb[ieb] = 0.; qcm[ieb] = 0.;
      for (int icdf = 0; icdf<12; icdf ++){
        qcb[ieb] += pcb[icdf%4]*_bm->BM[imp][2*(iys+1)][ieb][icdf];
        qcm[ieb] += pcm[icdf%4]*_bm->BM[imp][2*iys+1][ieb][icdf];
      }
      qy[ieb] += (yb-ya)/6.*(qca[ieb]+4*qcm[ieb]+qcb[ieb]);
      qca[ieb] = qcb[ieb];
    }
  }
      
  
  return 1;
}

int DYcalc::intYbins_Z(const int imp, double *qy){

  const double scale(PhysPar::mz);
  double qca[_nbins], qcb[_nbins], qcm[_nbins];
  for (int ib=0;ib<_nbins;ib++){
    qca[ib]=0.; qcb[ib]=0.; qcm[ib]=0.;
  }
  double xfxc[4] = {0.}; // pdf convolutions for dir=-1,+1 and fl=d,u
  //  dir \ flav  |   d |   u
  //  1           |   0 |   1
  //  -1          |   2 |   3
  for (int ib=0;ib<_nbins;ib++){
    (_pc->*(_pc->getPDFconv))(imp, 2*_ssb[ib], 1., xfxc[0], xfxc[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*_ssb[ib], -1., xfxc[2], xfxc[3]);
    for (int icdf = 0; icdf<12; icdf ++){
      qca[ib] += xfxc[icdf%4]*_bm->BM[imp][2*_ssb[ib]][0][icdf];
    }
    qy[ib] = 0.;
  }

  for (int iys=0;iys<_nys;iys++){
    int ib;
    for (ib = _nbins-1; ib>=0; ib--){
      if ( iys >= _ssb[ib]) break;
    }
    double ya = _ysteps[iys];
    double yb = _ysteps[iys+1];
    double pcb[4]={0.}, pcm[4]={0.};
    (_pc->*(_pc->getPDFconv))(imp, 2*(iys+1)  ,  1., pcb[0], pcb[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*(iys+1)  , -1., pcb[2], pcb[3]);
    (_pc->*(_pc->getPDFconv))(imp, 2*iys+1,  1., pcm[0], pcm[1]);
    (_pc->*(_pc->getPDFconv))(imp, 2*iys+1, -1., pcm[2], pcm[3]);
    qcb[ib] = 0.; qcm[ib] = 0.;
    for (int icdf = 0; icdf<12; icdf ++){
      qcb[ib] += pcb[icdf%4]*_bm->BM[imp][2*(iys+1)][0][icdf];
      qcm[ib] += pcm[icdf%4]*_bm->BM[imp][2*iys+1][0][icdf];
    }
    qy[ib] += (yb-ya)/6.*(qca[ib]+4*qcm[ib]+qcb[ib]);
    qca[ib] = qcb[ib];
  }
  
  return 1;
}

void DYcalc::getCalcRes(double *calc_res)
{
  for (int ib=0; ib<_nbins; ib++){
    calc_res[ib] = _bin_int[ib];
  }
}
