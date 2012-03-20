/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "aqcd.h"

// ===============================================
void AlphaS::bas_init(int ord, double Qc2, double Qb2, double Qt2) {
  int nf;
  if(ord >= 0 && ord <= MAX_QCD_ORDER) QCDord = ord;
  else throw Fiasco("Illegal QCD order: %d", ord);
  for(nf=3; nf <= MAX_FLAVORS; nf++) {
    beta[0][nf] = 11.0 -2.0/3.0*nf;
    beta[1][nf] = 102 -38.0/3.*nf;
    brat[nf] = beta[1][nf]/beta[0][nf];
  }

  TranScale[3] = TranScale[4] = Qc2;
  TranScale[5] = Qb2;
  TranScale[6] = Qt2;
  // for(int n=4; n <= MAX_FLAVORS; n++) cout <<"bas_init "<< n <<": " << sqrt(TranScale[n]) << endl;
}

// ----------   RG (beta fcn) mode  -------------------------

// ===============================================
void AlphaS_RG::Initialize(double a0, double QQ0, int ord,
                        double Qc2, double Qb2, double Qt2) {
  int nf, nf0;
  // cout.precision(10);
  DBG_SHOW(a0)
  DBG_SHOW(sqrt(QQ0))
  bas_init(ord, Qc2, Qb2, Qt2);
  nf0 = Afl(QQ0);
  etaTran[nf0] = alphainv(TranScale[nf0], 4*M_PI/a0, QQ0, nf0);
  DBG_SHOW(4*M_PI/etaTran[nf0])
  for(nf=nf0+1; nf <= MAX_FLAVORS; nf++)
    etaTran[nf] = alphainv(TranScale[nf], etaTran[nf-1], TranScale[nf-1], nf-1);
  for(nf=nf0-1; nf > 3; nf--)
    etaTran[nf] = alphainv(TranScale[nf], etaTran[nf+1], TranScale[nf+1], nf);
  etaTran[3] = etaTran[4];
  DBG_SHOW(nf0)
  // cout << "as= " << setprecision(10) << 4*M_PI/alphainv(QQ0, etaTran[nf0], TranScale[nf0], nf0) << endl;
  // cout << "as val = " << setprecision(10) << val(QQ0) << endl;
}


//----------   Lambda mode  -------------------------

// ===============================================================
int AlphaS_Lam::FindZero(double xini, double acc, double *x0) {
  double x,y,xL=xini,xR=xini,yL,yR;
  const int nMax=64;
  while((yL=dAlphaNext(xL)) < 0) xL /= 2;
  while((yR=dAlphaNext(xR)) > 0) xR *= 2;
  int n=0;
  if(yL*yR > 0) return 2;
  if(yL == 0) {*x0 = xL; return 0;}
  if(yR == 0) {*x0 = xR; return 0;}
  do {
    //x = (yL*xR - xL*yR)/(yL-yR);
    //--- bisection
    x = (xR + xL)/2;
    if((y = dAlphaNext(x)) == 0) break;
    if(y*yL > 0) {yL=y; xL=x;} else {yR=y; xR=x;}
  } while(++n < nMax && fabs(xL-xR) > acc);
  //*x0 = (yL*xR - xL*yR)/(yL-yR);
  *x0 = x;
  return fabs(xL-xR) > acc;
}

//============================================
double AlphaS_Lam::alpha_s_o2pi(double t, int nf) {
  t += DlnLam2[nf];
  double a = 2.0/beta[0][nf]/t;
  if(QCDord) a *= (1 - brat[nf]/beta[0][nf]*log(t)/t);
  return a;
}

//=============================================
double AlphaS_Lam::Alpha(double mu2, int Nf, double Lam2) {
  Lambda2[Nf] = Lam2;
  DlnLam2[Nf] = -log(Lambda2[Nf]);
  return alpha_s_o2pi(log(mu2), Nf);
}

//=============================================
double AlphaS_Lam::dAlphaNext(double Lam2) {
  /*
  cout << Nflav << endl;
  cout << beta[0][Nflav] << endl;
  cout << beta[1][Nflav] << endl;
  cout << brat[Nflav] << endl;
  cout << 2*M_PI*Alf0 << endl;
  cout << Alpha(m1sq,Nflav,Lam2) << endl;
  cout << m1sq << endl;
  cout << Lam2 << endl;
  */
	return Alf0 - Alpha(m1sq,Nflav,Lam2);
}

//=============================================
void AlphaS_Lam::FixLambdas(double Lam4) {
  double acc=1e-8, Lam2;
  int rc;
  Lambda2[4] = Lam4*Lam4;
  for(int nf=4; nf < MAX_FLAVORS; nf++) {
    m1sq = TranScale[nf+1];
    //printf("m^2 = %lg\n", m1sq);
    Nflav = nf;
    Lam2 = Lambda2[Nflav];
    Alf0 = Alpha(m1sq,Nflav,Lam2);
    //printf("Lam2 %lg\n", Lam2);
    Nflav++;
    rc=FindZero(Lam2, acc, &Lambda2[Nflav]);
    if(rc) throw Fiasco("FixLambdas rc = %d\n", rc);
    printf("-> Lambda[%d] = %lg\n", Nflav, sqrt(Lambda2[Nflav]));
  }
  m1sq = TranScale[4];
  //printf("m^2 = %lg\n", m1sq);
  Nflav = 4;
  Lam2 = Lambda2[Nflav];
  Alf0 = Alpha(m1sq,Nflav,Lam2);
  Nflav = 3;
  if(rc=FindZero(Lam2, acc, &Lambda2[Nflav]))
    fprintf(stderr, "FixLambdas rc = %d\n", rc);
  printf("-> Lambda[%d] = %lg\n", Nflav, sqrt(Lambda2[Nflav]));
}

//=============================================
void AlphaS_Lam::FixLambdas(double alpha0, double QQ0, int nf0) {
  double acc=1e-8, Lam2;
  int rc,nf;
  Alf0 = alpha0/(2*M_PI);
  m1sq = QQ0;
  Nflav = nf0;
  //cout << "3*beta[0] = " << 3*beta[0][Nflav] << endl;
  rc=FindZero(0.2, acc, &Lambda2[Nflav]);
    //fprintf(stderr, "FixLambdasA1 rc = %d\n", rc);
  if(rc) throw Fiasco("FixLambdasA1 rc = %d\n", rc);
  printf("-> Lambda[%d] = %lg\n", Nflav, sqrt(Lambda2[Nflav]));
  for(nf=nf0; nf < MAX_FLAVORS; nf++) {
    m1sq = TranScale[nf+1];
    //printf("m^2 = %lg\n", m1sq);
    Nflav = nf;
    Lam2 = Lambda2[Nflav];
    Alf0 = Alpha(m1sq,Nflav,Lam2);
    //printf("Lam2 %lg\n", Lam2);
    Nflav++;
    rc=FindZero(Lam2, acc, &Lambda2[Nflav]);
    if(rc) throw Fiasco("FixLambdasA2 rc = %d\n", rc);
    printf("-> Lambda[%d] = %lg\n", Nflav, sqrt(Lambda2[Nflav]));
  }
  for(nf=nf0-1; nf >= 3; nf--) {
    m1sq = TranScale[nf+1];
    //printf("m^2 = %lg\n", m1sq);
    Nflav=nf;
    Lam2 = Lambda2[Nflav+1];
    Alf0 = Alpha(m1sq,Nflav+1,Lam2);
    rc=FindZero(Lam2, acc, &Lambda2[Nflav]);
    if(rc) throw Fiasco("FixLambdasA3 rc = %d\n", rc);
    printf("-> Lambda[%d] = %lg\n", Nflav, sqrt(Lambda2[Nflav]));
  }
}
