/*
  F2_heavy via
  gamma* parton -> Q + Qbar + X folded with the parton densities
*/

/*
  params to fill:
      read(11,*) ipdf  --- not used
      read(11,*) ifch
      read(11,*) nborn
      read(11,*) ngcorr
      read(11,*) nqcorr
      read(11,*) qfac
      read(11,*) mfac
c      read(11,*) mq
c      read(11,*) mq2fix
      read(11,*) nlf
      read(11,*) ichoose
c      read(11,*) xbjfix
      read(11,87) q
 87   format(a10)

c      print*,' choice of distribution function: ipdf = ',ipdf
c      print*,' choice of f_l (0) f_t (1) f_2 (2): ifch = ',ifch
c      print*,' include the born terms (1) not (0): nborn = ',nborn
c      print*,' include the gluon corrections (1) not (0): ngcorr = '
c     #     ,ngcorr
c      print*,' include the quark corrections (1) not (0): nqcorr = '
c     #     ,nqcorr
c      print*,' factor multiplying Q^2 to give MF scale = ',qfac
c      print*,' factor multiplying m^2 to give MF scale = ',mfac
c      print*,' the mass of the produced heavy quark: mq = ',mq
c      print*,' the virtuality of the photon: mq2 = ', mq2fix
c      print*,' the number of light flavors: nlf = ', nlf
c      print*,' ichoose = (0) loop (1) one pt (2) mf dep = ',ichoose
c      print*,' the bjorken x fixed value: xbjfix = ', xbjfix
c      print*,' the produced quark is ',q

      common/include/nborn, ngcorr, nqcorr, ifch
      common/mass/m2, mq2, xbj
      common/coupling/scale, alphas, nlf, eh2
      COMMON/JSPARAM/qfac,mfac

*/

typedef struct {
  int nborn, ngcorr, nqcorr, ifch;
} COM_include;
extern COM_include include_;

typedef struct {
  double m2,QQ,x;
} COM_mass;
extern COM_mass mass_;

typedef struct {
  double scale, alphas, nlf, eh2;
} COM_coupling;
extern COM_coupling coupling_;

typedef struct {
  double qfac, mfac;
} COM_JSPARAM;
extern COM_JSPARAM jsparam_;

//===================================================
void JSsetFact(double qfac, double mfac) {
  jsparam_.qfac = qfac;
  jsparam_.mfac = mfac;
}

//===================================================
void JSsetContr(int born, int gluon, int quark) {
  //--- each born, gluon, quark = 0 or 1
  include_.nborn = born;
  include_.ngcorr = gluon;
  include_.nqcorr = quark;
}

//===================================================
void JSsetSF(int Ftype, int q_id, double mq2) {
  fixquark_(&q_id);
  mass_.m2 = mq2;
  include_.ifch = Ftype;
  coupling_.nlf = 3;
}

//===================================================
void JSsetParams(int Ftype, int q_id, double mq, double qfac, double mfac) {
  //--- qch = quark charge
  //ichoose = 1; fixed in jacksmisth.f
  JSsetSF(Ftype, q_id,mq*mq);
  JSsetFact(qfac, mfac);
  JSsetContr(1,1,1);
}

//===================================================
double JackSmith(double x, double QQ) {
  double fhad;
  if(QQ*(1/x-1) <= 4.00001*mass_.m2) return 0;
  //printf("JS x= %g, QQ = %g\n", x, QQ);
  js_(&x, &QQ, &fhad);
  return fhad;
}

