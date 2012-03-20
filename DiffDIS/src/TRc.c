#include <string.h>
//      COMMON/TRSFPAR/alambda,flavor,qsct,qsdt,iord
typedef struct {
  double Lambda, Nflavors, QQb4, QQc4;
  int ord;
} COM_TRSFPAR;
extern COM_TRSFPAR trsfpar_;

//      COMMON/GCONTR/F2pg,FLpg,F2cg,FLcg,F2bg,FLbg
typedef struct {
  double F2pg,FLpg,F2cg,FLcg,F2bg,FLbg;
} COM_GCONTR;
extern COM_GCONTR gcontr_;

void wate96z_();
void TRinit(double mc, double mb, int Nf, double Lam) {
  wate96z_();
  trsfpar_.ord = 1;
  trsfpar_.Nflavors = Nf;
  trsfpar_.QQc4 = 4*mc*mc;
  trsfpar_.QQb4 = 4*mb*mb;
  trsfpar_.Lambda = Lam;
}

//==============================================
//      subroutine sfun(x,q2,f2p,flp,f2c,flc,f2b,flb)
void sfunws_(double* x, double* q2,
  double* f2p, double* flp,
  double* f2c, double* flc,
  double* f2b, double* flb);

//==============================================
void TRsfun(double x, double QQ, double SF[]) {
  //printf("TR x = %g, QQ = %g\n", x, QQ);
  memset(SF, 0, 6*sizeof(SF[0]));
  if(x > 0.99999) return;
  sfunws_(&x, &QQ,
  //f2p,flp
  SF, SF+1,
  //f2c,flc,
  SF+2, SF+3,
  //f2b,flb
  SF+4, SF+5
  );

  /*

      subroutine sfun_wrap(x,q2,f2p,flp,f1p,rp,f2n,fln,f1n,rn,
     $     f2c,flc,f1c,f2b,flb,f1b,iflag,index,f2QCDNUM,flQCDNUM
     $     ,usekfactors)

     set usekfactors to false
     set f2QCDNUM,flQCDNUM to 0
     set index to 0
     set iflag to 0  
  */

}

//==============================================
void TRgetGcontr(double SFg[]) {
  SFg[0] = gcontr_.F2pg;
  SFg[1] = gcontr_.FLpg;
  SFg[2] = gcontr_.F2cg;
  SFg[3] = gcontr_.FLcg;
  SFg[4] = gcontr_.F2bg;
  SFg[5] = gcontr_.FLbg;
}
