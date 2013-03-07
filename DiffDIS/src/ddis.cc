/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2013
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 2.11.01
  
  Changes:
  2012-02-24
    Pflux, Rflux set to fitted ones upon reading a grid.
_____________________________________________________________*/

#include "ddis.h"

static pdf_base_t* g_PomPdf;
static double g_ReggeonID;
static double g_PionOrder;

// --- global link for StrFns

  AlphaS *alphas; // for StrFns but not set!!!

  //=========================================
  void pi0xf(double x, double QQ, double f[]) {
  //--- f[-6:6]
    int k;
    //--- pion
    double xF[7]; //--- g, uval=dval, usea=dsea, s,c,b
    if(g_PionOrder) GRVpi_HO(x, QQ, xF); else  GRVpi_LO(x, QQ, xF);
    f[0] = xF[0];
    f[1] = f[2] = xF[1]/2 + xF[2];
    for(k=3; k <= 5; k++) f[k] = xF[k];
    f[6] = 0;
    for(k=-6; k < 0; k++) f[k] = f[-k];
  }

  //=========================================
  void Pomxf(double x, double QQ, double f[]) {
  //--- f[-6:6]
    int k;
    //--- Pomeron
    g_PomPdf->Vals(x, QQ, f, 1);
    for(k=-6; k < 0; k++) f[-k] = f[k];
    // f[6] = f[-6] =0;
  }


//================================================
void dDIS_t::Configure(const string& gridName) {
  dPdf = new hdpdf_t(gridName);
  PomPdf = dPdf;
  Pflux = dPdf->Pflux;
  Rflux = dPdf->Rflux;

  PionOrder = dPdf->GetPionOrder();
  Reggeon_factor = dPdf->GetReggeon_factor();
  MaxFlavor = dPdf->GetMaxFlavor();
  
  int order = 1;
  #define GET_MASS(jq) dPdf->Qgrid.GetFixedInd(jq) < 0 ? -1. : sqrt(dPdf->Qgrid[dPdf->Qgrid.GetFixedInd(jq)])
  Phys.set(0, order); //--- set masses to default values
  Phys.set(0, //--- VFNS == 0
    order,
    GET_MASS(c_quark),
    GET_MASS(b_quark),
    GET_MASS(t_quark)
  );
  #undef GET_MASS
  
  config_F();
}

//================================================
void dDIS_t::Configure(pdf_base_t& pdf, int pi_ord, int mass_corr) {
  PomPdf = &pdf;
  PionOrder = pi_ord;
  MassCorrection = mass_corr;
  MaxFlavor = PomPdf->GetMaxFlavor();
  config_F();
}

//================================================
void dDIS_t::config_F() {
  if(MaxFlavor < 3) throw "MaxFlavor must be at least 3";
  if(!Phys.nFixedFlavors && MaxFlavor < 5) throw "MaxFlavor must be at least 5";
  g_PomPdf = PomPdf;
  g_PionOrder = PionOrder;
  // Phys.isPhoton = 0;

  FFSF = new FFStrFns_t(Pomxf, Phys.QCDorder);

  //-------- Lambda4 for Thorne-Roberts
  AlphaS_Lam as;
  as.Initialize(ALPHAS_Z0, MASS_Z0*MASS_Z0, Phys);
  Lambda4tr = as.GetLambda(4);
}

//=========================================
dDIS_t::dDIS_t(const string& gridName) {
  MassCorrection = 0;
  MaxFlavor = 5;

  GRVHO.mc = 1.5;
  GRVHO.mb = 4.5;
  GRVHO.mt = 100;
  GRVHO.Lam4 = 0.2;
  
  UseFitFlux = 1;

  // --- ZEUS
  Pflux.SetParams(1.11, 0.0, 7.0, -0.55, -0.09);
  Rflux.SetParams(0.70, 0.9, 2.0, -0.55, -0.09);
  // cout << "Exp. Pomeron flux: " << Pflux;
  // cout << "Exp. Reggeon flux: " << Rflux;

  if(!gridName.empty()) Configure(gridName);
  
  //JSsetFact(1, 4);
  JSsetFact(1, 0);
  JSsetContr(1, 1, 1); //--- born, gluon, quark, i.e. LO + NLO

  //--- TEST -- passed
}

//===========================================
void dDIS_t::FillSF(double x, double QQ) {
  //--- SF = F2a, FLa, F2c, FLc, F2b, FLb, jsF2c, jsFLc, jsF2b, jsFLb
  //---       0    1    2    3    4    5     6      7      8      9
  g_ReggeonID = 0;
  TRinit(Phys.m[c_quark], Phys.m[b_quark], 3, Lambda4tr);
  memset(SF, 0, sizeof(SF));
  if(Phys.nFixedFlavors) {
    FFSF->SetPDF(Pomxf);
    SF[0] = FFSF->F2(x, QQ, Phys.nFixedFlavors);
    if(MaxFlavor >= c_quark) SF[0] += SF[6] = jsF2(c_quark, x, QQ);
    if(MaxFlavor >= b_quark) SF[0] += SF[8] = jsF2(b_quark, x, QQ);
    //--- FL:
    SF[1] = FFSF->FL(x, QQ, Phys.nFixedFlavors);
    if(MaxFlavor >= c_quark) SF[1] += SF[7] = jsFL(c_quark, x, QQ);
    if(MaxFlavor >= b_quark) SF[1] += SF[9] = jsFL(b_quark, x, QQ);
  }
  else {
    DBG_SHOW(MassCorrection)
    TRsfun(x, QQ, SF);
    if(MassCorrection) SF[0] +=
        (SF[6] = jsF2(c_quark, x, min(QQ, Phys.m2[c_quark])))
      + (SF[8] = jsF2(b_quark, x, min(QQ, Phys.m2[b_quark])));
    if(MassCorrection > 1) SF[1] +=
        (SF[7] = jsFL(c_quark, x, min(QQ, Phys.m2[c_quark])))
      + (SF[9] = jsFL(b_quark, x, min(QQ, Phys.m2[b_quark])));
  }
}

//===========================================
double dDIS_t::SigmaRedPi(double y, double x, double QQ) {
  g_ReggeonID = 1;
  if(PionOrder) {
    if(!Phys.nFixedFlavors) {
      double piSF[10];
      TRinit(GRVHO.mc, GRVHO.mb, 3, GRVHO.Lam4);
      TRsfun(x, QQ, piSF);
      if(MassCorrection) piSF[0] +=
          (piSF[6] = jsF2(c_quark, x, min(QQ, Phys.m2[c_quark])))
        + (piSF[8] = jsF2(b_quark, x, min(QQ, Phys.m2[b_quark])));
      if(MassCorrection > 1) piSF[1] +=
          (piSF[7] = jsFL(c_quark, x, min(QQ, Phys.m2[c_quark])))
        + (piSF[9] = jsFL(b_quark, x, min(QQ, Phys.m2[b_quark])));
      return piSF[0] - y*y/(1+(1-y)*(1-y))*piSF[1];
    } else {
      FFSF->SetPDF(pi0xf);
      double F2pi = FFSF->F2(x, QQ, Phys.nFixedFlavors);
      if(MaxFlavor >= c_quark) F2pi += jsF2(c_quark, x, QQ);
      if(MaxFlavor >= b_quark) F2pi += jsF2(b_quark, x, QQ);
      double FLpi = FFSF->FL(x, QQ, Phys.nFixedFlavors);
      if(MaxFlavor >= c_quark) FLpi += jsFL(c_quark, x, QQ);
      if(MaxFlavor >= b_quark) FLpi += jsFL(b_quark, x, QQ);
      return F2pi - y*y/(1+(1-y)*(1-y))*FLpi;
    }
  } else {//--- LO Reggeon contr.
    double xF[6]; //--- g, uval=dval, usea=dsea, s,c,b
    //const double
    GRVpi_LO(x, QQ, xF);
    double S = xF[1]/2 + xF[2];
    double F2pi9 = 5*S + xF[3];
    if(MaxFlavor >= c_quark) F2pi9 += 4*xF[c_quark];
    if(MaxFlavor >= b_quark) F2pi9 += xF[b_quark];
      //(5*xF[1] + 2*(5*xF[2] + xF[3] + 4*xF[4] + xF[5]))/9;
    //JSsetContr(1, 0, 0); //--- born, i.e. LO
    return F2pi9*2./9;
  }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FORTRAN link for
  trs/jacksmith.f
  trs/rt.f
Uses global: hdpdf_t* g_dPdf
  g_dPdf is defined by the dDIS_t class
  and is automatically set to the most recent dDIS_t object
*/

extern "C" {

  // //=========================================
  // char* f2cstr(char* name, int len) {
    // char *cp = name + len-1;
    // while(len && *cp == ' ') {cp--; len--;}
    // char* nn = (char*)malloc(len+1);
    // strncpy(nn, name, len);
    // nn[len] =0;
    // return nn;
  // }

  //=========================================
  void xpdfvs_(double* px, double* pQ2,
        double* upv, double* dnv, double* usea, double* dsea,
        double* str, double* chm, double* bot, double* glu) {
    double x = *px;
    double QQ = *pQ2;
    double xF[7]; //--- g, uval=dval, usea=dsea, s,c,b
    *upv = *dnv = 0;
    if(g_ReggeonID) {//--- pion
      GRVpi_HO(x, QQ, xF);
      *usea = *dsea = xF[1]/2 + xF[2];
      *str = xF[3];
      *chm = xF[4];
      *bot = xF[5];
      *glu = xF[0];
    }
    else {//--- Pomeron
      g_PomPdf->Vals(x, QQ, xF, 1);
      *str = *usea = *dsea = xF[1];
      *chm = xF[4];
      *bot = xF[5];
      *glu = xF[0];
    }
  }

}

#undef _DEBUG_
