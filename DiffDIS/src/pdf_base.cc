/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

// #include <iomanip>
#include "pdf_base.h"
#include <fstream>
#include "DataTable.h"

/// @cond NIEWAZNE
// #define THE_CLASS pdf_base_t
// #define CAT(a,b) a##b
// #define QUAL_NAME(n) CAT(pdf_basic_t,NN##n)
// #define QUAL_NAME(n) pdf_basic##_t:##:n
#define QUAL_NAME(n) pdf_base_t:##:n
/// @endcond

// #define QCD_ORDER Phys.order

//--- If needed, define these in the makefile
//#define DE_BUG 1
//#define X1_TEST
//#define TEST_WGT
//#define TRACE_FLAVORS
//#define SHOW_INITIAL


  //===================================================================
  void QUAL_NAME(AllocWksp_1)(Distr_t fi[], int pn) {
    // --- new size is [Qgrid.npt +?1][X_GRD.npt]
    if(fi[pn] && Qgrid.GetN() == M_nQ && X_GRD.GetN() == M_nX) return; 
    M_nQ = Qgrid.GetN();
    M_nX = X_GRD.GetN();
    int nq = M_nQ + M_nQextra;
    DBG_SHOW(pn)
    try {
      // for(pn=0; pn < n_distr; pn++) {
        fi[pn] = new RealPtr_t[nq];
        //cout << fi[pn][0] << endl;
        for(int iq = 0; iq < nq; iq++) {
          fi[pn][iq] = new real_type[M_nX];
          memset(fi[pn][iq],0,M_nX*sizeof(real_type));
        }
        // M_used_distr_mem[pn] = 1;
      // }
    }
    catch(...) { emsg(1, "AllocWksp_1: could not allocate."); }
    // M_nAllocDistr++;
  }

  //===================================================================
  void QUAL_NAME(AllocWksp)(Distr_t* f) {
    // M_nAllocDistr = M_nFlavors + 2 + (M_FlavorSymm ? 0 : M_nFlavors); //--- NO space for q_sum
    if(Qgrid.GetN() != M_nQ || X_GRD.GetN() != M_nX) DeAllocWksp(f); 
    // M_nQ=Qgrid.GetN();
    for(int pn=M_pn_min; pn <= M_nFlavors; pn++) AllocWksp_1(f, pn);
  }

  //===================================================================
  void QUAL_NAME(DeAllocWksp_1)(Distr_t fi[], int pn) {
    DBG_SHOW(pn)
    // int iq;
    // for(pn=gluon; pn < n_distr; pn++) {
      for(int iq = 0; iq < M_nQ; iq++) delete[] fi[pn][iq];
      delete[] fi[pn];
      fi[pn] = NULL;
    // }
    //cout << "DeAlloc_q: 2" << endl;
  }

  //===================================================================
  void QUAL_NAME(DeAllocWksp)(Distr_t* f) {
    for(int pn=M_pn_min; pn <= M_nFlavors; pn++) if(f[pn]) DeAllocWksp_1(f, pn);
    if(f[q_sum]) DeAllocWksp_1(f, q_sum);
    // M_nAllocDistr = 0;
    M_nQ = 0;
    M_nX = 0;
    // memset(M_used_distr_mem, 0, sizeof(M_used_distr_mem));
  }

  //=================================
  void QUAL_NAME(ClearAll)() {
    int pn, iq;
    for(pn=M_pn_min; pn <= M_nFlavors; pn++) 
      if(isA(pn)) for(iq = 0; iq < M_nQ; iq++) memset(fi[pn][iq],0,X_GRD.GetN()*sizeof(real_type));
    if(isA(q_sum)) for(iq = 0; iq < M_nQ; iq++) memset(fi[q_sum][iq],0,X_GRD.GetN()*sizeof(real_type));
  }
  
/**
  \brief fill \c f with x^xpow_out * PDF(x,QQ).
  \internal
  \todo check M_PDFcurmode and thresholds
  \endinternal
*/
//===========================================================
void QUAL_NAME(Vals)(real_type xx, real_type QQ, real_type f[], int xpow_out) {
	Distr_t* distr=fi;
  real_type t, fac;
  int pn;

  t = Qgrid.tVal(QQ);
  int it = FindIndex(Qgrid.tval, M_nQ, t);
  if(it < 0 && fabs(Qgrid.tlo - t) < 1e-4) {it = 0; t = Qgrid.tlo;}
  require((it >= 0) && (it < M_nQ-1),
  //require(t >= Qgrid.tlo && t < Qgrid.thi,
    "pdf_basic_t::Val: Q^2 = %lf outside declared Q^2 range [%lf, %lf]",
    (double)QQ, (double)Qgrid[0], (double)Qgrid[M_nQ-1]
  );

  // xfac = pow(xx, xpow-Calc.x_exp);
  fac = pow(xx, xpow_out-M_Norm_xpow)/M_Norm_fac;
  int activeflav = Qgrid.GetNFlavors(it);
  //--- 2-dim interpolation, lin. in t
  for(pn = M_FlavorSymm ? gluon : -activeflav; pn <= activeflav; pn++)
    if(isA(pn)) f[pn] = fac*Ipol2(t, xx, Qgrid.tval, X_GRD.GetX(), distr[pn], M_nQ, M_nX, 1, 2);

  for(pn=activeflav+1; pn <= M_nFlavors; pn++) f[pn] = 0.0;
  if(!M_FlavorSymm) for(pn=activeflav+1; pn <= M_nFlavors; pn++) f[-pn] = 0.0;

  /*
  for(pn=gluon; pn <= activeflav; pn++) if(f[pn] < 0.0) {
  #ifdef NEG_PD_WARN
    printf("Negative PD(%d,%lg,%lg) = %lg\n",pn, xx, QQ, f[pn]);
  #endif
  #ifdef NONNEG_PD
    #ifdef NEG_PD_WARN
    printf("Set to 0.0\n");
    #endif
    f[pn] = 0.0;
  #endif
  }
  */
}

#ifndef doopa

/**
  \brief x^xpow_out * PDF(x,QQ)
  \internal
  \todo check M_PDFcurmode and thresholds
  \endinternal
*/
//===========================================================
real_type QUAL_NAME(Val)(real_type xx, real_type QQ, int pn, int xpow_out) {
  if(M_FlavorSymm && pn < 0) {
    if(M_PDFcurmode != mode_Quark) return 0;
    pn = -pn;
  }
	if(!isA(pn)) throw Fiasco("Val: PDF(%d) not available.", pn);
  real_type t, fac;

  t = Qgrid.tVal(QQ);
  int it = FindIndex(Qgrid.tval, Qgrid.npt, t);
  if(it < 0 && fabs(Qgrid.tlo - t) < 1e-4) {it = 0; t = Qgrid.tlo;}
  require((it >= 0) && (it < Qgrid.npt-1),
  //require(t >= Qgrid.tlo && t < Qgrid.thi,
    "Val: Q^2 = %lf outside declared Q^2 range [%lf, %lf]",
    (double)QQ, (double)Qgrid.qqVal(Qgrid.tlo), (double)Qgrid.qqVal(Qgrid.thi)
  );

  int activeflav = Qgrid.GetNFlavors(it);
  if(pn != q_sum && (pn > activeflav || pn < -activeflav)) return 0.0;
  fac = pow(xx, xpow_out-M_Norm_xpow)/M_Norm_fac;
  //--- 2-dim interpolation, lin. in t
  return fac*Ipol2(t, xx, Qgrid.tval, X_GRD.GetX(), fi[pn], M_nQ, M_nX, 1, 2);
  // return Ipol2(t, xx, Qgrid.tval, X_GRD.x, fi[pn], Qgrid.npt, X_GRD.nx, 1, 2)
			      // *xfac /(2*M_PI);
}
#endif

#ifndef NO_XDR
  // #define RW_ERR_MSG "pdf_base_t::RW error." 
  // #define XDR_R(a) require(xdr.RWdouble(&a),RW_ERR_MSG);
  // #define XDR_I(a) require(xdr.RWint(&a), RW_ERR_MSG);
  #define XDR_R(a) if(!xdr.RWdouble(&a)) return 0;
  #define XDR_I(a) if(!xdr.RWint(&a)) return 0;
    
  //==========================================================
  bool QUAL_NAME(RW1)(XDRio_t& xdr, int pn) {
    /*
      Format:
      pdf table [Qgrid.npt][X_GRD.npt]
    */
    int iq;
    // double x;
    //cout << "pn = "<< pn << endl;
    // XDR_I(iq)
    // require(iq==pn, "Qgrid read: bad parton index");
    if(xdr.isReading()) AllocWksp_1(fi,pn);
    for(iq = 0; iq < Qgrid.GetN(); iq++) {
      //cout << iq << endl;
      //cout << X_GRD.nx << endl;
      // require(xdr.RWdoublev(fi[pn][iq], X_GRD.GetN()), RW_ERR_MSG);
      if(!xdr.RWdoublev(fi[pn][iq], X_GRD.GetN())) return false;
      //cout << fi[pn][iq] << endl;
      //cout << fi[pn][iq][0] << endl;
    }
    // --- check stamp
    // XDR_R(x)
    // XDR_I(iq)
    // //require(iq==INTstamp && fabs(x-DBLstamp)<1e-10, "LoadGrids: chk failed.");
    // require(iq==INTstamp && x==DBLstamp, "LoadGrid: check failed.");
    return true;
  }
  
  /**
    \internal
    \todo save mode and normalization 
      PDFmode_t M_PDFmode, M_PDFcurmode;
      
  real_type M_Norm_fac;
  int M_Norm_xpow; 
  \endinternal
  */
  //==========================================================
  bool QUAL_NAME(RW)(XDRio_t& xdr) {
    /*
      Format:
      Qgrid
      Xgrid
      pdf tables:
        parton #
        fi[Qgrid.npt][X_GRD.npt]
        check-stamp
        ...
    */
    int pn;
    // require(Qgrid.RW(xdr), "Qgrid read failed.");
    // require(X_GRD.RW(xdr), "Xgrid read failed.");
    if(!Qgrid.RW(xdr)) return false;
    if(!X_GRD.RW(xdr)) return false;
    // if(xdr.isReading()) AllocWksp(fi);
    
    int iq = INTstamp;
    double x = DBLstamp;
    // string s1;
    // if(xdr.isReading()) s1 = "read"; else s1 = "written";
    // XDR_I(pn)
    if(xdr.isReading()) {
      // if(!Xgrid) Xgrid = new xgrid_base_t;
      // if(!X_GRD.RW(xdr)) return false;
      if(M_nQ && M_nQ != Qgrid.GetN()) DeAllocWksp(fi);
      if(M_nX && M_nX != X_GRD.GetN()) DeAllocWksp(fi);
      while(xdr.RWint(&pn)) {
        // if(!RW1(xdr, pn)) break;
        if(!RW1(xdr, pn)) return false;
        cout << "pdf(" << pn << ") read." << endl;
        // --- check stamp
        XDR_R(x)
        XDR_I(iq)
        // require(iq==INTstamp && x==DBLstamp, RW_ERR_MSG);
        if(iq != INTstamp || x != DBLstamp) return false;
      }
    }
    else {
      // int pn_min = M_FlavorSymm? 0 : -M_nFlavors;
      for(pn=M_pn_min; pn <= M_nFlavors; pn++) {
        if(!isA(pn)) continue;
        XDR_I(pn);
        if(!RW1(xdr, pn)) return false;
        cout << "pdf(" << pn << ") written." << endl;
        // --- check stamp
        XDR_R(x)
        XDR_I(iq)
      }
    }
    return true;
  }
  
  #undef XDR_I
  #undef XDR_R
#endif

/**
  \internal
  \todo save value of \c M_PDFcurmode
*/
//==========================================================
void QUAL_NAME(Save)(const char* dname, int xpow_out) {
  int npx = Xgrid.GetN(), Nq = Qgrid.GetN();
  int pn0 = M_pn_min;
  int iq, pn, ix;
  double fac;
  //#define PREC 6
  //#define OWID 9+PREC
  char ComStr[] = "# ";

  cout << "Saving to '" << dname <<"'..." << endl;
  ofstream vfile(dname);
  if(!vfile.good()) throw Fiasco("Cannot write to '%s'", dname);
  // cout << "vfile " << vfile.good() << endl;
  //cout << "vfile " << vfile.good() <<", " << fileno(ofil) << endl;
  //vfile << setprecision(PREC) << scientific;
  vfile << ComStr << "Made by pdf_t::Save" << endl;
  // vfile.TimeStamp();
  vfile << "# nx = " << npx <<endl;
  vfile << "# Q^2 grid (" << Nq << " points):" << endl;
  // for(iq=0; iq < Nq; iq++)
    // vfile << "# " << Qgrid.getQQ(iq) << "  "
          // << Qgrid.GetNFlavors(iq) << endl;
  Qgrid.Show(vfile);

  //--- x dep. on  Q^2 grid
  //==============================

  DataTable_t tbl;
  string Tit("x");
  if(!M_FlavorSymm) Tit.append(",tb,bb,cb,sb,ub,db");
  Tit.append(",g,d,u,s,c,b,t,qsum");
  tbl.Create("QQ,"+Tit);
  // vfile.ColTitles(Tit.c_str());
  int row = 0;
  for(iq=0; iq < Nq; iq++) {
    // vfile << "#$QQg = " << Qgrid.getQQ(iq) << endl;
    for(ix=0; ix < npx; ix++) {
      row = tbl.AddRow();
      tbl["QQ"][row] = Qgrid.getQQ(iq);
      tbl["x"][row] = Xgrid[ix];
      fac = pow(Xgrid[ix], xpow_out-M_Norm_xpow)/M_Norm_fac;
      int icol;
      for(pn=pn0, icol=2; pn <= M_nFlavors; pn++, icol++) if(isA(pn))
        // vfile.wr(fi[pn][iq][ix]/(2*M_PI*pow(Xgrid[ix], Calc.x_exp-1)));
        tbl[icol][row] = fac*fi[pn][iq][ix];
      if(isA(q_sum)) tbl["qsum"][row] = fac*fi[q_sum][iq][ix];
      // vfile << endl;
    }
    // vfile << endl;
    // vfile << endl;
  }
  tbl.SetPrecision(6);
  tbl.SplitOutput("QQ");
  tbl.Show(vfile, Tit);
  vfile.close(); 

#if 0
  //--- Q^2 dependence on x sub-grid
	//==============================

  string fn(dname);
  fn.append("t");
  cout << "Saving to '" << fn <<"'..." << endl;
  oFDstream vfile(fn.c_str(), 6);

  vfile.comment("Made by pdf_t::Save");
  vfile.TimeStamp();
  vfile << "# nq = " << Nq <<endl;
  //vfile << "# Delta t = " << Qgrid.tStep <<endl;
  int dix=15;
  vfile << "# x values:";
  for(ix=0; ix < npx; ix +=dix) vfile <<"# "<< Xgrid[ix] << endl;

  string Tit("QQ");
  if(!M_FlavorSymm) Tit.append(",tb,bb,cb,sb,ub,db");
  Tit.append(",g,d,u,s,c,b,t,qsum");
  vfile.ColTitles(Tit.c_str());
  for(ix=0; ix < npx; ix +=dix) {
    vfile << "#$x = " << Xgrid[ix] << endl;
    fac = pow(Xgrid[ix], xpow_out-M_Norm_xpow)/M_Norm_fac;
	  for(iq=0; iq < Nq; iq++) {
      vfile.wr(Qgrid.getQQ(iq));
      for(pn=pn0; pn <= M_nFlavors; pn++) if(isA(pn))
        // vfile.wr(fi[pn][iq][ix]/(2*M_PI*pow(Xgrid[ix], Calc.x_exp-1)));
        vfile.wr(fac*fi[pn][iq][ix]);
      if(isA(q_sum)) vfile.wr(fac*fi[q_sum][iq][ix]);
      vfile << endl;
	  }
    vfile << endl;
    vfile << endl;
  }
  //vfile.close();
#endif
  cout << "done." << endl;
}

#undef QUAL_NAME
