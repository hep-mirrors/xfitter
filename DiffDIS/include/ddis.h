/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _DDIS_H_
#define _DDIS_H_

#include <cmath>

#include "genut.h"
#include "dbgtools.h"
#include "physparams.h"
#include "aqcd.h"

#include "hdpdf.h"
#include "trs.h"
#include "grvpi.h"

#include "FFstrFns.h"

#include "JScif.h"

#ifdef HI_TWISTS
#include "mosat.h"
#endif

using namespace std;

// extern hdpdf_t* g_dPdf;
// extern double g_ReggeonID;
extern AlphaS *alphas; // for StrFns but not set!!!

/**
  \brief diffractive DIS

  Inclusive diffractive F2 and FL:
  \li FFNS
  \li ZM-VFNS
  \li GM-VFNS ??? Thorne-Roberts
  
  Reduced cross-section --- various contributions:
  \arg Pomeron, Reggeon and full
  \arg for Pomeron: F2a, FLa, F2c, FLc, F2b, FLb, jsF2c, jsFLc, jsF2b, jsFLb
  
  Reggeon modelled ~ pion; GRV LO/HO --> PionOrder.
  
  Simple P/R flux.
  dpdf3(xi,beta,Q^2) = flux_IP(xi) * dpdf_IP(beta,Q^2)
  + Reggeon_factor*flux_IR(xi) * dpdf_IR(beta,Q^2) 
  
  \internal    zdpdf needed for grid input
    
    \f$ \def\Pom{{I\!\!P}} \def\Reg{{I\!\!R}}
    f^{D(3)}(\xi,\beta,Q^2) = \Phi_\Pom(\xi)\, f_\Pom(\beta,Q^2)
    + N_\Reg \Phi_\Reg(\xi)\, f_\Reg(\beta,Q^2) 
    \f$\n
    but this is NOT used here.
    
    \f$ \def\Pom{{I\!\!P}} \def\Reg{{I\!\!R}}
    F_{2/\rm L}^{D(3)}(\xi,\beta,Q^2) = \Phi_\Pom(\xi)\, F_{2/\rm L}^\Pom(\beta,Q^2)
    + N_\Reg \Phi_\Reg(\xi)\, F_{2/\rm L}^\Reg(\beta,Q^2) 
    \f$
    
    \f$ \def\Pom{{I\!\!P}} \def\Reg{{I\!\!R}}
    F_{2/\rm L}^{\Pom/\Reg}(\beta,Q^2) \f$ are calculated from the respective PDFs and
    the fluxes can be taken from the grid file or calculated independently.
    
    Input:\n
    \arg f_Pom
    \arg P and R fluxes, 5 params. each
    \arg N_Reg = Reggeon_factor
    \arg PionOrder, MassCorrection
*/
// oooooooooooooooooooooooooooooooooooooooooooo
class dDIS_t {
	// AlphaS *alphas;
	FFStrFns_t* FFSF;
	// int nExpPoints;
	// int nVarPars, nVarPars1, nVarPars2;
	// int ReggeonID; //--- 0 = Pomeron, 1 = pion
	double Lambda4tr;
  int MaxFlavor;
	struct GRVHO_t {
		double mc,mb,mt,Lam4;
	} GRVHO; // = {1.5, 4.5, 100, 0.2};
	// bool GridsLoaded;
	// bool CorrelatedChisq;
	#ifdef HI_TWISTS
    MoSat_t mosat;
		// double lam2p1, AFL4;
		// int Twist4Model;
	#endif

protected:	
	hdpdf_t* dPdf;
	
public:
	PhysParams_t Phys;
	double Reggeon_factor;
	int PionOrder, MassCorrection;
  pdf_base_t* PomPdf;
	Flux_t Pflux,Rflux; ///< for \c UseFitFlux == false
  
	bool UseFitFlux;
	double SF[10];
	
	// double cfPstd, cfRstd; //--- xIP*flux [*Reggeon_factor] using fitted fluxes for -1 < t < tmax_kin
	// double cfP, cfR; //--- xIP*flux [*Reggeon_factor] using local fluxes
	double xsecP, xsecR; ///< x-sec w/o x*flux factor
	// double xsredPstd, xsredRstd, xsredP, xsredR; //--- xIP*x-secs
	
	double CoeffP[2], CoeffR[2]; //--- xIP*flux [*Reggeon_factor] 
	double xSigRedP[2], xSigRedR[2]; //--- xIP*x-secs
	/*---------------------------------------
	index 0/1 means using 
	0: local fluxes
	1: fitted fluxes for -1 < t < tmax_kin
	----------------------------------------*/
	#define cfPstd CoeffP[1]
	#define cfRstd CoeffR[1]
	#define cfP CoeffP[0]
	#define cfR CoeffR[0]
	#define xsredPstd xSigRedP[1]
	#define xsredRstd xSigRedR[1]
	#define xsredP xSigRedP[0]
	#define xsredR xSigRedR[0]
	
protected:	
	//===========================================
	void FillSF(double x, double QQ);

private:	
	// void Pomxf(double x, double QQ, double f[]);

	//=========================================
	// inline double Sqr(double x) {return x*x;}

	//===========================================
	double jsF2(int q_id, double x, double QQ) {
		JSsetSF(2, q_id, Phys.m2[q_id]);
		return JackSmith(x, QQ);
	}

	//===========================================
	double jsFL(int q_id, double x, double QQ) {
		JSsetSF(0, q_id, Phys.m2[q_id]);
		return JackSmith(x, QQ);
	}

	//===========================================
	double SigmaRedPom(double y, double x, double QQ) {
		FillSF(x, QQ);
		return SF[0] - y*y/(1+(1-y)*(1-y))*SF[1];
	}

	//===========================================
	double SigmaRedPi(double y, double x, double QQ);

  /*
	//=============================================================
	double SigmaRedFit(double y, double xP, double beta, double QQ) {
		return dPdf->Pflux.f(xP)*SigmaRedPom(y, beta, QQ)
				+ dPdf->GetReggeon_factor()*dPdf->Rflux.f(xP)*SigmaRedPi(y, beta, QQ);
	}

	//=============================================================
	// double FluxVal(bool isFit, bool isReggeon, double xP) {
		// if(isFit) return (isReggeon ? dPdf->Rflux.f(xP) : dPdf->Pflux.f(xP));
	double FluxVal(bool isReggeon, double xP) {
		if(UseFitFlux) return (isReggeon ? dPdf->Rflux.f(xP) : dPdf->Pflux.f(xP));
		return (isReggeon ? Rflux.f(xP) : Pflux.f(xP));
	}
  */

public:

	//=========================================
  void Configure(const string& gridName);
  void Configure(pdf_base_t& pdf, int pi_ord, int mass_corr=0);
  void config_F();
	dDIS_t(const string& gridName = "");

	// double cfPstd, cfRstd; //--- x*flux [*Reggeon_factor] using fitted fluxes for -1 < t < tmax_kin
	// double cfP, cfR; //--- x*flux [*Reggeon_factor] using local fluxes
	// double xsecP, xsecR; //--- x-sec w/o x*flux factor
	// double xsredPstd, xsredRstd, xsredP, xsredR; //--- xIP*x-secs
	
	//=============================================================
	double SigmaRed(double y, double xP, double beta, double QQ) {
		return Pflux.f(xP)*SigmaRedPom(y, beta, QQ)
				// + dPdf->GetReggeon_factor()*Rflux.f(xP)*SigmaRedPi(y, beta, QQ);
				+ Reggeon_factor*Rflux.f(xP)*SigmaRedPi(y, beta, QQ);
	}

	//=============================================================
	void FillSigmaCoeffs(double xP) {
		cfPstd = xP*dPdf->Pflux.f(xP);
		cfRstd = xP*Reggeon_factor*dPdf->Rflux.f(xP);
		cfP = xP*Pflux.f(xP);
		cfR = xP*Reggeon_factor*Rflux.f(xP);
	}
	
	//=============================================================
	void FillSigmaRed(double y, double xP, double beta, double QQ) {
		FillSigmaCoeffs(xP);
		xsecP = SigmaRedPom(y, beta, QQ);
		xsecR = SigmaRedPi(y, beta, QQ);
		xsredPstd = cfPstd*xsecP;
		xsredRstd = cfRstd*xsecR;
		xsredP = cfP*xsecP;
		xsredR = cfR*xsecR;
	}

	//=============================================================
	double xSigmaRed(double y, double xP, double beta, double QQ) {
		FillSigmaRed(y, xP, beta, QQ);
		return xSigRedP[UseFitFlux]+xSigRedR[UseFitFlux];
	}

	#ifdef HI_TWISTS
    //===========================================
    double xSigmaRed_HiTwist(double y, double xP, double beta, double QQ) {
      return mosat.HiTwist(y, QQ, xP, beta);
    }
	
    // =============================================
    void SetHTxi0G(double xi0=MoSat_t::x0_GBW) {
      mosat.Setx0(xi0); 
    }
	
    // =============================================
    void SetHTxi0(double xi0=MoSat_t::x0_GBW) {
      mosat.Setx0(xi0, xi0); 
    }

    // =============================================
    void SetHTas0(double as0) {
      mosat.SetAlphaQCD(as0); 
    }
	#endif

  friend class plug_DDIS_t;

};

#endif
