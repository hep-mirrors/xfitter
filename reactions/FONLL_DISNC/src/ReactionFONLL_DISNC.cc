/*
   @file ReactionFONLL_DISNC.cc
   @date 2017-11-29
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-11-29
*/

#include <iostream>
#include "ReactionFONLL_DISNC.h"
#include "xfitter_cpp.h"

// APFEL C++ interface header
#include "APFEL/APFEL.h"

// The class factory
extern "C" ReactionFONLL_DISNC* create()
{
  return new ReactionFONLL_DISNC();
}

extern "C" {
  void APFEL_set_pdfs( pXFXlike xfx);  //! Set PDFs
}


// Initialize at the start of the computation
int ReactionFONLL_DISNC::initAtStart(const string &s)
{
  // Retrieve parameters needed to initialize APFEL.
  const double MCharm     = GetParam("mch");
  const double MBottom    = GetParam("mbt");
  const double MTop       = GetParam("mtp");
  const double Mz         = GetParam("Mz");
  const double Mw         = GetParam("Mw");
  const int    PtOrder    = OrderMap(GetParamS("Order")) - 1;
  const double sin2thw    = GetParam("sin2thW");
  const double Vud        = GetParam("Vud");
  const double Vus        = GetParam("Vus");
  const double Vub        = GetParam("Vub");
  const double Vcd        = GetParam("Vcd");
  const double Vcs        = GetParam("Vcs");
  const double Vcb        = GetParam("Vcb");
  const double Vtd        = GetParam("Vtd");
  const double Vts        = GetParam("Vts");
  const double Vtb        = GetParam("Vtb");
  const double gf         = GetParam("gf");
  const double Q_ref      = Mz;
  const double Alphas_ref = GetParam("alphas");

  // FONLL-specific settings
  const string scheme     = "FONLL-" + GetParamS("FONLLVariant");
  const string MassScheme = GetParamS("MassScheme");
  const bool   runm       = ( GetParamI("Running") == 0 ? false : true);

  if (PtOrder == 0)
    {
      const string msg = "F: FONLL at LO not available. Use the ZM-VFNS instead.";
      hf_errlog_(17120601,msg.c_str(), msg.size());
    }
  else if (PtOrder == 1 && scheme == "FONLL-C")
    {
      const string msg = "F: At NLO only the FONLL-A and FONLL-B schemes are possible";
      hf_errlog_(17120602,msg.c_str(), msg.size());
    }
  else if (PtOrder == 3 && (scheme == "FONLL-C" ||  scheme == "FONLL-B"))
    {
      const string msg = "F: At NNLO only the FONLL-C scheme is possible";
      hf_errlog_(17120603,msg.c_str(), msg.size());
    }

  // If the MSbar masses are being used check that APFEL is used also
  // for the evolution.
  if (MassScheme == "MSbar")
    {
      const string msg = "F: When using MSbar heavy quark masses it is necessarey to use APFEL for the DGLAP evolution";
      hf_errlog_(17120604,msg.c_str(), msg.size());
    }

  // Set Parameters
  APFEL::SetZMass(Mz);
  APFEL::SetWMass(Mw);
  APFEL::SetSin2ThetaW(sin2thw);
  APFEL::SetGFermi(gf);
  APFEL::SetCKM(Vud, Vus, Vub,
		Vcd, Vcs, Vcb,
		Vtd, Vts, Vtb);
  APFEL::EnableDynamicalScaleVariations(true);
  if (MassScheme == "Pole")
    APFEL::SetPoleMasses(MCharm, MBottom, MTop);
  else if (MassScheme == "MSbar")
    {
      APFEL::SetMSbarMasses(MCharm, MBottom, MTop);
      APFEL::EnableMassRunning(runm);
    }
  APFEL::SetAlphaQCDRef(Alphas_ref, Q_ref);
  APFEL::SetPerturbativeOrder(PtOrder);
  APFEL::SetMassScheme(scheme);
  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1, 50, 3, 9.8e-7);
  APFEL::SetGridParameters(2, 40, 3, 1e-2);
  APFEL::SetGridParameters(3, 40, 3, 7e-1);

  // MSbar mass settings (Hard-coded for now).
  // They essentially set the scale at which MSbar masses are set.
  // The scale scale does not have to do the same of the mass itself.
  const double Qc = MCharm;
  const double Qb = MBottom;
  const double Qt = MTop;
  if (Qc != MCharm || Qb != MBottom || Qt !=MTop)
    {
      APFEL::InitializeAPFEL();
      double McQ = Qc;
      double MbQ = Qb;
      double MtQ = Qt;
      if (Qc != MCharm)
	McQ = APFEL::HeavyQuarkMass(4, Qc);
      if(Qb != MBottom)
	MbQ = APFEL::HeavyQuarkMass(5, Qb);
      if(Qt != MTop)
	MtQ = APFEL::HeavyQuarkMass(6, Qt);
      if (MassScheme == "Pole")
	APFEL::SetPoleMasses(McQ, MbQ, MtQ);
      else if (MassScheme == "MSbar")
	{
	  APFEL::SetMSbarMasses(McQ, MbQ, MtQ);
	  APFEL::SetMassScaleReference(Qc, Qb, Qt);
	}
    }

  // Initialize the APFEL DIS module
  APFEL::InitializeAPFEL_DIS();

  return 0;
}

// Compute all predictions in here and store them to be returned
// by the specific functions.
void ReactionFONLL_DISNC::initAtIteration()
{
  // VB: With the following command, APFEL will be calling the "ExternalSetAPFEL1"
  // routine in FONLL/src/FONLL_wrap.f. This is not optimal but until that routine is
  // there, I cannot find a way to override it.
  APFEL::SetPDFSet("external1");

  // Also make sure that proper PDFs are taken by external1 function which is located in FONLL/src directory
  APFEL_set_pdfs( getXFX() );

  APFEL::SetProcessDIS("NC");
  // Loop over the data sets.
  for ( auto dataSetID : _dsIDs)
    {
      // Charge of the projectile.
      // This does not have any impact because the difference between
      // "electron" and "positron" for a NC cross section is relevant
      // only when constructing the reduced cross section. But we keep
      // it here for clarity.
      const double charge = GetCharge(dataSetID);

      if (charge < 0)
      	APFEL::SetProjectileDIS("electron");
      else
       	APFEL::SetProjectileDIS("positron");

      // Get x,Q2 arrays.
      auto *q2p = GetBinValues(dataSetID,"Q2");
      auto *xp  = GetBinValues(dataSetID,"x");
      auto q2   = *q2p;
      auto x    = *xp;

      const size_t Np = GetNpoint(dataSetID);
      // Resize arrays.
      _f2fonll[dataSetID].resize(Np);
      _flfonll[dataSetID].resize(Np);
      _f3fonll[dataSetID].resize(Np);

      double Q2save = 0;
      for (size_t i = 0; i < Np; i++)
	{
	  // Skip all points with Q2 < 1 GeV^2.
	  if (q2[i] < 1)
	    continue;

	  // Recompute structure functions only if the value of Q2
	  // changes.
	  if (q2[i] != Q2save)
	    {
	      const double Q = sqrt(q2[i]);
	      APFEL::ComputeStructureFunctionsAPFEL(Q,Q);
	    }

	  // Compute structure functions by interpolation in x for the
	  // appropriate component (total, charm, or bottom).
	  switch (GetDataFlav(dataSetID))
	    {
	    case dataFlav::incl:
	      _f2fonll[dataSetID][i] = APFEL::F2total(x[i]);
	      _flfonll[dataSetID][i] = APFEL::FLtotal(x[i]);
	      _f3fonll[dataSetID][i] = - charge * APFEL::F3total(x[i]);
	      break;
	    case dataFlav::c:
	      _f2fonll[dataSetID][i] = APFEL::F2charm(x[i]);
	      _flfonll[dataSetID][i] = APFEL::FLcharm(x[i]);
	      _f3fonll[dataSetID][i] = - charge * APFEL::F3charm(x[i]);
	      break;
	    case dataFlav::b:
	      _f2fonll[dataSetID][i] = APFEL::F2bottom(x[i]);
	      _flfonll[dataSetID][i] = APFEL::FLbottom(x[i]);
	      _f3fonll[dataSetID][i] = - charge * APFEL::F3bottom(x[i]);
	      break;
	    }
	  Q2save = q2[i];
	}
    }
}

// FONLL structure functions
void ReactionFONLL_DISNC::F2 BASE_PARS 
{
  val = _f2fonll[dataSetID];
}

void ReactionFONLL_DISNC::FL BASE_PARS 
{
  val = _flfonll[dataSetID];
}

void ReactionFONLL_DISNC::xF3 BASE_PARS 
{
  val = _f3fonll[dataSetID];
}
