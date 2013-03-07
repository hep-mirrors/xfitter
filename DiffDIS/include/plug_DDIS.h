/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 1.1.3
_____________________________________________________________*/

#ifndef PLUG_DDIS_H_
#define PLUG_DDIS_H_

#include "ddis.h"
#include "TheModel_base.h"

/**
  Parameters of the model (besides PDFs):
  \arg Flux_tmin, Flux_tmax — t limits for the integrated flux
  \arg Pomeron_tslope — Pomeron flux t-slope
  \arg Pomeron_a0 — Pomeron intercept
  \arg Pomeron_a1 — Pomeron slope
  \arg Reggeon_tslope — Reggeon flux t-slope
  \arg Reggeon_a0 — Reggeon intercept
  \arg Reggeon_a1 — Reggeon slope
  \arg Reggeon_factor
  \arg Pion_order — QCD order of the pion PDFs
  \arg Mass_corr — type of fixed order mass corrections
  \arg PDFgrid — name of grids to read
*/
  // \arg PhysParams_t Phys

// oooooooooooooooooooooooooooooooooooooooo
class plug_DDIS_t : public TheModel_base_t {
protected:
  dDIS_t M_theory;
  bool hiTwist;
  
public:
  // ====================================================
  plug_DDIS_t() : TheModel_base_t("DDIS-P+R") {
    Params.SetName("DDIS");
    Params.SetQuoteChars();
    Params.SetVerbose(2);
    
    Params.Create("Flux_tmin", -1.0);
    Params.Create("Flux_tmax", 0.0);  //--- 0 means kin. max = -{\xi^2 \over 1-\xi} m_p^2
    Params.Create("Pomeron_tslope", 7.0);
    Params.Create("Pomeron_a0", 1.11);
    Params.Create("Pomeron_a1", 0.0);
    Params.Create("Reggeon_tslope", 2.0);
    Params.Create("Reggeon_a0", 0.70);
    Params.Create("Reggeon_a1", 0.9);
    Params.Create("Reggeon_factor", 2.70);
    Params.Create("Pion_order", 0);
    Params.Create("Mass_corr", 0);
    Params.Create("PDFgrid", "");
    Params.Create("HigherTwists", 0);
    
    // hiTwist = 0;
    hiTwist = Params.GetInt("HigherTwists");
  }
  
  // ====================================================
  void Initialize() {
    // cout << "==> plug_DDIS_t::Initialize" << endl;
    if(Params.GetString("PDFgrid").empty()) {
      // --- prepare for fitting new PDFs
      if(!Phys) throw Fiasco("plug_DDIS_t::Configure ERROR: Phys not set.");
      if(!Pdf) throw Fiasco("plug_DDIS_t::Configure ERROR: Pdf not set.");
      M_theory.Phys = *Phys;
    }
    else {
      // --- use PDFs grids
      M_theory.Configure(Params.GetString("PDFgrid"));
      Phys = &M_theory.Phys;
      Pdf = M_theory.dPdf;
    }
    isInitialized = true;
    Configure();
  }
  
  // ====================================================
  void Configure() {
    if(!isInitialized) Initialize();
    else {
      if(Params.GetString("PDFgrid").empty()) {
        M_theory.Configure(*Pdf, Params.GetInt("Pion_order"), Params.GetInt("Mass_corr"));
        M_theory.Pflux.SetParams(Params.GetReal("Pomeron_a0"), Params.GetReal("Pomeron_a1"), Params.GetReal("Pomeron_tslope")); 
        M_theory.Rflux.SetParams(Params.GetReal("Reggeon_a0"), Params.GetReal("Reggeon_a1"), Params.GetReal("Reggeon_tslope")); 
        M_theory.Reggeon_factor = Params.GetReal("Reggeon_factor");
      }
      hiTwist = Params.GetInt("HigherTwists");
      #ifndef HI_TWISTS
        if(hiTwist) throw "Higher twists constributions not available";
      #endif
      cout << "Pomeron flux: " << M_theory.Pflux;
      cout << "Reggeon flux: " << M_theory.Rflux;
    }
  }
  
  // // ====================================================
  // void SetHiTwists(bool ht=true) {
  // #ifdef HI_TWISTS
    // hiTwist = ht;
  // #else
    // throw "Higher twists construbutions not available";
  // #endif
  // }
  
};

#endif
