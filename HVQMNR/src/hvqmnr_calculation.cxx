#include <fortran_wrapper.h>
#include <hvqmnr_grid.h>
#include <hvqmnr_mnr.h>
#include <hvqmnr_frag.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2D.h>


namespace HVQMNR
{
  // Structure to store steering parameters
  struct Steering
  {
    double ptmin;
    double ptmax;
    int    npt;
    int    nptsm;
    double ymin;
    double ymax;
    int    ny;
    int    nsfnb;
    int    nx3;
    int    nx4;
    int    nbz;
  };

  // Initialise calculation with default parameters
  void DefaultInit(Steering& steer, double mq, MNR& mnr, Frag& frag, Grid& grid, Grid& grid_smoothed)
  {
    // MNR (parton level cross sections)
    mnr.bFS_Q = true;
    mnr.bFS_A = true;
    // x3 and x4 binning
    mnr.fBn_x3 = steer.nx3;
    mnr.fBn_x4 = steer.nx4;
    mnr.fSF_nb = steer.nsfnb;
    mnr.CalcBinning();
    // Number of flavours
    mnr.fC_nl = 3;
    // Parton level pT-y grids
    grid.SetL(steer.npt, steer.ptmin, steer.ptmax, mq);
    grid.SetY(steer.ny, steer.ymin, steer.ymax);
    grid.SetW(1);
    grid_smoothed.SetL(steer.nptsm, steer.ptmin, steer.ptmax, mq);
    grid_smoothed.SetY(steer.ny, steer.ymin, steer.ymax);
    grid_smoothed.SetW(1);
    // Fragmentation
    frag.SetNz(steer.nbz);
  }
  
  // Fill output arrays
  int FillOutput(TH1* hcs, double* xsec, double* ymin, double* ymax, double* ptmin, double* ptmax)
  {
    int bin = 0;
    for(int bpt = 0; bpt < hcs->GetNbinsX(); bpt++)
    {
      for(int by = 0; by < hcs->GetNbinsY(); by++) 
      {
        xsec[bin]  = hcs->GetBinContent(bpt+1, by+1);
        ymin[bin]  = hcs->GetYaxis()->GetBinLowEdge(by+1);
        ymax[bin]  = hcs->GetYaxis()->GetBinUpEdge(by+1);
        ptmin[bin] = hcs->GetXaxis()->GetBinLowEdge(bpt+1);
        ptmax[bin] = hcs->GetXaxis()->GetBinUpEdge(bpt+1);
        bin++;
      }
    }
    return bin;
  }
  
  // Transform FORTRAN string to TString
  TString MakeTStringFromFORTRAN(char* for_str, int len)
  {
    TString str(for_str, len);
    str = TString(str.Data(), str.First(' '));
    return str;
  }
}

extern "C"
{
  // HVQMNR common block defined in hvqmnr.f
  struct COMMON_HVQMNR_Pars
  {
    double sqrt_s;
    double mc;
    double mb;
    double mf_A_c;
    double mf_B_c;
    double mf_C_c;
    double mr_A_c;
    double mr_B_c;
    double mr_C_c;
    double mf_A_b;
    double mf_B_b;
    double mf_C_b;
    double mr_A_b;
    double mr_B_b;
    double mr_C_b;
    double fragpar_c;
    double fragpar_b;
    bool debug;
  };
  extern COMMON_HVQMNR_Pars hvqmnr_pars_;
  
  
  // Calculation for LHCb 7 TeV charm (arXiv:1302.2864)
  void hvqmnr_lhcb_7tev_charm_(int* IfcnCount, char* XSecType, char* FinalState, 
                               int* nbin, double* xsec, double* ymin, double* ymax, double* ptmin, double* ptmax,
                               int len_XSecType, int len_FinalState)
  {
    // Static instances
    static int IfcnCount_Last = 0;
    static HVQMNR::MNR mnr;
    static HVQMNR::Frag frag;
    static HVQMNR::Grid grid, grid_smoothed;
    const int nfinalstate = 5;
    static TH2D hcs[nfinalstate];
    static TH2D* hcs_ptr[nfinalstate];
    
    // Transform input FORTRAN strings to ROOT TStrings
    TString xsectype = HVQMNR::MakeTStringFromFORTRAN(XSecType, len_XSecType);
    TString str_finalstate = HVQMNR::MakeTStringFromFORTRAN(FinalState, len_FinalState);

    // If this is first call at this fit iteration, perform calculation
    if(IfcnCount_Last != *IfcnCount)
      {
      // Update fcn counter
      IfcnCount_Last = *IfcnCount;
        
      // If this is first call at first iteration, perform initialisation and pre-calculation
      if(IfcnCount_Last == 1)
      {
        // Stereing parameters for this calculation (modify only if you understand what you are doing)
        HVQMNR::Steering steer;
        steer.ptmin = 0.001;
        steer.ptmax = 20.0;
        steer.npt   = 25;
        steer.nptsm = 500;
        steer.ymin  = 2.0;
        steer.ymax  = 6.0;
        steer.ny    = 40;
        steer.nsfnb = 500;
        steer.nx3   = 25;
        steer.nx4   = 125;
        steer.nbz   = 50;
          
        HVQMNR::DefaultInit(steer, hvqmnr_pars_.mc, mnr, frag, grid, grid_smoothed);
        // MNR (parton-level calculation)
        mnr.SetDebug(hvqmnr_pars_.debug);
        mnr.fC_sh = TMath::Power(hvqmnr_pars_.sqrt_s, 2.0); // centre-of-mass energy squared
        mnr.CalcConstants();
        // Fragmentation
        frag.SetDebug(hvqmnr_pars_.debug);
        // Add needed final states
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "dzero", hvqmnr_pars_.fragpar_c), HVQMNR::Frag::fM_dzero);
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "",      hvqmnr_pars_.fragpar_c), HVQMNR::Frag::fM_dstar);
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "dch",   hvqmnr_pars_.fragpar_c), HVQMNR::Frag::fM_dch);
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "",      hvqmnr_pars_.fragpar_c), HVQMNR::Frag::fM_ds);
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "",      hvqmnr_pars_.fragpar_c), HVQMNR::Frag::fM_lambdac);

        // Set binning for cross-section histograms
        int nbin_y = 5;
        double bin_y[6] = { 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 };
        int nbin_pt = 8;
        double bin_pt[9] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
        // D0, D+, D*+, D_s
        for(int f = 0; f < 4; f++) 
        {
          hcs_ptr[f] = &hcs[f];
          hcs_ptr[f]->SetBins(nbin_pt, bin_pt, nbin_y, bin_y);
        }
        // ... and Lambda_c
        if(xsectype.EqualTo("norm_y"))
        {
          int nbin_pt_lambdac = 1;
          double bin_pt_lambdac[2] = { 2.0, 8.0 };
          hcs_ptr[4] = &hcs[4];
          hcs_ptr[4]->SetBins(nbin_pt_lambdac, bin_pt_lambdac, nbin_y, bin_y);
        }
        else 
        {
          int nbin_y_lambdac = 1;
          double bin_y_lambdac[2] = { 2.0, 4.5 };
          hcs_ptr[4] = &hcs[4];
          hcs_ptr[4]->SetBins(nbin_pt, bin_pt, nbin_y_lambdac, bin_y_lambdac);
        }
      }

      // Update parameters and perform calculation
      mnr.SetScaleCoef(hvqmnr_pars_.mf_A_c, hvqmnr_pars_.mf_B_c, hvqmnr_pars_.mf_C_c, hvqmnr_pars_.mr_A_c, hvqmnr_pars_.mr_B_c, hvqmnr_pars_.mr_C_c);
      mnr.CalcXS(&grid, hvqmnr_pars_.mc);
      HVQMNR::Grid::InterpolateGrid(&grid, &grid_smoothed, hvqmnr_pars_.mc);
      for(int f = 0; f < nfinalstate; f++) frag.GetFF(f)->SetParameter(1, hvqmnr_pars_.fragpar_c);
      frag.CalcCS(&grid_smoothed, hvqmnr_pars_.mc, hcs_ptr);
    }

    // Get histogram corresponding for final state cross-section and fill output arrays
    TH1* finalstate_hcs = NULL;
    if     (str_finalstate == "dzero")    finalstate_hcs = hcs_ptr[0];
    else if(str_finalstate == "dstar")    finalstate_hcs = hcs_ptr[1];
    else if(str_finalstate == "dch")      finalstate_hcs = hcs_ptr[2];
    else if(str_finalstate == "ds")       finalstate_hcs = hcs_ptr[3];
    else if(str_finalstate == "lambdac")  finalstate_hcs = hcs_ptr[4];
    else
    {
      printf("ERROR in hvqmnr_lhcb_7tev_charm(): unknown FinalState %s\n", str_finalstate.Data());
      hf_stop_();
    }
    *nbin = HVQMNR::FillOutput(finalstate_hcs, xsec, ymin, ymax, ptmin, ptmax);
  }


  // Calculation for LHCb 7 TeV beauty (arXiv:1306.3663)
  void hvqmnr_lhcb_7tev_beauty_(int* IfcnCount, char* FinalState, 
                               int* nbin, double* xsec, double* ymin, double* ymax, double* ptmin, double* ptmax,
                               int len_FinalState)
  {
    // Static instances
    static int IfcnCount_Last = 0;
    static HVQMNR::MNR mnr;
    static HVQMNR::Frag frag;
    static HVQMNR::Grid grid, grid_smoothed;
    const int nfinalstate = 3;
    static TH2D hcs[nfinalstate];
    static TH2D* hcs_ptr[nfinalstate];
    
    // Transform input FORTRAN strings to ROOT TStrings
    TString str_finalstate = HVQMNR::MakeTStringFromFORTRAN(FinalState, len_FinalState);

    // If this is first call at this fit iteration, perform calculation
    if(IfcnCount_Last != *IfcnCount)
      {
      // Update fcn counter
      IfcnCount_Last = *IfcnCount;
        
      // If this is first call at first iteration, perform initialisation and pre-calculation
      if(IfcnCount_Last == 1)
      {
        // Stereing parameters for this calculation (modify only if you understand what you are doing)
        HVQMNR::Steering steer;
        steer.ptmin = 0.001;
        steer.ptmax = 70.0;
        steer.npt   = 35;
        steer.nptsm = 500;
        steer.ymin  = 2.0;
        steer.ymax  = 5.0;
        steer.ny    = 50;
        steer.nsfnb = 500;
        steer.nx3   = 125;
        steer.nx4   = 125;
        steer.nbz   = 100;
          
        HVQMNR::DefaultInit(steer, hvqmnr_pars_.mb, mnr, frag, grid, grid_smoothed);
        // MNR (parton-level calculation)
        mnr.SetDebug(hvqmnr_pars_.debug);
        mnr.fC_sh = TMath::Power(hvqmnr_pars_.sqrt_s, 2.0); // centre-of-mass energy squared
        mnr.CalcConstants();
        // Fragmentation
        frag.SetDebug(hvqmnr_pars_.debug);
        // Add needed final states
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "", hvqmnr_pars_.fragpar_b), HVQMNR::Frag::fM_bch);
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "", hvqmnr_pars_.fragpar_b), HVQMNR::Frag::fM_bzero);
        frag.AddOut(HVQMNR::Frag::GetFragFunction(0, "", hvqmnr_pars_.fragpar_b), HVQMNR::Frag::fM_bs);

        // Set binning for cross-section histograms
        for(int f = 0; f < nfinalstate; f++) hcs_ptr[f] = &hcs[f];
        int nbin_y = 5;
        double bin_y[6] = { 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 };
        // B+
        int nbin_pt_bch = 27;
        double bin_pt_bch[28] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 
                                  7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.5, 12.5, 14.0, 16.5, 23.5, 40.0 };
        hcs_ptr[0]->SetBins(nbin_pt_bch, bin_pt_bch, nbin_y, bin_y);
        // B0
        int nbin_pt_bzero = 19;
        double bin_pt_bzero[20] = { 0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 
                                    7.0, 7.5, 8.0, 9.0, 10.0, 11.5, 14.0, 19.5, 40.0 };
        hcs_ptr[1]->SetBins(nbin_pt_bzero, bin_pt_bzero, nbin_y, bin_y);
        // Bs
        int nbin_pt_bs = 15;
        double bin_pt_bs[16] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
                                 9.0, 10.0, 11.0, 13.0, 15.0, 19.0, 40.0 };
        hcs_ptr[2]->SetBins(nbin_pt_bs, bin_pt_bs, nbin_y, bin_y);
      }

      // Update parameters and perform calculation
      mnr.SetScaleCoef(hvqmnr_pars_.mf_A_b, hvqmnr_pars_.mf_B_b, hvqmnr_pars_.mf_C_b, hvqmnr_pars_.mr_A_b, hvqmnr_pars_.mr_B_b, hvqmnr_pars_.mr_C_b);
      mnr.CalcXS(&grid, hvqmnr_pars_.mb);
      HVQMNR::Grid::InterpolateGrid(&grid, &grid_smoothed, hvqmnr_pars_.mb);
      for(int f = 0; f < nfinalstate; f++) frag.GetFF(f)->SetParameter(1, hvqmnr_pars_.fragpar_b);
      frag.CalcCS(&grid_smoothed, hvqmnr_pars_.mb, hcs_ptr);
    }

    // Get corresponding for final state cross-section histogram and fill output arrays
    TH1* finalstate_hcs = NULL;
    if     (str_finalstate == "bch")   finalstate_hcs = hcs_ptr[0];
    else if(str_finalstate == "bzero") finalstate_hcs = hcs_ptr[1];
    else if(str_finalstate == "bs")    finalstate_hcs = hcs_ptr[2];
    else
    {
      printf("ERROR in hvqmnr_lhcb_7tev_beauty(): unknown FinalState %s\n", str_finalstate.Data());
      hf_stop_();
    }
    *nbin = HVQMNR::FillOutput(finalstate_hcs, xsec, ymin, ymax, ptmin, ptmax);
  }
}
