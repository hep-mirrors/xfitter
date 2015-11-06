#include <fortran_wrapper.h>
#include <hvqmnr_grid.h>
#include <hvqmnr_mnr.h>
#include <TMath.h>
#include <TF1.h>


namespace HVQMNR
{
  // Interface to original MNR routines (hvqcrsx.f)
  extern "C" {
    // gg
    double hqh0gg_(double* t1, double* rho);
    double hqh0gg_opt_(double* t1, double* rho);
    double hqhpgg_(double* tx, double* t1, double* rho);
    double hqbpgg_(double* tx, double* t1, double* rho);
    double hqhlgg_(double* tx, double* t1, double* rho);
    double hqhdgg_(double* t1, double* rho, int* nl);
    double hqbdgg_(double* t1, double* rho, int* nl);
    // qa
    double hqh0qa_(double* t1, double* rho);
    double hqhpqa_(double* tx, double* t1, double* rho);
    double hqbpqa_(double* tx, double* t1, double* rho);
    double hqhlqa_(double* tx, double* t1, double* rho);
    double hqhdqa_(double* t1, double* rho, int* nl);
    double hqbdqa_(double* t1, double* rho, int* nl);
    // ... and asym.
    double ashpqa_(double* tx, double* t1, double* ro);
    double ashdqa_(double* t1t, double* rot, int* nl);
    // gq
    double hqhpqg_(double* tx, double* t1, double* rho);
    double hqbpqg_(double* tx, double* t1, double* rho);
    double hqhlqg_(double* tx, double* t1, double* rho);
    // ... and asym.
    double ashpqg_(double* tx, double* t1, double* ro);
    double ashdqg_(double* t1t, double* rot, int* nl);
  }

  MNR::MNR()
  {
    fSF_pdf = NULL;
    fC_mem  = NULL;
    fBc_x3  = NULL;
    fBw_x3  = NULL;
    fBc_x4  = NULL;
    fBw_x4  = NULL;
    bDebug = false;
    bFirst = true;
    fNRecalc = 0;
  }

  MNR::~MNR() 
  {
    if(fSF_pdf) delete fSF_pdf;
    if(fC_mem)  delete fC_mem;
    if(fBc_x3)  delete fBc_x3;
    if(fBw_x3)  delete fBw_x3;
    if(fBc_x4)  delete fBc_x4;
    if(fBw_x4)  delete fBw_x4;
  }

  void MNR::CalcConstants() 
  {
    fC_sqrt_sh = TMath::Sqrt(fC_sh);
    fC_xnorm = fC_hc2 * fC_pi / fC_sh;
    fC_b0 = (11 * fC_vca - 4 * fC_vtf * fC_nl) / (12 * fC_pi);
  }

  void MNR::CalcBinning() 
  {
    // PDFs in x for certain factorisation scale
    fSF_pdf = new double[fSF_nb * fSF_npart];
    fSF_log10_min_x = TMath::Log10(fSF_min_x);
    fSF_log10_max_x = TMath::Log10(fSF_max_x);
    fSF_step_log10_x = TMath::Log10(fSF_max_x / fSF_min_x) / (fSF_nb - 1);
    fSF_min_adoptx = (fSF_log10_min_x - fSF_log10_min_x) / fSF_step_log10_x;
    fSF_max_adoptx = (fSF_log10_max_x - fSF_log10_min_x) / fSF_step_log10_x;

    // Rapidity in parton centre-of-mass system
    double bb_x3[fBn_x3+1];
    fBc_x3 = new double[fBn_x3];
    fBw_x3 = new double[fBn_x3];
    double minx3 = -6.0;
    double maxx3 = 6.0;
    double stepx3 = (maxx3 - minx3) / fBn_x3;
    for(int b = 0; b < fBn_x3 + 1; b++) 
    {
      bb_x3[b] = minx3 + stepx3 * b;
      if(b == 0) continue;
      fBc_x3[b-1] = (bb_x3[b] + bb_x3[b-1]) / 2;
      fBw_x3[b-1] =  bb_x3[b] - bb_x3[b-1];
    }

    // t3 (3 body variable): use specially optimised binning
    double fbx4_12 = 0.7;
    TF1 fbx4_1("fbx4_1", "0.973 * TMath::Power(x / 0.7, 2)", 0.0, fbx4_12);
    TF1 fbx4_2("fbx4_1", "1 + TMath::Power(x - 1, 3)", fbx4_12, 1.0);
    double bb_x4[fBn_x4+1];
    fBc_x4 = new double[fBn_x4];
    fBw_x4 = new double[fBn_x4];
    for(int b = 0; b < fBn_x4 + 1; b++) 
    {
      bb_x4[b] = 1. * b / fBn_x4;
      if(bb_x4[b] <= fbx4_12) bb_x4[b] = fbx4_1.Eval(bb_x4[b]);
      else bb_x4[b] = fbx4_2.Eval(bb_x4[b]);
      if(b == 0) continue;
      fBc_x4[b-1] = (bb_x4[b] + bb_x4[b-1]) / 2;
      fBw_x4[b-1] = bb_x4[b] - bb_x4[b-1];
    }
  }

  void MNR::SetScaleCoef(double mf_a, double mf_b, double mf_c, double mr_a, double mr_b, double mr_c) 
  {
    // Check for possible nan 
    // (it happens if corresponding parameters were not provided in steering.txt via ExtraMinimisationParameters)
    if(mf_a != mf_a || mf_b != mf_b || mf_c != mf_c || mr_a != mr_a || mr_b != mr_b || mr_c != mr_c)
    {
      printf("ERROR in MNR::SetScaleCoef(): some parameters are nan\n");
      printf("Make sure that you provide all needed parameters in steering.txt (ExtraMinimisationParameters)\n");
      hf_stop_();
    }
    fMf_A = mf_a;
    fMf_B = mf_b;
    fMf_C = mf_c;
    fMr_A = mr_a;
    fMr_B = mr_b;
    fMr_C = mr_c;
  }

  double MNR::GetMf2(double xm2, double pt2) 
  {
    return fMf_A * pt2 + fMf_B * xm2 + fMf_C;
  }

  double MNR::GetMr2(double xm2, double pt2) 
  {
    return fMr_A * pt2 + fMr_B * xm2 + fMr_C;
  }

  void MNR::PrecalculatePDF(double mf2) 
  {
    if(mf2 < fSF_min_mf2 || mf2 > fSF_max_mf2) 
    {
      printf("WARNING in MNR::PrecalculatePDF(): mf2 %e out of range %e .. %e\n", mf2, fSF_min_mf2, fSF_max_mf2);
      printf("PDFs are set to 0\n");
      for(int b = 0; b < fSF_nb; b++) 
        for(int nf = -fC_nl; nf <= fC_nl; nf++) fSF_pdf[b*fSF_npart+nf] = 0.0;
      return;
    }
    for(int b = 0; b < fSF_nb; b++) 
    {
      double log10x = fSF_log10_min_x + b * fSF_step_log10_x;
      this->GetPDF(mf2, TMath::Power(10., log10x), fSF_pdf + fSF_npart * b);
    }
  }

  int MNR::GetSF(double& pdf_gg, double& pdf_qq, double& pdf_qq_a, double& pdf_qg, double& pdf_qg_a, 
                 double& pdf_qg_r, double& pdf_qg_a_r, double adoptx1, double adoptx2, double mf2/* = -1.0*/) 
  {
    pdf_gg = pdf_qq = pdf_qq_a = pdf_qg = pdf_qg_a = pdf_qg_r = pdf_qg_a_r = 0.0;
    if(mf2 > 0.0 && (mf2 < fSF_min_mf2 || mf2 > fSF_max_mf2)) 
    {
      printf("WARNING in MNR::GetSF(): mf2 %e out of range %e .. %e\n", mf2, fSF_min_mf2, fSF_max_mf2);
      return 1;
    }
    if(adoptx1 > fSF_max_adoptx || adoptx2 > fSF_max_adoptx) return 1;
    if(adoptx1<fSF_min_adoptx) adoptx1 = fSF_min_adoptx;
    if(adoptx2<fSF_min_adoptx) adoptx2 = fSF_min_adoptx;
    double pdf1[fSF_npart], pdf2[fSF_npart];
    if(mf2 < 0.0) 
    {
      double part1, part2;
      double delta1 = modf(adoptx1, &part1);
      double delta2 = modf(adoptx2, &part2);
      double one_m_delta1 = 1.0 - delta1;
      double one_m_delta2 = 1.0 - delta2;
      int offset1 = int(part1) * fSF_npart + 6;
      int offset2 = int(part2) * fSF_npart + 6;
      int one_p_offset1 = offset1 + fSF_npart;
      int one_p_offset2 = offset2 + fSF_npart;
      for(int nf = -fC_nl; nf <= fC_nl; nf++) 
      {
        int six_p_nf = 6 + nf;
        pdf1[six_p_nf] = fSF_pdf[offset1+nf] * one_m_delta1 + fSF_pdf[one_p_offset1+nf] * delta1;
        pdf2[six_p_nf] = fSF_pdf[offset2+nf] * one_m_delta2 + fSF_pdf[one_p_offset2+nf] * delta2;
      }
    }
    else 
    {
      double x1 = TMath::Power(10.0, fSF_step_log10_x * adoptx1 + fSF_log10_min_x);
      double x2 = TMath::Power(10.0, fSF_step_log10_x * adoptx2 + fSF_log10_min_x);
      this->GetPDF(mf2, x1, pdf1);
      this->GetPDF(mf2, x2, pdf2);
    }
    // Calculate gg
    pdf_gg = pdf1[6] * pdf2[6];
    // ... and the rest
    for(int nf = 1; nf <= fC_nl; nf++) 
    {
      int six_p_nf = 6 + nf;
      int six_m_nf = 6 - nf;
      double qa = pdf1[six_p_nf] * pdf2[six_m_nf]; // h1_p1(x1) x h2_p2(x2)
      double aq = pdf1[six_m_nf] * pdf2[six_p_nf]; // h1_a1(x1) x h2_a2(x2)
      double qg = pdf1[six_p_nf] * pdf2[6]; // h1_p1(x1) x h2_p2(x2)
      double ag = pdf1[six_m_nf] * pdf2[6]; // h1_a1(x1) x h2_a2(x2)
      double gq = pdf1[6] * pdf2[six_p_nf]; // h1_p2(x1) x h2_p1(x2)
      double ga = pdf1[6] * pdf2[six_m_nf]; // h1_a2(x1) x h2_a1(x2)
      // qq
      pdf_qq      +=  qa + aq;
      pdf_qq_a    +=  qa - aq;
      // qg
      pdf_qg      +=  qg + ag;
      pdf_qg_a    +=  qg - ag;
      pdf_qg_r    +=  gq + ga;
      pdf_qg_a_r  +=  gq - ga;
    }
    return 0;
  }

  void MNR::GetPDF(double mf2, double x, double pdf[13]) 
  {
    hf_get_pdfs_(&x, &mf2, pdf);
  }

  double MNR::GetAs(double mr2) 
  {
    return hf_get_alphas_(&mr2);
  }

  void MNR::Precalc(Grid* grid) 
  {
    fNRecalc++;
    printf("MNR::Precalc(): recalculation NO %d\n", fNRecalc);
    int n_l = grid->NL();
    double* p_l = grid->LPtr();
    int n3 = n_l * fBn_x3 * fBn_x4;
    int n2 = n_l * fBn_x3;
    // Calculate required memory (in MB)
    int ndouble = 17 * n2 + 16 * n3 + 2 * fBn_x3 + 6 * fBn_x4;
    double mb = sizeof(double) * ndouble / (1024. * 1024.);
    printf("MNR::Precalc(): required %.0f MB\n", mb);

    // Allocate memory in one place, because 
    // (1) one call to new is faster than multiple calls
    // (2) it is faster to access memory allocated in one place
    fC_mem = new double[ndouble];
    // Now go through all needed variables and set their pointers to corresponding offsets
    int mem_offset = 0;
    // LO gg
    fCh0_hqh0gg = fC_mem + mem_offset;
    mem_offset += n2;
    // LO qq
    fCh0_hqh0qa = fC_mem + mem_offset;
    mem_offset += n2;
    // NLO gg
    fCh3_hqhpgg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_hqbpgg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_hqhlgg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh2_hqhdgg = fC_mem + mem_offset;
    mem_offset += n2;
    fCh2_hqbdgg = fC_mem + mem_offset;
    mem_offset += n2;
    fCh2_hqh0gg = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_hqhpgg = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_hqbpgg = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_hqhlgg = fC_mem + mem_offset;
    mem_offset += n2;
    // NLO qq
    fCh3_hqhpqa = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_hqbpqa = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_hqhlqa = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_a_ashpqa = fC_mem + mem_offset;
    mem_offset += n3;
    fCh2_hqhdqa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh2_hqbdqa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh2_hqh0qa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh2_a_ashdqa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_hqhpqa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_hqbpqa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_hqhlqa = fC_mem + mem_offset;
    mem_offset += n2;
    fCh3c_a_ashpqa = fC_mem + mem_offset;
    mem_offset += n2;
    // NLO qg
    fCh3_hqhpqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_hqbpqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_hqhlqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_r_hqhpqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_r_hqbpqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_r_hqhlqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_a_ashpqg = fC_mem + mem_offset;
    mem_offset += n3;
    fCh3_a_r_ashpqg = fC_mem + mem_offset;
    mem_offset += n3;
    // normalisation coefficients
    fC_N = fC_mem + mem_offset;
    mem_offset += n2;
    fC_NN = fC_mem + mem_offset;
    mem_offset += n3;
    // kinematics
    fCk_t32 = fC_mem + mem_offset;
    mem_offset += fBn_x4;
    fCk_t33 = fC_mem + mem_offset;
    mem_offset += fBn_x4;
    fCk_tx = fC_mem + mem_offset;
    mem_offset += fBn_x4;
    fCk_lntx = fC_mem + mem_offset;
    mem_offset += fBn_x4;
    fCk_lntx_o_tx = fC_mem + mem_offset;
    mem_offset += fBn_x4;
    fCk_pxtcor = fC_mem + mem_offset;
    mem_offset += fBn_x4;
    fCk_t1t = fC_mem + mem_offset;
    mem_offset += fBn_x3;
    fCk_chyprim2 = fC_mem + mem_offset;
    mem_offset += fBn_x3;

    // Calculate x3 (binning in parton CMS rapidity) variables
    for(int c_x3 = 0; c_x3 < fBn_x3; c_x3++) 
    {
      double yprim = fBc_x3[c_x3];
      fCk_t1t[c_x3] = 0.5 * (1 - TMath::TanH(yprim));
      double chyprim = TMath::CosH(yprim);
      fCk_chyprim2[c_x3] = chyprim * chyprim;
    }
    fCk_sum_o_tx = 0.0;
    fCk_sum_lntx_o_tx = 0.0;

    // Calculate x4 (t3, 3 body variable) variables
    for(int c_x4 = 0; c_x4 < fBn_x4; c_x4++) 
    {
      double t3 = fBc_x4[c_x4];
      double t32 = t3 * t3;
      double t33 = t32 * t3;
      fCk_t32[c_x4] = t32;
      fCk_t33[c_x4] = t33;
      double tx = 1 - t3;
      double lntx = TMath::Log(tx);
      fCk_tx[c_x4] = tx;
      fCk_lntx[c_x4] = lntx;
      fCk_lntx_o_tx[c_x4] = lntx / tx;
      fCk_pxtcor[c_x4] = TMath::Log10(t3) / fSF_step_log10_x;
      fCk_sum_o_tx += fBw_x4[c_x4] / tx;
      fCk_sum_lntx_o_tx += fBw_x4[c_x4] * lntx/tx;
    }

    // Loop over L = xm^2 / (xm^2 + pT^2) bins
    for(int c_l = 0; c_l < n_l; c_l++) 
    {
      if(bDebug) if(c_l % 10 == 0) printf("MNR::Precalc(): 1st dimension: %3d from %3d\n", c_l, n_l);
      double l2 = p_l[c_l];
      // Loop over x3 bins
      for(int c_x3 = 0; c_x3 < fBn_x3; c_x3++) 
      {
        double yprim = fBc_x3[c_x3];
        double t1t = 0.5 * (1 - TMath::TanH(yprim));
        double chyprim = TMath::CosH(yprim);
        double chyprim2 = chyprim * chyprim;
        double rot = l2 / chyprim2;
        double tz = 0.0;
        int n2 = c_l * fBn_x3 + c_x3;
        fC_N[n2] = fC_xnorm * (1. / chyprim2) * fBw_x3[c_x3];
        // LO gg
        fCh0_hqh0gg[n2] = hqh0gg_opt_(&t1t, &rot);
        // LO qq
        fCh0_hqh0qa[n2] = hqh0qa_(&t1t, &rot);
        // NLO gg
        fCh2_hqhdgg[n2] = hqhdgg_(&t1t, &rot, &fC_nl);
        fCh2_hqbdgg[n2] = hqbdgg_(&t1t, &rot, &fC_nl);
        fCh2_hqh0gg[n2] = hqh0gg_(&t1t, &rot);
        fCh3c_hqhpgg[n2] = hqhpgg_(&tz, &t1t, &rot);
        fCh3c_hqbpgg[n2] = hqbpgg_(&tz, &t1t, &rot);
        fCh3c_hqhlgg[n2] = hqhlgg_(&tz, &t1t, &rot);
        // NLO qq
        fCh2_hqhdqa[n2] = hqhdqa_(&t1t, &rot, &fC_nl);
        fCh2_hqbdqa[n2] = hqbdqa_(&t1t, &rot, &fC_nl);
        fCh2_hqh0qa[n2] = hqh0qa_(&t1t, &rot);
        fCh2_a_ashdqa[n2] = ashdqa_(&t1t, &rot, &fC_nl);
        fCh3c_hqhpqa[n2] = hqhpqa_(&tz, &t1t, &rot);
        fCh3c_hqbpqa[n2] = hqbpqa_(&tz, &t1t, &rot);
        fCh3c_hqhlqa[n2] = hqhlqa_(&tz, &t1t, &rot);
        fCh3c_a_ashpqa[n2] = ashpqa_(&tz, &t1t, &rot);
        // Loop over x4 bins
        for(int c_x4 = 0; c_x4 < fBn_x4; c_x4++) 
        {
          double t3 = fBc_x4[c_x4];
          double t32 = t3 * t3;
          double t33 = t32 * t3;
          double tx = 1 - t3;
          double t1 = t1t * t3;
          double ro = rot * t32;
          int n3 = c_l * fBn_x3 * fBn_x4 + c_x3 * fBn_x4 + c_x4;
          // NLO n
          fC_NN[n3] = fC_N[n2] * t33 * fBw_x4[c_x4] / fC_2pi;
          // NLO gg
          fCh3_hqhpgg[n3] = hqhpgg_(&tx, &t1, &ro);
          fCh3_hqbpgg[n3] = hqbpgg_(&tx, &t1, &ro);
          fCh3_hqhlgg[n3] = hqhlgg_(&tx, &t1, &ro);
          // NLO qq
          fCh3_hqhpqa[n3] = hqhpqa_(&tx, &t1, &ro);
          fCh3_hqbpqa[n3] = hqbpqa_(&tx, &t1, &ro);
          fCh3_hqhlqa[n3] = hqhlqa_(&tx, &t1, &ro);
          fCh3_a_ashpqa[n3] = ashpqa_(&tx, &t1, &ro);
          // NLO qg
          double t1_r = t3 - t1;
          fCh3_hqhpqg[n3] = hqhpqg_(&tx, &t1, &ro);
          fCh3_hqbpqg[n3] = hqbpqg_(&tx, &t1, &ro);
          fCh3_hqhlqg[n3] = hqhlqg_(&tx, &t1, &ro);
          fCh3_r_hqhpqg[n3] = hqhpqg_(&tx, &t1_r, &ro);
          fCh3_r_hqbpqg[n3] = hqbpqg_(&tx, &t1_r, &ro);
          fCh3_r_hqhlqg[n3] = hqhlqg_(&tx, &t1_r, &ro);
          fCh3_a_ashpqg[n3] = ashpqg_(&tx, &t1, &ro);
          fCh3_a_r_ashpqg[n3] = ashpqg_(&tx, &t1_r, &ro);
        }
      }
    }
  }

  void MNR::CalcXS(Grid* grid, double xm) 
  {
    // First call: precalculate variables
    if(bFirst) {
      this->Precalc(grid);
      bFirst = false;
    }
    grid->Zero();
    
    // Check heavy-quark mass for nan
    if(xm != xm) {
      grid->NonPhys();
      return;
    }
    
    int ncontr = grid->GetNContr();
    double xm2 = xm * xm;
    int n_l = grid->NL();
    double* p_l = grid->LPtr();
    int n_y = grid->NY();
    double* p_y = grid->YPtr();
    int nbw = grid->NW();
    // Loop over pT (internally L)
    for(int c_l = 0; c_l < n_l; c_l++) 
    {
      if(bDebug) if( c_l % 10 == 0) printf("MNR::CalcXS(): 1st dimension: %3d from %3d\n", c_l, n_l);
      double l2 = p_l[c_l];
      double mt2 = xm2 / l2;
      double pt2 = mt2 - xm2;
      double pt = TMath::Sqrt(pt2);
      // Factorisation scale
      double mf2 = this->GetMf2(xm2, pt2);
      if(mf2 < fSF_min_mf2 || mf2 > fSF_max_mf2) 
      {
        grid->NonPhys(c_l);
        continue;
      }
      // Precalculate PDFs at the calculated factorisation scale
      this->PrecalculatePDF(mf2);
      // Renormalisation scale
      double mr2 = this->GetMr2(xm2, pt2);
      if(mr2 <= 0.0) 
      {
        grid->NonPhys(c_l);
        continue;
      }
      // Strong coupling and its powers
      double as = GetAs(mr2);
      double as2 = as * as;
      double as3 = as2 * as;
      // Ratios of scales
      double xmf = TMath::Log(mf2 / xm2);
      double xmr = 4 * fC_pi * fC_b0 * TMath::Log(mr2 / mf2);
      // Loop over rapidity
      for(int c_y = 0; c_y < n_y; c_y++) 
      {
        double y = p_y[c_y];
        // Loop over rapidity in parton CMS
        for(int c_x3 = 0; c_x3 < fBn_x3; c_x3++) 
        {
          double yprim = fBc_x3[c_x3];
          double chyprim2 = fCk_chyprim2[c_x3];
          double ycm = y - yprim;
          double eycm = TMath::Exp(ycm);
          double taut = 4 * chyprim2 * mt2 / fC_sh;
          double rtaut = TMath::Sqrt(taut);
          double x1t = rtaut * eycm;
          if(x1t > 1.0) continue; // outside allowed kinematic region
          double x2t = rtaut / eycm;
          if(x2t > 1.0) continue; // outside allowed kinematic region
          int n2 = c_l * fBn_x3 + c_x3;
          double kinN = (pt / mt2) * (1 / taut);
          double N = fC_N[n2] * as2 * kinN;
          double px1t = (TMath::Log10(x1t) - fSF_log10_min_x) / fSF_step_log10_x;
          double px2t = (TMath::Log10(x2t) - fSF_log10_min_x) / fSF_step_log10_x;
          // *** LO ***
          double pdf_gg_t, pdf_qq_t, pdf_qq_a_t, pdf_qg_t, pdf_qg_a_t, pdf_qg_r_t, pdf_qg_a_r_t;
          int kinok_t = this->GetSF(pdf_gg_t, pdf_qq_t, pdf_qq_a_t, pdf_qg_t, pdf_qg_a_t, pdf_qg_r_t, pdf_qg_a_r_t, px1t, px2t);
          if(kinok_t) continue;
          // gg
          double w_lo_gg = N * fCh0_hqh0gg[n2] * pdf_gg_t;
          // qq
          double w_lo_qq = N * fCh0_hqh0qa[n2] * pdf_qq_t;
          // X-section
          int bw_t = 0;
          if(nbw != 1) 
          {
            double what_t = taut * fC_sh - 4 * xm2;
            bw_t = grid->FindWBin(what_t);
          }
          // Add all LO contributions
          for(int c = 0; c < ncontr; c++) 
          {
            double& cs = grid->CS(c, c_l, c_y, bw_t);
            MNRContribution* contr = grid->GetContr(c);
            if(contr->fActive == 0) continue;
            if(contr->fLO) 
            {
              if(contr->fgg) cs += w_lo_gg;
              if(contr->fqq) cs += w_lo_qq;
            }
          }
          // *** NLO (CE) *** (counter events)
          double NN = fC_N[n2] * as3 * kinN / fC_2pi;
          //gg
          double me_gg_h2   =  (fCh2_hqhdgg[n2] + xmf * fCh2_hqbdgg[n2] + xmr * fCh2_hqh0gg[n2]);
          double me_gg_h3c  =  (fCh3c_hqhpgg[n2] + xmf * fCh3c_hqbpgg[n2]) * fCk_sum_o_tx + fCh3c_hqhlgg[n2] * fCk_sum_lntx_o_tx;
          double w_nlo_gg_c = NN*(me_gg_h2-me_gg_h3c)*pdf_gg_t;
          //qq
          double me_qq_h2    =  (fCh2_hqhdqa[n2] + xmf * fCh2_hqbdqa[n2] + xmr * fCh2_hqh0qa[n2]);
          double me_qq_h2_a  =  (bFS_Q&&bFS_A) ? 0.0 : fCh2_a_ashdqa[n2];
          double me_qq_h3c   =  (fCh3c_hqhpqa[n2] + xmf * fCh3c_hqbpqa[n2]) * fCk_sum_o_tx + fCh3c_hqhlqa[n2] * fCk_sum_lntx_o_tx;
          double me_qq_h3c_a =  (bFS_Q&&bFS_A) ? 0.0 : fCh3c_a_ashpqa[n2] * fCk_sum_o_tx;
          if(bFS_A&&!bFS_Q) pdf_qq_a_t *= -1;
          double w_nlo_qq_c = NN * ((me_qq_h2 - me_qq_h3c) * pdf_qq_t + (me_qq_h2_a - me_qq_h3c_a) * pdf_qq_a_t);
          // X-section
          for(int c = 0; c < ncontr; c++) 
          {
            double& cs_c = grid->CS(c, c_l, c_y, bw_t);
            MNRContribution* contr = grid->GetContr(c);
            if(contr->fActive == 0) continue;
            if(contr->fNLO) {
              if(contr->fgg) cs_c += w_nlo_gg_c;
              if(contr->fqq) cs_c += w_nlo_qq_c;
            }
          }
          // Loop over t3 (3 body variable)
          for(int c_x4 = fBn_x4 - 1; c_x4 >= 0; c_x4--) 
          {
            double px1 = px1t - fCk_pxtcor[c_x4];
            double px2 = px2t - fCk_pxtcor[c_x4];
            double pdf_gg, pdf_qq, pdf_qq_a, pdf_qg, pdf_qg_a, pdf_qg_r, pdf_qg_a_r;
            int kinok = this->GetSF(pdf_gg, pdf_qq, pdf_qq_a, pdf_qg, pdf_qg_a, pdf_qg_r, pdf_qg_a_r, px1, px2);
            if(kinok) break;
            int n3 = c_l * fBn_x3 * fBn_x4 + c_x3 * fBn_x4 + c_x4;
            double t32 = fCk_t32[c_x4];
            double NN = fC_NN[n3] * as3 * kinN;
            double tx = fCk_tx[c_x4];
            double lntx_o_tx = fCk_lntx_o_tx[c_x4];
            // gg
            double me_gg_h3   =  kinok ? 0.0 : (fCh3_hqhpgg[n3] + xmf * fCh3_hqbpgg[n3]) / tx + fCh3_hqhlgg[n3] * lntx_o_tx;
            double w_nlo_gg_d = NN * me_gg_h3 * pdf_gg;
            // qq
            double me_qq_h3   =  kinok ? 0.0 : (fCh3_hqhpqa[n3] + xmf * fCh3_hqbpqa[n3]) / tx + fCh3_hqhlqa[n3] * lntx_o_tx;
            double me_qq_h3_a =  (kinok || (bFS_Q && bFS_A)) ? 0.0 : fCh3_a_ashpqa[n3] / tx;
            if(bFS_A && !bFS_Q) pdf_qq_a *= -1;
            double w_nlo_qq_d = NN * (me_qq_h3 * pdf_qq + me_qq_h3_a * pdf_qq_a);
            double me_qg_h3     =  kinok ? 0.0 : (fCh3_hqhpqg[n3] + xmf * fCh3_hqbpqg[n3]) / tx + fCh3_hqhlqg[n3] * lntx_o_tx;
            double me_qg_h3_r   =  kinok ? 0.0 : (fCh3_r_hqhpqg[n3] + xmf * fCh3_r_hqbpqg[n3]) / tx + fCh3_r_hqhlqg[n3] * lntx_o_tx;
            double me_qg_h3_a   =  (kinok || (bFS_Q && bFS_A)) ? 0.0 : fCh3_a_ashpqg[n3] / tx;
            double me_qg_h3_a_r =  (kinok || (bFS_Q && bFS_A)) ? 0.0 : fCh3_a_r_ashpqg[n3] / tx;
            if(bFS_A && !bFS_Q) 
            {
              pdf_qg_a *= -1;
              pdf_qg_a_r *= -1;
            }
            double w_nlo_qg_d = NN * (me_qg_h3 * pdf_qg + me_qg_h3_a * pdf_qg_a + me_qg_h3_r * pdf_qg_r + me_qg_h3_a_r * pdf_qg_a_r);
            // X-section
            int bw = 0;
            // Determine W bins, if cross section in mupltiple W bins is needed
            if(nbw != 1) 
            {
              double tau = taut / t32;
              double what = tau * fC_sh - 4 * xm2;
              bw = grid->FindWBin(what);
            }
            // Add all NLO contributions
            for(int c = 0; c < ncontr; c++) 
            {
              double& cs = grid->CS(c, c_l, c_y, bw);
              MNRContribution* contr = grid->GetContr(c);
              if(contr->fActive == 0) continue;
              if(contr->fNLO) {
                if(contr->fgg) cs += w_nlo_gg_d;
                if(contr->fqq) cs += w_nlo_qq_d;
                if(contr->fqg) cs += w_nlo_qg_d;
              }
            }
          }
        }
      }
    }
    // Multiply all cross sections by 2, if both particle and particle final states are needed
    int n_w = grid->NW();
    if(bFS_Q && bFS_A)
      for(int c = 0; c < ncontr; c++)
        for(int bl = 0; bl < n_l; bl++)
          for(int by = 0; by < n_y; by++)
            for(int bw = 0; bw < n_w; bw++)
              grid->CS(c, bl, by, bw) *= 2;
  }
  
  // Constants
  const double MNR::fC_pi       = 3.14159265359e0;
  const double MNR::fC_2pi      = 6.28318530718e0;
  const double MNR::fC_hc2      = 3.8937966e+2;
  const double MNR::fC_vca      = 3.0e0;
  const double MNR::fC_vtf      = 0.5e0;
  const int    MNR::fSF_npart   = 13;
  const double MNR::fSF_min_x   = 1e-6;
  const double MNR::fSF_max_x   = 1e0;
  const double MNR::fSF_min_mf2 = 1e0;
  const double MNR::fSF_max_mf2 = 8e4;
}
