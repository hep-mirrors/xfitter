#include "EvolutionAPFELxx.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include "hf_errlog.h"

#include <apfel/alphaqcd.h>
#include <apfel/messages.h>
#include <apfel/rotations.h>
#include <apfel/evolutionbasisqcd.h>
#include <apfel/constants.h>
#include <apfel/tools.h>


namespace xfitter
{
  // the class factories
  extern "C" EvolutionAPFELxx*create(const char*name){
    return new EvolutionAPFELxx(name);
  }
  const char*EvolutionAPFELxx::getClassName()const{return "APFELxx";}


  //_________________________________________________________________________________
  // Moments {K_ns(1), K_qq(2), K_qg(2), K_gq(2), K_gg(2)} of the
  // scheme-change kernels at nf, obtained by applying the PDF
  // transformation to test sets with only valence, singlet or gluon.
  // A dedicated fine grid makes them accurate to ~1e-5.
  std::array<double, 5> SchemeChangeMoments(std::vector<double> const& Thresholds, int const& nf)
  {
    const apfel::Grid g{{apfel::SubGrid{200, 1e-7, 4}, apfel::SubGrid{150, 1e-1, 4}, apfel::SubGrid{150, 6e-1, 4}, apfel::SubGrid{150, 9e-1, 4}}};
    const auto K = SchemeChangeKernels(g, Thresholds);

    // int_xmin^1 dx x^(N-2) d(x) with Simpson in ln(x)
    const double xmin = g.GetJointGrid().GetGrid()[0];
    const auto Moment = [=] (apfel::Distribution const& d, int const& N) -> double
    {
      const int n = 20000;
      const double h = - log(xmin) / n;
      double s = 0;
      for (int i = 0; i <= n; i++)
        {
          const double x = std::min(xmin * exp(i * h), 1 - 1e-12);
          s += ( i == 0 || i == n ? 1 : ( i % 2 ? 4 : 2 ) ) * pow(x, N - 1) * d.Evaluate(x);
        }
      return s * h / 3;
    };

    // Shift (as/4pi = 1) of a test set in which the objects in "objs"
    // equal D and all others vanish. Heavy T's and V's follow the
    // singlet and the total valence, as for vanishing heavy quarks.
    const apfel::Distribution D{g, [] (double const& x) -> double { return x * x * pow(1 - x, 3); }};
    const apfel::Distribution Z{g, [] (double const&) -> double { return 0; }};
    const auto Shift = [&] (std::vector<int> const& objs) -> apfel::Set<apfel::Distribution>
    {
      std::map<int, apfel::Distribution> m;
      for (int i = 0; i <= 12; i++)
        m.insert({i, (std::find(objs.begin(), objs.end(), i) != objs.end() ? D : Z)});
      const apfel::Set<apfel::Distribution> f{apfel::EvolutionBasisQCD{nf}, m};
      return apfel::ChangeFactorisationSchemeMSbarToK(f, K, nf, apfel::FourPi) - f;
    };
    using EB = apfel::EvolutionBasisQCD;
    const apfel::Set<apfel::Distribution> dV = Shift({EB::VALENCE, EB::V15, EB::V24, EB::V35});
    const apfel::Set<apfel::Distribution> dS = Shift({EB::SIGMA, EB::T15, EB::T24, EB::T35});
    const apfel::Set<apfel::Distribution> dG = Shift({EB::GLUON});
    const double D1 = Moment(D, 1);
    const double D2 = Moment(D, 2);
    return {Moment(dV.at(EB::VALENCE), 1) / D1,
            Moment(dS.at(EB::SIGMA), 2) / D2, Moment(dG.at(EB::SIGMA), 2) / D2,
            Moment(dS.at(EB::GLUON), 2) / D2, Moment(dG.at(EB::GLUON), 2) / D2};
  }

  //_________________________________________________________________________________
  void EvolutionAPFELxx::atStart()
  {
    // APFEL++ banner
    apfel::Banner();

    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    _inPDFs=XFITTER_PARS::getInputDecomposition(yamlNode);
    // Retrieve parameters needed to initialize APFEL++.
    const double* MCharm   = XFITTER_PARS::getParamD("mch");
    const double* MBottom  = XFITTER_PARS::getParamD("mbt");
    const double* MTop     = XFITTER_PARS::getParamD("mtp");
    const YAML::Node xGrid = yamlNode["xGrid"];

    vector<apfel::SubGrid> sgv;
    for(auto const& sg : xGrid)
      sgv.push_back(apfel::SubGrid{sg[0].as<int>(), sg[1].as<double>(), sg[2].as<int>()});

    // x-space grid (the grid parameters should be in parameters.yaml
    // in APFELxx/yaml. I need to find a clever way to retrieve them)
    _Grid = std::unique_ptr<const apfel::Grid>(new apfel::Grid(sgv));

    // Vectors of masses and thresholds
    _Masses = {0, 0, 0, *MCharm, *MBottom};//, *MTop};
    _Thresholds = _Masses;

    // Initialize QCD evolution objects
    //_DglapObj = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses, _Thresholds);
    // APFELxx >= cereal branch: exact 4-loop pure-singlet splitting function
    // (arXiv:2609.37563), only qg, gq, gg remain approximate.
    std::vector<int>  IMod = {0, 0, 0};
    //imod 1 : A
    //imod 2 : B
    //imod 0 : (A+B)/2
    if (yamlNode["P3SGps"] && yamlNode["P3SGps"].as<int>() != 0)
      hf_errlog(26100301, "W: P3SGps ignored, the N3LO pure-singlet splitting function is exact");
    IMod[0] = yamlNode["P3SGqg"].as<int>();
    IMod[1] = yamlNode["P3SGgq"].as<int>();
    IMod[2] = yamlNode["P3SGgg"].as<int>();
    bool n3lo = OrderMap(XFITTER_PARS::getParamS("Order"))-1 >= 3;

    //IMod[7] = yamlNode["aHg"].as<int>();
    //IMod[8] = yamlNode["aggH"].as<int>();
    //IMod[9] = yamlNode["aqqHns"].as<int>();

    _DglapObj = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses, _Thresholds, false, 1e-5, n3lo, IMod);

    _Masses3F = {0, 0, 0};
    _Thresholds3F = _Masses3F;

    _Masses4F = {0, 0, 0, *MCharm};
    _Thresholds4F = _Masses4F;
    
    _DglapObj3F = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses3F, _Thresholds3F, false, 1e-5, n3lo, IMod);
    _DglapObj4F = apfel::InitializeDglapObjectsQCD(*_Grid, _Masses4F, _Thresholds4F, false, 1e-5, n3lo, IMod);

    // Change of factorisation scheme
    const auto K = SchemeChangeKernels(*_Grid, _Thresholds);
    if (!K.empty())
      {
        _DglapObj   = apfel::ChangeFactorisationSchemeMSbarToK(_DglapObj,   K);
        _DglapObj3F = apfel::ChangeFactorisationSchemeMSbarToK(_DglapObj3F, SchemeChangeKernels(*_Grid, _Thresholds3F));
        _DglapObj4F = apfel::ChangeFactorisationSchemeMSbarToK(_DglapObj4F, SchemeChangeKernels(*_Grid, _Thresholds4F));

        // The input PDFs are f^K = (1 + a0 K) f^MSbar at Q0, with f^MSbar
        // obeying the standard sum rules: in moment space this gives the
        // factor c_ns of the valence sum rules and the weights (w_q, w_g)
        // of the momentum sum rule, (w_q, w_g) = (1, 1) (1 + a0 K(2))^-1.
        const double Q0 = *XFITTER_PARS::getParamD("Q0");
        const int nf0 = apfel::NF(Q0, _Thresholds);
        const std::array<double, 5> Km = SchemeChangeMoments(_Thresholds, nf0);
        std::cout << "[INFO] Scheme-change kernel moments (as/4pi = 1) at nf = " << nf0 << ": K_ns(1) = " << Km[0]
                  << ", K_qq(2) = " << Km[1] << ", K_qg(2) = " << Km[2] << ", K_gq(2) = " << Km[3] << ", K_gg(2) = " << Km[4] << std::endl;
        const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
        const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
        const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
        const std::vector<double> Masses = _Masses, Thresholds = _Thresholds;
        XFITTER_PARS::gSchemeSumRuleFactors = [=] () -> std::array<double, 3>
        {
          const double a0  = apfel::AlphaQCD{*Alphas_ref, *Q_ref, Masses, Thresholds, PtOrder}.Evaluate(Q0) / apfel::FourPi;
          const double det = ( 1 + a0 * Km[1] ) * ( 1 + a0 * Km[4] ) - a0 * a0 * Km[2] * Km[3];
          return {1 + a0 * Km[0], ( 1 + a0 * ( Km[4] - Km[3] ) ) / det, ( 1 + a0 * ( Km[1] - Km[2] ) ) / det};
        };
        const std::array<double, 3> sr = XFITTER_PARS::gSchemeSumRuleFactors();
        std::cout << "[INFO] Sum-rule factors at start: c_ns = " << sr[0] << ", w_q = " << sr[1] << ", w_g = " << sr[2] << std::endl;
      }

    atConfigurationChange();
  }

  //_________________________________________________________________________________
  void EvolutionAPFELxx::atIteration()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    // Retrieve the relevant parameters needed to compute the evolutions
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
    const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
    const YAML::Node QGrid   = yamlNode["QGrid"];

    // Reinitialise and tabulate the running coupling at every
    // iteration. This is fast enough and allows for the reference
    // value to be fitted.
    apfel::AlphaQCD a{*Alphas_ref, *Q_ref, _Masses, _Thresholds, PtOrder};
    const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
    _AlphaQCD = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

    // 3F and 4F couplings, matched to the VFNS coupling at Q0 (assumes Q0 < mc)
    apfel::AlphaQCD a3F{_AlphaQCD(*Q0), *Q0, _Masses3F, _Thresholds3F, PtOrder};
    const apfel::TabulateObject<double> Alphas3F{a3F, 100, 0.9, 1001, 3};
    _AlphaQCD3F = [=] (double const& mu) -> double{ return Alphas3F.Evaluate(mu); };

    apfel::AlphaQCD a4F{_AlphaQCD(*Q0), *Q0, _Masses4F, _Thresholds4F, PtOrder};
    const apfel::TabulateObject<double> Alphas4F{a4F, 100, 0.9, 1001, 3};
    _AlphaQCD4F = [=] (double const& mu) -> double{ return Alphas4F.Evaluate(mu); };
    
    // Construct the DGLAP objects
    const auto Dglap = BuildDglap(_DglapObj,
      [=] (double const& x, double const&)->std::map<int,double>{
        return apfel::PhysToQCDEv(_inPDFs->xfxMap(x));
      },
    *Q0, PtOrder, _AlphaQCD);

    const auto Dglap3F = BuildDglap(_DglapObj3F,
      [=] (double const& x, double const&)->std::map<int,double>{
        return apfel::PhysToQCDEv(_inPDFs->xfxMap(x));
      },
    *Q0, PtOrder, _AlphaQCD3F);

    const auto Dglap4F = BuildDglap(_DglapObj4F,
      [=] (double const& x, double const&)->std::map<int,double>{
        return apfel::PhysToQCDEv(_inPDFs->xfxMap(x));
      },
    *Q0, PtOrder, _AlphaQCD4F);
    
    // Tabulate PDFs (ideally the parameters of the tabulation should
    // be read from parameters.yaml).
    _TabulatedPDFs = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>
      (new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*Dglap,
          QGrid[0].as<int>(),
          QGrid[1].as<double>(),
          QGrid[2].as<double>(),
          QGrid[3].as<int>()});

    _TabulatedPDFs3F = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>
      (new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*Dglap3F,
          QGrid[0].as<int>(),
          QGrid[1].as<double>(),
          QGrid[2].as<double>(),
          QGrid[3].as<int>()});

    _TabulatedPDFs4F = std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>
      (new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*Dglap4F,
          QGrid[0].as<int>(),
          QGrid[1].as<double>(),
          QGrid[2].as<double>(),
          QGrid[3].as<int>()});
    
  }
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> EvolutionAPFELxx::GetTabulatedPDFs() {return *_TabulatedPDFs;};
  std::function<double(double const& Q)>                 EvolutionAPFELxx::GetAlphaQCD() {return _AlphaQCD;};
  
  std::map<int,double>EvolutionAPFELxx::xfxQmap(double x,double Q){
    return apfel::QCDEvToPhys(_TabulatedPDFs->EvaluateMapxQ(x,Q));
  }
  double EvolutionAPFELxx::xfxQ(int i,double x,double Q){
    return _TabulatedPDFs->EvaluatexQ(i,x,Q);
  }
  void EvolutionAPFELxx::xfxQarray(double x,double Q,double*pdfs){
    // Get map of PDFs
    const std::map<int,double> fset = apfel::QCDEvToPhys(_TabulatedPDFs->EvaluateMapxQ(x, Q));
    pdfs[0] =fset.at(-6);
    pdfs[1] =fset.at(-5);
    pdfs[2] =fset.at(-4);
    pdfs[3] =fset.at(-3);
    pdfs[4] =fset.at(-2);
    pdfs[5] =fset.at(-1);
    pdfs[6] =fset.at(0);
    pdfs[7] =fset.at(1);
    pdfs[8] =fset.at(2);
    pdfs[9] =fset.at(3);
    pdfs[10]=fset.at(4);
    pdfs[11]=fset.at(5);
    pdfs[12]=fset.at(6);
  }
  double EvolutionAPFELxx::getAlphaS(double Q){
    return _AlphaQCD(Q);
  }

  std::map<int,double>EvolutionAPFELxx::xfxQmapNF(int nf,double x,double Q){
    switch(nf){
    case 3: return apfel::QCDEvToPhys(_TabulatedPDFs3F->EvaluateMapxQ(x,Q));
    case 4: return apfel::QCDEvToPhys(_TabulatedPDFs4F->EvaluateMapxQ(x,Q));
    case 5: return xfxQmap(x,Q);
    default:
      hf_errlog(26101301,"F: EvolutionAPFELxx::xfxQmapNF: nf must be 3, 4 or 5");
      return {};
    }
  }
  double EvolutionAPFELxx::getAlphaSNF(int nf,double Q){
    switch(nf){
    case 3: return _AlphaQCD3F(Q);
    case 4: return _AlphaQCD4F(Q);
    case 5: return _AlphaQCD(Q);
    default:
      hf_errlog(26101302,"F: EvolutionAPFELxx::getAlphaSNF: nf must be 3, 4 or 5");
      return 0;
    }
  }

  vector<double> EvolutionAPFELxx::getXgrid() {
    return _Grid->GetJointGrid().GetGrid();
  }

  vector<double> EvolutionAPFELxx::getQgrid() {
    return _TabulatedPDFs->GetQGrid();
  }

  //_________________________________________________________________________________
  std::map<int, std::map<int, apfel::Set<apfel::Operator>>> SchemeChangeKernels(apfel::Grid const& g, std::vector<double> const& Thresholds)
  {
    const YAML::Node node = XFITTER_PARS::rootNode["FactorisationScheme"];
    const std::string scheme = (node ? node.as<std::string>() : "MSbar");
    if (scheme == "MSbar")
      return {};
    if (OrderMap(XFITTER_PARS::getParamS("Order")) - 1 > 1)
      hf_errlog(26102401, "F: FactorisationScheme " + scheme + " is only available up to NLO");
    if (scheme == "Krk")       return apfel::InitializeSchemeChangeKernelsKrk(g, Thresholds);
    if (scheme == "PHYS")      return apfel::InitializeSchemeChangeKernelsPHYS(g, Thresholds);
    if (scheme == "DIS")       return apfel::InitializeSchemeChangeKernelsDIS(g, Thresholds);
    if (scheme == "POS")       return apfel::InitializeSchemeChangeKernelsPOS(g, Thresholds);
    if (scheme == "MPOS")      return apfel::InitializeSchemeChangeKernelsMPOS(g, Thresholds);
    if (scheme == "MPOSDelta") return apfel::InitializeSchemeChangeKernelsMPOSDelta(g, Thresholds);
    if (scheme == "Aversa")    return apfel::InitializeSchemeChangeKernelsAversa(g, Thresholds);
    hf_errlog(26102402, "F: Unknown FactorisationScheme " + scheme);
    return {};
  }
}

