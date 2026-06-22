#pragma once
class TermData;

// DIS nuclear corrections according to 
// S.A. Kulagin and R. Petti, Nucl.Phys.A 765 (2006) 126-187 [hep-ph/0412425]
class DIS_NUKE {
  // TODO check if this can be used at the level of BaseDIS reaction classes
  friend class ReactionBaseFFABM;
  friend class ReactionFFABM_DISNC;
  friend class ReactionFFABM_DISCC;
  public:
    DIS_NUKE(TermData* td);
  private:
    int _ftyp = 0; // nuclear correction: 1 Carbon (A=12, Z=6), 2 Iron (A=56, Z=26), 3 Lead (A=207, Z=82), 4 Emulsion (mixture of nuclei) 
    int _kint = 0; // nuclear interaction: 1 charged lepton, 2 neutrino (-2 antineutrino) CC, 3 neutrino (-3 antineutrino) NC, 4 neutrino CC charm production (-4 antineutrino), 5 neutrino CC not-charm (-5 antineutrino) [4, 5 not implemented]
    int _kord = 0; // order: 1 LO, 2 NLO, 3 NNLO [not implemented]
    bool _need_f3bar = false;
    // private methods can be used by friends
    bool need_f3bar() const {return _need_f3bar;};
    double apply(const double q2, const double x, double& f2, double& fl, double& f3, double const* f3_bar_ptr=nullptr) const;
};
