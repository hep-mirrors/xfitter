#include "ReactionTheory.h"
#include "ReactionBaseDISNC.h"
#include <IntegrateDIS.h>
#include <hoppet_v1.h>

/**
  @class' ReactionHOPPET_DISNC

  @brief HOPPET DISNC reaction

  @version 0.1
  @date 2024-08-13
  */

class ReactionHOPPET_DISNC : public ReactionBaseDISNC
{
private:
    typedef ReactionBaseDISNC Super;

public:
    ReactionHOPPET_DISNC(){};

public:
    virtual string getReactionName() const override { return "HOPPET_DISNC"; };
    void virtual atStart() override final;
    virtual void initTerm(TermData *td) override final;
    void atIteration();

protected:
	  void F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal);
	  void FL(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal);
	  void xF3(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal);

private:
    void calcF2FLF3(unsigned dataSetID);

    map<int, valarray<double>> _f2;
    map<int, valarray<double>> _fl;
    map<int, valarray<double>> _f3;

    int _order;
    double _xmuR;
    double _xmuF;
    double _muR_Q;
    double _dy;
    int _param_coefs;
    double* _convfac;
    double* _alphaem;
    double* _Mz;
    double* _Mw;
    double* _sin2thetaW;
    // temporary, needed for alphaS evolution
    double* _alphas;
    double _Q0;
};
