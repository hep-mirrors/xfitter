#include "ReactionTheory.h"
#include "ReactionBaseDISCC.h"
#include <IntegrateDIS.h>
#include <hoppet_v1.h>

/**
  @class' ReactionHOPPET_DISCC

  @brief HOPPET DISCC reaction

  @version 0.1
  @date 2024-08-13
  */

class ReactionHOPPET_DISCC : public ReactionBaseDISCC
{
private:
    typedef ReactionBaseDISCC Super;

public:
    ReactionHOPPET_DISCC(){};

public:
    virtual string getReactionName() const override { return "HOPPET_DISCC"; };
    void virtual atStart() override final;
    virtual void initTerm(TermData *td) override final;
    void atIteration();

protected:
    virtual valarray<double> F2(TermData *td) override final;
    virtual valarray<double> FL(TermData *td) override final;
    virtual valarray<double> xF3(TermData *td) override final;

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

    // temporary: allow different orders in evolution and DIS SFs
    int _order_HOPPET_Evolution;
};
