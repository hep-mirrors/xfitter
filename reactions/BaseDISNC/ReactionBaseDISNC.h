#pragma once

#include "ReactionTheory.h"
#include <IntegrateDIS.h>

/**
  @class' ReactionBaseDISNC

  @brief A wrapper class for BaseDISNC reaction

  Based on the ReactionTheory class.

  @version 0.1
  @date 2017-04-08
  */

// Define standard parameters used by SF and x-sections:
#define BASE_PARS (TermData * td, valarray<double> & val, map<string, valarray<double>> & err)

class ReactionBaseDISNC : public ReactionTheory
{
public:
   ReactionBaseDISNC(){};

public:
   virtual string getReactionName() const override { return "BaseDISNC"; };
   virtual void atStart() override;
   virtual void initTerm(TermData *td) override;
   virtual void reinitTerm(TermData *td) override; //! allow for polarisation update.

   //!< Initialize all EWK couplings here:
   virtual void atIteration() override;
   virtual void compute(TermData *, valarray<double> &val, map<string, valarray<double>> &err) override;

protected:
   enum class dataType
   {
      signonred,
      sigred,
      f2,
      fl,
      f3,
      sigred_nof3
   }; //!< Define compute output.
   enum class dataFlav
   {
      incl,
      c,
      b
   }; //!< Define final state.

   /*
       A few methods specific for DIS NC process.
    */

   virtual void F2gamma BASE_PARS;
   virtual void F2gammaZ BASE_PARS;
   virtual void F2Z BASE_PARS;

   //! compute full F2
   virtual void F2 BASE_PARS;

   virtual void FLgamma BASE_PARS;
   virtual void FLgammaZ BASE_PARS;
   virtual void FLZ BASE_PARS;

   //!< compute full FL
   virtual void FL BASE_PARS;

   virtual void xF3gammaZ BASE_PARS;
   virtual void xF3Z BASE_PARS;

   //!< compute full xF3
   virtual void xF3 BASE_PARS;

   //! reduced cross section
   virtual void sred BASE_PARS;

   // Helper functions:
   void kappa(TermData *td, valarray<double> &k);

private:
   map<unsigned, int> _npoints;         //!< Number of points in a dataset.
   map<unsigned, double> _polarisation; //!< longitudinal polarisation
   map<unsigned, double> _charge;       //!< lepton beam charge
   map<unsigned, dataType> _dataType;   //!< cross section (reduced, F2, FL)
   map<unsigned, dataFlav> _dataFlav;   //!< flavour (incl, c, b)

protected:
   // some parameters which may change from iteration to iteration:
   vector<unsigned> _dsIDs; //!< list of termIDs managed by the reaction.
   map<unsigned, TermData*> _tdDS; //! store term data for later access. 
   double _alphaem;
   double _Mz;
   double _Mw;
   double _sin2thetaW; 
   double _ae, _ve;
   double _au, _ad;
   double _vu, _vd;

   // conversion constant factor
   double _convfac;

protected:
   const int GetNpoint(unsigned termID) { return _npoints[termID]; }
   const double GetPolarisation(unsigned termID) { return _polarisation[termID]; }
   const double GetCharge(unsigned termID) { return _charge[termID]; }
   const dataType GetDataType(unsigned termID) { return _dataType[termID]; }
   const dataFlav GetDataFlav(unsigned termID) { return _dataFlav[termID]; }
   TermData* GetTermData(unsigned termID) { return _tdDS[termID];}

   // Another decomposition:
   virtual void GetF2ud(TermData *td, valarray<double> &f2u, valarray<double> &f2d);
   virtual void GetFLud(TermData *td, valarray<double> &flu, valarray<double> &fld);
   virtual void GetxF3ud(TermData *td, valarray<double> &xf3u, valarray<double> &xf3d);

private:
   // Some buffering mechanism to avoid double calls
   map<unsigned, valarray<double>> _f2u; //!< F2 for u-type quarks
   map<unsigned, valarray<double>> _f2d; //!< F2 for d-type quarks
   map<unsigned, valarray<double>> _flu; //!< FL for u-type quarks
   map<unsigned, valarray<double>> _fld; //!< FL for d-type quarks
   map<unsigned, valarray<double>> _xf3u;
   map<unsigned, valarray<double>> _xf3d;
   // Store QCDNUM PDF ID
   map<unsigned,int> _ipdfSet;
protected:
   // for integrated cross sections
   // method is based on legacy subroutine GetIntegratedDisXsection
   map<unsigned, IntegrateDIS *> _integrated;
   virtual const valarray<double> *GetBinValues(TermData *td, const string &binName); //! interface for integerated sigma
   // higher twist
   map<unsigned, bool> _flag_ht;
   std::vector<double> _ht_x;
   std::vector<double> _ht_2;
   std::vector<double> _ht_t;
   double _ht_alpha_2;
   double _ht_alpha_t;
   void ApplyHigherTwist(TermData *td, const int f_type, valarray<double>& val, map<string, valarray<double>>& err);
};
