/*
   @file Reaction_DISNC_Hoppet.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/

#include "Reaction_DISNC_Hoppet.h"
#include <iostream>
#include <cstdio>
#include "ReactionBaseDISNC.h"
#include "hoppet_v1.h" // Include the HOPPET header
//#include <xfitter_cpp.h>
#include "xfitter_pars.h"
//#include "xfitter_cpp_base.h"
#include <BaseEvolution.h>
#include <hf_errlog.h>



using namespace hoppetv1;
// the class factories
extern "C" Reaction_DISNC_Hoppet *create()
{
  return new Reaction_DISNC_Hoppet();
}

// Initialize at the start of the computation
//void Reaction_DISNC_Hoppet::atStart()
//{
//}

// Main function to compute results at an iteration
void Reaction_DISNC_Hoppet::F2(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
	const double xmuR = 1;
    const double xmuF = 1;
	const double mc = 1.414213563;
    const double mb = 4.5;
    const double mt = 175.0;
    
    auto &x = *GetBinValues(td, "x");
    auto &Q2 = *GetBinValues(td, "Q2");
    auto &y = *GetBinValues(td, "y");
    const double pi = 3.1415926535897932384626433832795029;
    valarray<double> yplus = 1.0 + (1.0 - y) * (1.0 - y);
    valarray<double> factor = 2 * pi * _alphaem * _alphaem * yplus / (Q2 * Q2 * x) * _convfac;
  //  val *= factor;
    const bool param_coefs = true;
    int	 order_max = 4;

    hoppetInitStrFct(order_max, param_coefs, xmuR, xmuF);

  
    //std::vector<double> StrFct(6); // Adjust size if needed
    
   const auto _convfac = *XFITTER_PARS::getParamD("convFac");
   const auto _alphaem = *XFITTER_PARS::getParamD("alphaem");
   const auto MZ = *XFITTER_PARS::getParamD("Mz");
   const auto MW = *XFITTER_PARS::getParamD("Mw");
   const auto _sin2thetaW = *XFITTER_PARS::getParamD("sin2thW");
  
  const double MW2 = MW*MW; 
  const double MZ2 = MZ*MZ;

    // Construct structure functions
  const double s2tw = 1 - MW2 / MZ2;
  const double VD   = - 0.5 + 2 * s2tw / 3;
  const double VU   = + 0.5 - 4 * s2tw / 3;
  const double AD   = - 0.5;
  const double AU   = + 0.5;
  const double Ve   = - 0.5 + 2 * s2tw;
  const double Ae   = - 0.5;
  
  std::vector<double> StrFct(12);
 
   for (int y = 0; y <= x.size(); y++)
   {    const double Q = std::sqrt(Q2[y]);
	    const double PZ  = Q / ( Q + MZ2 ) / ( 4 * s2tw * ( 1 - s2tw ) );
        const double PZ2 = PZ * PZ;
	   // const double x = 1 / ( 1 + exp(y) );
  hoppetStrFct(x[y],Q, xmuR * Q, xmuF * Q, &StrFct[0]);

  //  const double F1NCh = StrFct[iF1EM] + StrFct[iF1Z] * (Ve * Ve + Ae * Ae) * PZ2 - StrFct[iF1gZ] * Ve * PZ;
    const double F2NCh = StrFct[iF2EM] + StrFct[iF2Z] * (Ve * Ve + Ae * Ae) * PZ2 - StrFct[iF2gZ] * Ve * PZ;
   // const double F3NCh = 2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ;
    // const double FLNCh = F2NCh - 2 * x * F1NCh;
    
    valExternal[y] = F2NCh;
    
}

}



void Reaction_DISNC_Hoppet::atIteration()
{
  // Make sure to call the parent class initialization:
  //super::atIteration();
   super::atIteration();

  

  _ve = -0.5 + 2. * _sin2thetaW; // !
  _ae = -0.5;                    // !
  _au = 0.5;
  _ad = -0.5;
  _vu = _au - (4. / 3.) * _sin2thetaW;
  _vd = _ad + (2. / 3.) * _sin2thetaW;

  //  print (_Mz);

  // Re-set internal maps (faster access):
  for (auto ds : _dsIDs)
  {
    (_f2u[ds])[0] = -100.;
    (_flu[ds])[0] = -100.;
    (_xf3u[ds])[0] = -100.;
  }

}

//
void Reaction_DISNC_Hoppet::initTerm(TermData *td)
{
  unsigned termID = td->id;

  {
    const string &name = getReactionName();
    if (name == "BaseDISNC" or name == "RT_DISNC")
    {
      xfitter::BaseEvolution *pdf = td->getPDF();
      if (pdf->getClassName() != string("HOPPET"))
      {
        std::cerr << "[ERROR] " << getReactionName() << " can only work with HOPPET evolution; got evolution \"" << pdf->_name << "\" of class \"" << pdf->getClassName() << "\" for termID=" << termID << std::endl;
        hf_errlog(19052311, "F: Chosen DISNC reaction can only work with HOPPET evolution, see stderr for details");
      }

      //xfitter::requireZMSTF();
    }
  }

  _dsIDs.push_back(termID);

  _tdDS[termID] = td;
  _polarisation[termID] = (td->hasParam("epolarity")) ? *td->getParamD("epolarity") : 0;
  _charge[termID] = (td->hasParam("echarge")) ? *td->getParamD("echarge") : 0;

  _dataType[termID] = dataType::sigred;
  _dataFlav[termID] = dataFlav::incl;
  string msg = "I: Calculating DIS NC reduced cross section";

  if (td->hasParam("type"))
  {
    string type = td->getParamS("type");
    if (type == "signonred")
    {
      _dataType[termID] = dataType::signonred;
      msg = "I: Calculating DIS NC double-differential (non-reduced) cross section";
    }
    else if (type == "sigred")
    {
      _dataType[termID] = dataType::sigred;
      msg = "I: Calculating DIS NC reduced cross section";
    }
   /* else if (type == "sigred_noF3")
    {
      _dataType[termID] = dataType::sigred_nof3;
      msg = "I: Calculating DIS NC reduced cross section w/o F3";
    }*/
    else if (type == "F2")
    {
      _dataType[termID] = dataType::f2;
      msg = "I: Calculating DIS NC F2";
    }
    else if (type == "FL")
    {
      _dataType[termID] = dataType::fl;
      msg = "I: Calculating DIS NC FL";
    }
   /* else if (type == "F3")
    {
      _dataType[termID] = dataType::f3;
      msg = "I: Calculating DIS NC F3";
    }*/
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown type = %s", termID, type.c_str());
      string str = buffer;
      hf_errlog(17101901, str);
    }
  }

  if (td->hasParam("flav"))
  {
    string flav = td->getParamS("flav");
    if (flav == "incl")
    {
      _dataFlav[termID] = dataFlav::incl;
      msg += " inclusive";
    }
    else if (flav == "c")
    {
      _dataFlav[termID] = dataFlav::c;
      msg += " charm";
    }
    else if (flav == "b")
    {
      _dataFlav[termID] = dataFlav::b;
      msg += " beauty";
    }
    else
    {
      char buffer[256];
      sprintf(buffer, "F: dataset with id = %d has unknown flav = %s", termID, flav.c_str());
      string str = buffer;
      hf_errlog(17101902, str);
    }
  }

  double s = -1.0;
  if (td->hasParam("energy"))
    s = pow(*td->getParamD("energy"), 2.0);

  auto *q2minp = td->getBinColumnOrNull("Q2min");
  auto *q2maxp = td->getBinColumnOrNull("Q2max");
  if (!q2minp)
    q2minp = td->getBinColumnOrNull("q2min");
  if (!q2maxp)
    q2maxp = td->getBinColumnOrNull("q2max");
  auto *yminp = td->getBinColumnOrNull("ymin");
  auto *ymaxp = td->getBinColumnOrNull("ymax");
  auto *xminp = td->getBinColumnOrNull("xmin");
  auto *xmaxp = td->getBinColumnOrNull("xmax");

  /*if (q2minp && q2maxp)
  {
    if (s < 0)
      hf_errlog(18060100, "F: centre-of-mass energy is required for integrated DIS dataset " + std::to_string(termID));
    if (_dataType[termID] == dataType::f3)
      hf_errlog(18060101, "F: F3 is not supported for integrated DIS datasets " + std::to_string(termID));

    _integrated[termID] = new IntegrateDIS(*q2minp, *q2maxp, *yminp, *ymaxp, *xminp, *xmaxp, s);
  }*/

 // hf_log(msg);
};

void Reaction_DISNC_Hoppet::F2gamma BASE_PARS
{
  valarray<double> f2u, f2d;
}

void Reaction_DISNC_Hoppet::F2gammaZ BASE_PARS
{
  valarray<double> f2u, f2d;
}

void Reaction_DISNC_Hoppet::F2Z BASE_PARS
{
  valarray<double> f2u, f2d;
}

void Reaction_DISNC_Hoppet::FLgamma BASE_PARS
{
  valarray<double> flu, fld;
}

void Reaction_DISNC_Hoppet::FLgammaZ BASE_PARS
{
  valarray<double> flu, fld;
}

void Reaction_DISNC_Hoppet::FLZ BASE_PARS
{
  valarray<double> flu, fld;
}

void Reaction_DISNC_Hoppet::xF3gammaZ BASE_PARS
{
  valarray<double> xf3u, xf3d;
}

void Reaction_DISNC_Hoppet::xF3Z BASE_PARS
{
  valarray<double> xf3u, xf3d;
}

     

/*
int main()
{
	  const int nflav     = -5;
      const int sc_choice = 1 ;
      const double zmass  = apfel::ZMass;
      const double wmass  = apfel::WMass;
      const double s2tw   = 1 - pow(wmass / zmass, 2);
      const double Ve     = - 0.5 + 2 * s2tw;
      const double Ae     = - 0.5;
	
	  const double xbmin = 0.00001;
      const double xbmax = 0.95;
      const double ybmin = log(( 1 - xbmax ) / xbmax);
      const double ybmax = log(( 1 - xbmin ) / xbmin);
      const double ybstp = ( ybmax - ybmin ) / nyb;
      const std::vector<double> muv{2, 5, 10, 50, 100};
      const double pdf[13];
      const double StrFct[14];
      
      
	for (double Q : muv)
	{
	  const double PZ  = pow(Q, 2) / ( pow(Q, 2) + pow(zmass, 2) ) / ( 4 * s2tw * ( 1 - s2tw ) );
	  const double PZ2 = PZ * PZ;
	
	
	 for (double y = ybmin; y <= 1.0000001 * ybmax; y += ybstp)
	    {
	      const double x = 1 / ( 1 + exp(y) );
	      hoppetStrFct(x, Q, xmuR * Q, xmuF * Q, StrFct);

	      // Construct structure functions
	      const double F1NCh = StrFct[iF1EM] + StrFct[iF1Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF1gZ] * Ve * PZ;
	      const double F2NCh = StrFct[iF2EM] + StrFct[iF2Z] * ( Ve * Ve + Ae * Ae ) * PZ2 - StrFct[iF2gZ] * Ve * PZ;
	      const double F3NCh = 2 * StrFct[iF3Z] * Ae * Ve * PZ2 - StrFct[iF3gZ] * Ae * PZ;
	      const double FLNCh = F2NCh - 2 * x * F1NCh;

	     }
   }
	}
*/
