#ifndef __INTEGRATEDIS_H
#define __INTEGRATEDIS_H
#include <valarray>

// Class used by ReactionBaseDISNC and ReactionBaseDISCC for integrating DIS cross sections
// and providing them over Q2, y, x ranges.
// Method is based on legacy subroutine GetIntegratedDisXsection from dis_sigma.f

class IntegrateDIS
{
  private:
    // number of Q2, x subbins are hardcoded now
    // in the future they could be specified optionally in data files and passed here from DIS reactions classes
    const int _nsplit_x = 25;
    const int _nsplit_q2 = 25;
    //const int _nsplit_x = 100;
    //const int _nsplit_q2 = 100;
    //const int _nsplit_x = 5;
    //const int _nsplit_q2 = 5;
    std::valarray<int> _nSubBins;
    std::valarray<double> _q2;
    std::valarray<double> _x;
    std::valarray<double> _y;
    std::valarray<double> _deltaq2;
    std::valarray<double> _deltax;

  public:
    //IntegrateDIS();

    // initialise integrated cross section for one dataset, return number of subbins
    int init(const double s,
             const std::valarray<double>* q2minp, const std::valarray<double>* q2maxp,
             const std::valarray<double>* yminp, const std::valarray<double>* ymaxp,
             const std::valarray<double>* xminp, const std::valarray<double>* xmaxp);

    // calculate integrated cross sections by integrating over subbins
    std::valarray<double> compute(const std::valarray<double>& val);

    // get bin values
    std::valarray<double>* getBinValuesQ2() { return &_q2; }
    std::valarray<double>* getBinValuesX() { return &_x; }
    std::valarray<double>* getBinValuesY() { return &_y; }
};
#endif
