#ifndef __INTEGRATEDIS_H
#define __INTEGRATEDIS_H
#include <valarray>
#include <string>
class TermData;

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

    int _init(const std::valarray<double>* sp,
             const std::valarray<double>* q2minp, const std::valarray<double>* q2maxp,
             const std::valarray<double>* yminp, const std::valarray<double>* ymaxp,
             const std::valarray<double>* xminp, const std::valarray<double>* xmaxp);

  public:
    //IntegrateDIS();

    // try to initialise integrated cross section for one dataset, return true on success
    bool init_from_td(TermData* td, bool isReduced, std::string& msg);

    // calculate integrated cross sections by integrating over subbins
    std::valarray<double> compute(const std::valarray<double>& val);

    // get bin values
    std::valarray<double>* getBinValuesQ2() { return &_q2; }
    std::valarray<double>* getBinValuesX() { return &_x; }
    std::valarray<double>* getBinValuesY() { return &_y; }

    // get number of bins
    size_t getNPoints() {return _q2.size();};
};
#endif
