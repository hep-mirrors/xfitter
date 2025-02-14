#include <IntegrateDIS.h>
#include <xfitter_cpp_base.h>
#include <cassert>
#include <iostream>
#include <TermData.h>
using std::cerr;
using std::endl;

//IntegrateDIS* IntegrateDIS::init_from_td(TermData* td, string& msg) {
bool IntegrateDIS::init_from_td(TermData* td, bool isReduced, std::string& msg) {
  // bins
  // if Q2min, Q2max, ymin and ymax (and optionally xmin, xmax) are provided, integrated cross sections are calculated
  auto *q2minp = td->getBinColumnOrNull("Q2min");
  auto *q2maxp = td->getBinColumnOrNull("Q2max");
  // also try small first letter for Q2 (for backward compatibility)
  if (!q2minp)
    q2minp = td->getBinColumnOrNull("q2min");
  if (!q2maxp)
    q2maxp = td->getBinColumnOrNull("q2max");
  auto *yminp = td->getBinColumnOrNull("ymin");
  auto *ymaxp = td->getBinColumnOrNull("ymax");
  // optional xmin, xmax for integrated cross sections
  auto *xminp = td->getBinColumnOrNull("xmin");
  auto *xmaxp = td->getBinColumnOrNull("xmax");
  // check if centre-of-mass energy is provided
  //double s = -1.0;
  auto *sp = td->getBinColumnOrNull("s");
  bool sp_delete = 0;
  if (!sp) {
    if (td->hasParam("energy"))
    {
      double energy = *td->getParamD("energy");
      sp = new valarray<double>(energy * energy, td->getNbins());
      sp_delete = 1;
    }
    else if (td->hasParam("mh") && td->hasBinColumn("beam_energy")) {
      sp = new valarray<double>(2.*(*td->getParamD("mh"))*(*td->getBinColumnOrNull("beam_energy")));
      sp_delete = 1;
    }
  }

  if (q2minp)
  {
    // integrated cross section
    if (!sp)
      hf_errlog(18060100, "F: centre-of-mass energy is required for integrated DIS term " + std::to_string(td->id));
    if (isReduced)
      hf_errlog(18060200, "F: integrated DIS can be calculated only for non-reduced cross sections, term " + std::to_string(td->id));
    if(!q2maxp) {
      q2maxp = new valarray<double>((*sp)); // TODO fix memory leak
    }
    this->_init(sp, q2minp, q2maxp, yminp, ymaxp, xminp, xmaxp);
    if(sp_delete) {
      delete sp;
    }
    msg += " (integrated)";
    return true;
  }
  return false;
}

int IntegrateDIS::_init(const std::valarray<double>* sp,
                        const std::valarray<double>* q2minp, const std::valarray<double>* q2maxp,
                        const std::valarray<double>* yminp, const std::valarray<double>* ymaxp,
                        const std::valarray<double>* xminp, const std::valarray<double>* xmaxp)
{
  for (size_t i = 0; i < sp->size(); i++) {printf("s[%ld] = %f\n", i, (*sp)[i]);}
  // prepare number of bins and arrays
  const int npoints = q2minp->size();
  _nSubBins.resize(npoints);
  int nSubBins = _nsplit_x * _nsplit_q2;
  _q2.resize(npoints * nSubBins);
  _deltaq2.resize(npoints * nSubBins);
  _x.resize(npoints * nSubBins);
  _deltax.resize(npoints * nSubBins);
  _y.resize(npoints * nSubBins);
  int currentSubBin = 0;

  printf("SZ q2minp->size() = %ld\n", q2minp->size());
  for(size_t i = 0; i < q2minp->size(); i++)
  {
    const double s = (*sp)[i];
    const double q2min = (*q2minp)[i];
    const double q2max = (*q2maxp)[i];
    // if ymin is not provided, assume there is no lower boundary on y
    double ymin = (yminp) ? (*yminp)[i] : (q2min / s);
    // if ymax is not provided, assume there is no upper boundary on y
    const double ymax = (ymaxp) ? (*ymaxp)[i] : 1.;

    // just copy previous entry, if binning is the same
    // Actually, what is the point of having two points with the same binning? Why do we have this check? --Ivan
    if(i > 0
      && q2min == (*q2minp)[i - 1]
      && q2max == (*q2maxp)[i - 1]
      && (!yminp || ymin == (*yminp)[i - 1])
      && (!ymaxp || ymax == (*ymaxp)[i - 1])
      //Check that xmin and xmax are the same as in previous bin
      //But only if xmin and xmax columns exist
      && (
        xminp==nullptr ||
        ( (*xminp)[i] == (*xminp)[i - 1] )
      )
      && (
        xmaxp==nullptr ||
        ( (*xmaxp)[i] == (*xmaxp)[i - 1] )
      )
      && s == (*sp)[i-1]
    )
    {
      // -1 means that the value from previous bin will be copied in compute()
      _nSubBins[i] = -1;
      hf_errlog(19091101, "I: IntegrateDIS: two datapoints have identical bin, will only calculate cross-section once for both");
      continue;
    }

    // avoid division by zero
    if (ymin == 0.) {
      cerr<<"[ERROR] In IntegrateDIS: ymin==0 would lead to division by zero, please edit your datafile."
      " Bin id="<<i<<"; q2min="<<q2min<<"; q2max="<<q2max<<";";
      if(xminp)cerr<<" xmin="<<(*xminp)[i]<<";";
      else cerr<<" no xmin given;";
      if(xmaxp)cerr<<" xmax="<<(*xmaxp)[i]<<";";
      else cerr<<" no xmax given;";
      cerr<<" s="<<s<<";"<<endl;
      hf_errlog(19091100, "F: ymin==0 in IntegrateDIS leads to division by zero, see stderr");
    }

    // determine x range
    double xmin = q2min / (s * ymax);
    double xmax = q2max / (s * ymin);
    if(xminp)
      xmin = std::max(xmin, (*xminp)[i]);
    if(xmaxp)
      xmax = std::min(xmax, (*xmaxp)[i]);
    if(xmax > 1.0)
      xmax = 0.999999;
      //hf_errlog(18060101, "xmaxCalc > 1.0");

    //printf("DIS integrated q2min = %f q2max = %f xmin = %f xmax = %f ymin = %f ymax = %f\n", q2min, q2max, xmin, xmax, ymin, ymax);
    // do integration in log space for Q2 and x
    int j = 0;
    double q2_1 = -1.0;
    double q2_2 = -1.0;
    double x_1 = -1.0;
    double x_2 = -1.0;
    for(int iq2 = 0; iq2 <= _nsplit_q2; iq2++)
    {
      q2_1 = q2_2;
      q2_2 = exp(log(q2min) + (log(q2max) - log(q2min)) / _nsplit_q2 * iq2);
      if(iq2 > 0)
      {
        for(int ix = 0; ix <= _nsplit_x; ix++)
        {
          x_1 = x_2;
          x_2 = exp(log(xmin) + (log(xmax) - log(xmin)) / _nsplit_x * ix);
          if(ix > 0)
          {
            j = j + 1;
            _deltaq2[currentSubBin] = q2_2 - q2_1;
            double q2 = exp(log(q2_1) + 0.5 * (log(q2_2) - log(q2_1)));
            _q2[currentSubBin] = q2;
            _deltax[currentSubBin] = x_2 - x_1;
            double x = exp(log(x_1) + 0.5 * (log(x_2) - log(x_1)));
            _x[currentSubBin] = x;
            double y = q2 / (s * x);
            // check that calculated y values agree with limits given in data file
            // TODO: this is not optimal, better would be to skip such points
            // however, then arrays of subbins will have different size
            if(y < ymin || y > ymax)
              y = 0.0;
            _y[currentSubBin] = y;
            currentSubBin++;
          }
          //printf("bin %ld  _q2 = %f  x = %f  y = %f\n", i, _q2[i], _x[i], _y[i]);
        }
      }
    }
    assert(nSubBins == j);
    _nSubBins[i] = nSubBins;
    printf("bin %ld  q2 = %f %f  x = %f %f  y = %f %f\n", i, q2min, q2max, xmin, xmax, ymin, ymax);
  }
  return currentSubBin;
}

std::valarray<double> IntegrateDIS::compute(const std::valarray<double>& val)
{
  std::valarray<double> valIntegrated(_nSubBins.size());
  unsigned int nSubBins = 0;
  unsigned int nBins = 0;
  for(auto& n : _nSubBins)
  {
    if(n == -1)
    {
      // copy cross section from previous bin
      valIntegrated[nBins] = valIntegrated[nBins - 1];
    }
    else
    {
      // integrate
      double xsec = 0.0;
      for(size_t i = nSubBins; i < nSubBins + n; i++)
      {
        // skip points out of kinematic range
        if(_y[i] == 0.0)
          continue;
        double dxsec = val[i] * _deltaq2[i] * _deltax[i];
        printf("SZ Q2,x = %e,%e -> xsec = %e\n", _q2[i], _x[i], val[i]);
        xsec += dxsec;
        //printf("%f %f %f: %f += %f [ %f * %f * %f ]\n", _q2[i], _x[i], _y[i],
        //       xsec, dxsec, val[i], _deltaq2[i], _deltax[i]);
      }
      valIntegrated[nBins] = xsec;
      // increase number of processed bins
      nSubBins += n;
    }
    nBins++;
  }
  fflush(stdout);
  //throw 24;
  assert(valIntegrated.size() == nBins);
  return valIntegrated;
}
