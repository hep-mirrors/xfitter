#include <IntegrateDIS.h>
#include <xfitter_cpp_base.h>
#include <cassert>

int IntegrateDIS::init(const double s,
                        std::valarray<double>* q2minp, std::valarray<double>* q2maxp,
                        std::valarray<double>* yminp, std::valarray<double>* ymaxp,
                        std::valarray<double>* xminp, std::valarray<double>* xmaxp)
{
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

  for(size_t i = 0; i < q2minp->size(); i++)
  {
    const double q2min = (*q2minp)[i];
    const double q2max = (*q2maxp)[i];
    const double ymin = (*yminp)[i];
    const double ymax = (*ymaxp)[i];

    //if(ymin == 0.)

    // just copy previous entry, if binning is the same
    //bool flagCopyValue = false;
    if(i > 0 && q2min == (*q2minp)[i - 1] && q2max == (*q2maxp)[i - 1] && ymin == (*yminp)[i - 1] && ymax == (*ymaxp)[i - 1])
      //flagCopyValue = true;
      //if(flagCopyValue)
    {
      // -1 means that the value from previous bin will be copied in compute()
      _nSubBins[i] = -1;
      continue;
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
        }
      }
    }
    assert(nSubBins == j);
    _nSubBins[i] = nSubBins;
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
        xsec += dxsec;
        printf("%f %f: %f += %f [ %f * %f * %f ]\n", _q2[i], _x[i],
               xsec, dxsec, val[i], _deltaq2[i], _deltax[i]);
      }
      valIntegrated[nBins] = xsec;
      // increase number of processed bins
      nSubBins += n;
    }
    nBins++;
  }
  assert(valIntegrated.size() == nBins);
  return valIntegrated;
}
