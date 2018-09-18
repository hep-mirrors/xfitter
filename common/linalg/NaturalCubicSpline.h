#ifndef NATURALCUBICSPLINE_H
#define NATURALCUBICSPLINE_H

// code mostly from https://stackoverflow.com/questions/1204553/are-there-any-good-libraries-for-solving-cubic-splines-in-c

#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
using namespace std;

typedef vector<double> vec;

class NaturalCubicSpline
{
  public:
    NaturalCubicSpline(vec &x, vec &y)
    {
      _vectorSplineSet = Spline(x, y);
    }

    double Eval(double x, const bool flagDerivative = false)
    {
      for(size_t i = 1; i < _vectorSplineSet.size(); i++)
      {
        if(_vectorSplineSet[i].x > x || i == (_vectorSplineSet.size() - 1))
        {
          //printf("Eval: x = %f in section %lu [%f %f]\n", x, i, _vectorSplineSet[i - 1].x, _vectorSplineSet[i].x);
          if( (i == 1 && x < _vectorSplineSet[0].x) || (i == (_vectorSplineSet.size() - 1) && x > _vectorSplineSet[_vectorSplineSet.size() - 1].x) )
            printf("Warning: x = %f outside spline range [%f %f]\n", x, _vectorSplineSet[0].x, _vectorSplineSet[_vectorSplineSet.size() - 1].x);
          const SplineSet& s = _vectorSplineSet[i - 1];
          double y = 0.0;
          double dx = x - s.x;
          if(flagDerivative)
            y = 3 * s.d * dx * dx + 2 * s.c * dx + s.b;
          else
            y = s.d * dx * dx * dx + s.c * dx * dx + s.b * dx + s.a;
          //printf("y = %f\n", y);
          return y;
        }
      }
      // should not be here
      throw;
    }

  private:
    struct SplineSet
    {
        double a;
        double b;
        double c;
        double d;
        double x;
    };
    std::vector<SplineSet> _vectorSplineSet;

    std::vector<SplineSet> Spline(vec &x, vec &y)
    {
      // check input
      // at least 4 points
      if(x.size() < 4)
        hf_errlog(18091000, "F: Natural cubic spline needs at least 4 input points");
      // x in ascending order
      for(size_t s = 1; s < x.size(); s++)
        if(x[s] <= x[s - 1])
          hf_errlog(18091001, "F: Natural cubic spline needs x points in accessing order");
      // x and y have same size
      if(x.size() != y.size())
        hf_errlog(18091002, "F: Natural cubic spline needs same number of x and y points");

      int n = x.size()-1;
      vec a;
      a.insert(a.begin(), y.begin(), y.end());
      vec b(n);
      vec d(n);
      vec h;

      for(int i = 0; i < n; ++i)
        h.push_back(x[i+1]-x[i]);

      vec alpha;
      for(int i = 0; i < n; ++i)
        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

      vec c(n+1);
      vec l(n+1);
      vec mu(n+1);
      vec z(n+1);
      l[0] = 1;
      mu[0] = 0;
      z[0] = 0;

      for(int i = 1; i < n; ++i)
      {
        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
      }

      l[n] = 1;
      z[n] = 0;
      c[n] = 0;

      for(int j = n-1; j >= 0; --j)
      {
        c[j] = z [j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
        d[j] = (c[j+1]-c[j])/3/h[j];
      }

      vector<SplineSet> output_set(n);
      for(int i = 0; i < n; ++i)
      {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
      }
      return output_set;
    }
};

/*int main()
{
    vec x(11);
    vec y(11);
    for(int i = 0; i < x.size(); ++i)
    {
        x[i] = i;
        y[i] = sin(i);
    }

    vector<SplineSet> cs = spline(x, y);
    for(int i = 0; i < cs.size(); ++i)
        cout << cs[i].d << "\t" << cs[i].c << "\t" << cs[i].b << "\t" << cs[i].a << endl;
}*/

#endif // NATURALCUBICSPLINE_H
