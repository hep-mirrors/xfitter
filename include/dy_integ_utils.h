#ifndef UTILS_NAMESPACE
#define UTILS_NAMESPACE 1


namespace Utils
{
// calculates numerical integral using simpson rule
int simpson( double (*F)(double, void*), void *, const double , 
            const double , const int , double &, double &);

// simpson with binning adjusted for W-mass integration
int simpsonW( double (*F)(double, void*), void *, const double , 
            const double , const int , double &, double &);

// simpson with binning adjusted for Z-mass integration
int simpsonZ( double (*F)(double, void*), void *, const double , 
            const double , const int , double &, double &);

// transforms cosine theta to a frame defined by b and g
// arguments are betta, gamma and cosine theta
double costh_LT(const double &, const double &, const double &);

}

#endif
