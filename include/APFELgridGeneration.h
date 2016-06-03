//
// APFELgridGeneration: FastKernel tables generation
//

#include <math.h>
#include <string>

using namespace std;

namespace APFELgridGen {
  // Functions to generate the x-space grid a la APPLgrid
  static double m_transvar = 6.0;
  static double appl_fy(double x) { return - log(x) + m_transvar * ( 1 - x ); }
  static double appl_fx(double y);

  // Generate the FK table
  void generateFK(const string &agfile, const double &Q0, const double &mc, const double &mb, const double &mt, const double &AsRef, const double &QRef, const int &pto);
}
