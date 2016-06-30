//
// APFELgridGeneration: FastKernel tables generation
//

#include "APFELgridGeneration.h"

#include <iostream>
#include <fstream>

#include "APFEL/APFEL.h"
#include "APFELgrid/APFELgrid.h"
#include "APFELgrid/fastkernel.h"
#include "appl_grid/appl_grid.h"

namespace APFELgridGen {
  // Generate the FK table
  void generateFK(const string &fkfile, const double &Q0, const double &mc, const double &mb, const double &mt, const double &AsRef, const double &QRef, const int &pto) {

    // Read the APPLgrid through the usual procedure.
    // It is assumed that the APPLgrid file to be used has the same name of the FK table
    // with the only difference that the extension is .root rather that .fk.
    // Replace the extension .fk with .root in input file name of the FK table.
    string agfile = fkfile;
    agfile.replace(agfile.end()-3, agfile.end(), ".root");
    const appl::grid g(agfile);

    // Generate the external x-space grid
    const int nx = 30;
    const double ymin = appl_fy(0.99*APFELgrid::get_appl_Xmin(g, true));
    const double ymax = appl_fy(1.0);
      
    // Populate grid
    double *xg = new double[nx+1];
    for (int i=0; i<=nx; i++) xg[i] = appl_fx( ymin + ( ( ymax - ymin ) / ( (double) nx ) ) * i );

    // Initialize APFEL with the input parameters
    APFEL::SetTheory("QCD");
    APFEL::SetPoleMasses(mc,mb,mt);
    APFEL::SetAlphaQCDRef(AsRef,QRef);
    APFEL::SetPerturbativeOrder(pto);
    APFEL::SetNumberOfGrids(1);
    APFEL::SetExternalGrid(1,nx,5,xg);

    // Compute the FK table
    //NNPDF::FKTable<double>* FK = APFELgrid::computeFK(Q0,fkfile,g,agfile);
    NNPDF::FKTable<double>* FK = APFELgrid::computeFK(Q0,"No name",g,agfile);

    // Add the relevant tags to the FK table
    FK->AddTag(NNPDF::FKHeader::THEORYINFO, "PerturbativeOrder", pto);
    FK->AddTag(NNPDF::FKHeader::THEORYINFO, "MCharm", mc);
    FK->AddTag(NNPDF::FKHeader::THEORYINFO, "MBottom", mb);
    FK->AddTag(NNPDF::FKHeader::THEORYINFO, "MTop", mt);
    FK->AddTag(NNPDF::FKHeader::THEORYINFO, "AlphasRef", AsRef);
    FK->AddTag(NNPDF::FKHeader::THEORYINFO, "QRef", QRef);

    // Write the FK table to file
    ofstream outfile;
    outfile.open(fkfile.c_str());
    FK->Print(outfile);
    outfile.close();

    // Cleaning
    APFEL::CleanUp();
    delete[] xg;
    return;
  }

  // Functions to generate the x-space grid a la APPLgrid
  static double appl_fx(double y) {
    // use Newton-Raphson: y = ln(1/x)
    // solve   y - yp - a(1 - exp(-yp)) = 0
    // deriv:  - 1 -a exp(-yp)

    if ( m_transvar==0 )  return exp(-y);
    
    const double eps  = 1e-12;  // our accuracy goal
    const int    imax = 100;    // for safety (avoid infinite loops)
    
    double yp = y;
    double x, delta, deriv;
    for ( int iter=imax ; iter-- ; ) {
      x = exp(-yp);
      delta = y - yp - m_transvar*(1-x);
      if ( fabs(delta)<eps ) return x; // we have found good solution
      deriv = -1 - m_transvar*x;
      yp -= delta / deriv;
    }
    // exceeded maximum iterations
    cerr << "_fx2() iteration limit reached y = " << y << endl;
    cout << "_fx2() iteration limit reached y = " << y << endl;
    return exp(-yp);
  }
}
