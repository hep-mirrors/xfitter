/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/// Trapezoidal rule integration.
// ============================================
double TblInt1(
                double xh,   //!< step 
                double yi[], //!< values
                int np)      //!< # points
{
 double ys;
 int ii;

 if(np<2) return(0.0);
 ys=0.5*(yi[0]+yi[np-1]);
 for(ii=1;ii<np-1;ii++) ys +=yi[ii];
 ys *=xh;
 return(ys);
}

// ============================================
double TblInt2(double xi[], double yi[], int np)
{
 double ys;
 int ii;

 if(np<2) return(0.0);
 ys=0;
 for(ii=1; ii < np; ii++) ys += (yi[ii]+yi[ii-1])*(xi[ii]-xi[ii-1]);
 return 0.5*ys;
}
