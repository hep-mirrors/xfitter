//Support of chi2 priors for metamorphCollection
#include "metamorphCollection.h"

using namespace std;

// Function definitions
inline double ReLU(double x){
  return max(0.0, x);
}//ReLU ----------------------------------------------------------------------

double metamorphCollection::Chi2prior()
/* calculate the extra chi2 penalty associated with the metamorphs or their
   parameters. The blocks below include penalties of various types that can
   be uncommented depending on the user's need
*/
{
  double chi2=0.0; 
  //Uncomment any of the blocks below or add your own block

  //Prior contributions from individual metamorphs
  /* for (int ifl =0; ifl < NMeta; ifl++)
    chi2 += MetaVector[ifl].Chi2prior();
  */
  
  //Penalties on the (quadrature) sum of values of control points, i.e.,
  //deviations from modulators
  /* double w1=0.1; double chi2tmp=0.0;
  for (int ifl =0; ifl < NMeta; ifl++)
    for (int j = maxSc; j < maxSc + Nm[ifl]+1; j++){
      //uncomment one:
      //chi2tmp += abs(Scm[i][j]); //ridge or
      chi2tmp += Scm[ifl][j]*Scm[ifl][j]; //lasso
    }//for int j
  
  chi2 += w1*chi2tmp;
  */
  
  //Penalties on Bernstein coefficients
  /*double w2=1.0; double chi2tmp=0.0;
  for (int ifl =0; ifl < NMeta; ifl++)
    for (int j = 0; j < Nm[ifl]+1; j++){
      double Ci = MetaVector[ifl].Cs(j);       
      //uncomment one:
      //chi2tmp += ReLU(abs(Ci)-10);      
      chi2tmp += abs(log(1+abs(Ci)));     
    }//for int j #2
  chi2 +=w2*chi2tmp;
  */

  return chi2;
  
}//Chi2prior-------------------------------------------------------------------

void metamorphCollection::Chi2Diagnostics(const string& outputDirectory)
{
  ofstream fantoCout;
  
  string fname=outputDirectory+"fantomas_chi2diagnostics.txt";
  fantoCout.open(fname, ofstream::out);
  
  if (!fantoCout.is_open()) {
    cerr << "Failed to open file " << fname << "\n";
    return;
  }                                                                                                  
  
  fantoCout << "Write your diagnostics for Chi2prior to this file" << endl;
  
  fantoCout.close();
} // metamorphCollection::Chi2Diagnostics ------------------------------------

