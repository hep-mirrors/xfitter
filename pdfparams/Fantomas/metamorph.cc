#include "metamorph.h"
#include <vector>

using namespace std;

double DefCarrier(double x, double *a)
// The default carrier function. The input parameters are the momentum fraction 
// x and a pointer to the array a with control parameters.
{
  return a[0] * pow(x, a[1]) * pow(1.0 - x, a[2]);
} //DefCarrier------------------------------------------------

long int bic(int Nm, int k)
// Binomial coefficient C^Nm_k
{
  if (k > Nm)
    return 0;
  if (k == 0 || k == Nm)
    return 1;
  
  // Recur
  return bic(Nm - 1, k - 1) + bic(Nm - 1, k);
} //bic---------------------------------------------------------

int ControlPoint::SetBoundary(const int MappingModeIn, const double fmIn, const double fpIn) 
{
  MappingMode = MappingModeIn;
  
  fm = fmIn;
  fp = fpIn;
  if (fm >= fp) 
  {
    cout << "STOP: wrong boundaries in ControlPoint. fm = " << fm 
              << " > fp = " << fp << endl;
    exit(3);
  }
  
  return 1;
} //ControlPoint::SetBoundary----------------------------------

ControlPoint::ControlPoint(const double xIn, double &psIn, const int MappingModeIn, const double fmIn, const double fpIn) 
{
  x = xIn;
  ps = &psIn;
  SetBoundary(MappingModeIn, fmIn, fpIn);
  
} //constructor ControlPoint::ControlPoint ->

double ControlPoint::s() 
{
  return *ps;
}//ControlPoint::s-------------------------------------------

double ControlPoint::f() 
{
  double sf = *ps, ftmp;
  
  switch (MappingMode) 
  {
    case 0: // no scaling
      ftmp = sf;
      break;
    //lk25 MappingMode != 0 is currently not supported
    /*case 1: //bounded linear scaling
      ftmp = (fp+fm)/2.0 + (fp-fm)*sf/2.0;
      break;
    case 2: //bounded softsign scaling
      // Given a variable s in the interval (-infinity, +infinity), returns
      // a value in the interval f1m <= f <= fp. Softsign is approximately
      // linear in the interval -1 < s < 1 and is equal to (fp+fm)/2 for s=0.
      ftmp = (fp+fm)/2.0 + (fp-fm)*sf/(1.0+fabs(sf))/2.0;
      break;
    */
    default:
      cout << "STOP in ControlPoint: mapping mode " << MappingMode 
                << " is not implemented" << endl;
      cout << "x = " << x << ", s = " << sf << endl;
      exit(2);
  } //switch MappingMode
 
  if (ftmp < fm || ftmp > fp) 
  {
      cout << endl << "The returned function is out of bounds in ControlPoint::f"
	   << endl;
      cout << "The functions and the bounds are " << fm << " < " << ftmp 
                << " < " << fp << endl;
  }

  return ftmp;
  
} //ControlPoint::f----------------------------------------------

metamorph::metamorph(const int NmIn, const double XsIn[], double SmIn[], double ScIn[], const double xPowerIn, const vector <double>& vstretch, double (*CarrierIn)(double, double *))
// Constructor: construct a metamorph class by reading Nm values from an 
// input array Xs and setting control point parameters to point to an array 
// of external variables Sm.
// Also provide a pointer to the array Sc containing the parameters for the 
// carrier function. DefCarrier is the 
// default carrier function. An optional input parameter is a pointer to an 
// alternative carrier function.   
{

  pSc = ScIn;
  Carrier = CarrierIn; // set the pointer to the carrier function, given by DefCarrier
                      //unless specified otherwise;

  if (xPower <=0 || xPower >= 1 || NmIn < 0){
    cerr << "STOP: unacceptable input parameters for a metamorph:"<< endl
	 <<"xPower, Nm = "<< xPowerIn << ", "<< NmIn << endl;
    exit(1);
  }

  Nm = NmIn;
  xPower = xPowerIn;
  if (vstretch.size()>=2)
    xstrmax=vstretch[1];
  if (vstretch.size()>=1)
    xstrmin=vstretch[0];

  //Fill in vector CP with control points at stretched x values (Ys)
  vector<double> Ys(Nm + 1);     // TH added VLA fix, default initialized to 0

  for (int icp = 0; icp < Nm+1; icp++){
    Ys[icp] = yx(XsIn[icp]); //stretched X-value for icp-th CP 
    CP.push_back(ControlPoint(Ys[icp], SmIn[icp]));
  }
              
  //factorial and binomial coefficeint calculation
  //matrices T, M
  M = new cl2DArray<double>(Nm+1, Nm+1);
  InvM = new cl2DArray<double>(Nm+1, Nm+1);
  T = new cl2DArray<double>(Nm+1, Nm+1);
  InvT = new cl2DArray<double>(Nm+1, Nm+1);
  multi = new cl2DArray<double>(Nm+1, Nm+1);
  
  for (int i = 0; i < Nm+1; i++)
    for (int j = 0; j < Nm+1; j++) 
    {
      (*T)(i,j) = pow(Ys[i], j);
      if (i <= j)
        (*M)(j,i) = pow(-1, j-i) * bic(Nm, i) * bic(Nm-i, Nm-j);
      else
        (*M)(j,i) = 0.0;
    } //for (int j=0...
  
  C = new double[Nm+1];
  P = new double[Nm+1];
  Pdec = new int[Nm+2];
  LUPDecompose(M, Nm+1, 1.0e-10, Pdec);
  LUPInvert(M, Pdec, Nm+1, InvM); //Calculate Inverse of M

  LUPDecompose(T, Nm+1, 1.0e-10, Pdec);
  LUPInvert(T, Pdec, Nm+1, InvT); //Calculate Inverse of T

} // constructor metamorph::metamorph-------------------------------

double metamorph::yx(const double x){
  double base = 1.0 / (1.0 / (pow(x, stretchPower)+
			      pow(xstrmin, stretchPower))
		       + (1.0 - pow(xstrmax, stretchPower))
		       / pow(xstrmax, stretchPower));
  double xstretched = pow(base, 1.0/stretchPower);
  return pow(xstretched, xPower);
}//metamorph::yx-----------------------------------------------------

void metamorph::SetXstretching(const double xstrminIn, const double xstrmaxIn, const double stretchPowerIn){
  xstrmin=xstrminIn; xstrmax=xstrmaxIn; stretchPower=stretchPowerIn;
}//metamorph::SetXstretching-------------------------------

void metamorph::GetXstretching(double xstrminOut, double xstrmaxOut, double stretchPowerOut){
  xstrminOut=xstrmin; xstrmaxOut=xstrmax; stretchPowerOut=stretchPower;
}//metamorph::GetXstretching-------------------------------

int metamorph::SetBoundary(const int MappingModeIn, double fmIn[], double fpIn[]) 
{
  int j = 0;
  for (vector<ControlPoint>::iterator it = CP.begin(); it != CP.end(); ++it) 
  {
    it->SetBoundary(MappingModeIn, fmIn[j], fpIn[j]);
    j++;
  }

  return 1;
} //metamorph::SetBoundary -------------------------------------

void metamorph::UpdateModulator() 
{
  double Ci;
  //The value of each vector P element is the f() value of each element of the 
  // CP vector divided by the carrier function at the x() value
  for (int ix = 0; ix < Nm+1; ix++) 
    P[ix] = CP[ix].f();
 
  for (int i = 0; i < Nm+1; i++)
    for (int j = 0; j < Nm+1; j++) 
    {
      (*multi)(i,j) = 0.0;
      for (int k = 0; k < Nm+1; k++)
        (*multi)(i,j) += (*InvM)(i,k) * (*InvT)(k,j);
    } //InvM.InvT ->

  for (int i = 0; i < Nm+1; i++) 
  {
    Ci = 0.0;
    for (int k = 0; k < Nm+1; k++)
      Ci += (*multi)(i,k) * P[k];
    C[i] = Ci;
  } //C=InvM.InvT.P->

  ModulatorSet = true;

} //UpdateModulator----------------------------------------------

double metamorph::Cs(const int i) 
// Return the value for the ith Bezier coefficient
{
  return C[i];
} //double C-----------------------------------------------------

double metamorph::Modulator(const double x) 
// Return the value of the Modulator function for the momentum fraction x
{
  // lk22 changed modulator function from B_(n,l)(x) to 1 + B_(n,l)(x) with B = 0 for Nm=0 
  double y = yx(x); //a stretched x parameter of the Bezier curve
  double mod = 1.0;
  if (Nm != 0) 
    for (int im = 0; im < Nm+1; im++)
      mod += C[im] * bic(Nm, im) * pow(y, im) * pow(1 - y, Nm-im);

  return mod;
} //double Modulator->

double metamorph::f(const double x) 
// Return the value of the metamorph for the momentum fraction x
{
  double ftmp;
  
  //double fp = 100;
  //double fm = 0;
  double smod = Modulator(x);
  double modtmp = smod;
  //double modtmp = (fp+fm)/2.0 + (fp-fm)*smod/(1.0+fabs(smod))/2.0;
  //double modtmp = exp(smod);
  ftmp = (*Carrier)(x, pSc) * modtmp;
  return ftmp;
} //f-----------------------------------------------------------------

double metamorph::GetMellinMoment(double MellinPower, int npts) 
// metamorph::GetMellinMoment() returns <x^(n+1) f>, i.e., the  integral of x^(n+1) * f(x) * dx
// over 0 < x < 1, where n=MellinPower. Integration is implemented based on
//function adxmoment by Jon Pumplin. Employs mapping z = (x + c*x**p)/(1+c)
//and simple integration method in z. 
{
  int maxpts = 20000;
  double c = 0.2;
  double p = 0.2;
  static int nold = 0.0;
  double one = 1.0;
  double x0, z0, u, z, wt, x, y, ypri, sum;
  static double xvec[20001], wvec[20001]; // arrays of size maxpts+1 for simplicity
  static int iMellinerr0=0, iMellinerr1=0;
  const double oneminuseps = -0.93; //Produce a warning if the integrand diverges as 

  // With the default carrier, check that the Mellin moment is integrable. Warn
  // if the integration is slowly converging
  if (Carrier==DefCarrier){
    
    if (pSc[1] + MellinPower < -1 || pSc[2] < -1) {
      cout << "STOP: a non-converging Mellin moment in metamorph "<< ID << endl;
      cout << "MellinPower, Sc[1], Sc[2] = " << MellinPower << ", "<< pSc[1] << ", " << pSc[2] << endl;
      exit(3);
    }
    
    if (pSc[1] + MellinPower < oneminuseps){
      if (iMellinerr0 <= 5){
	cout << "WARNING: unstable numerical integration for a Mellin moment in metamorph " << ID << endl;
	cout << "At x->0, the integrand grows as x^"<< pSc[1] + MellinPower << endl;
      }
      if (iMellinerr0 == 5)
	cout << "Too many warnings; printing to cout is suppressed." << endl;
      iMellinerr0++;
    } // if(pSc[1] + MellinPower < oneminuseps)
    
    if (pSc[2] < oneminuseps){
      if (iMellinerr1 <= 5){
	cout << "WARNING: unstable numerical integration for a Mellin moment in metamorph " << ID << endl;
	cout << "At x->1, the integrand grows as (1-x)^"<< pSc[2] << endl;
      }
      if (iMellinerr1 == 5)
	cout << "Too many warnings; printing to cout is suppressed." << endl;
      iMellinerr1 ++;
    } // if(pSc[2] < oneminuseps)
  }//if (Carrier=DefCarrier)
  
  // compute the points and weights on first call, or if npts has changed...
  if (npts != nold) 
  {
    if (npts < 10 || npts > maxpts) 
    {
      cout << "metamorph::GetMellinMoment:  fatal npts= " << npts << endl;
      exit(3);
    }
    nold = npts;
    // find z0 as point where x and c*x**p are equal...
    x0 = pow(c, (one / (one - p)));
    z0 = (x0 + c * pow(x0, p)) / (one + c);
    for (int i = 1; i <= npts; i++) 
    {
      // formerly used equal spacing in z...
      //            z = (i - 0.5d0)/float(npts)
      //            wt = one/float(npts)
      // changed to equal spacing in u = sqrt(1-z), which should allow 
      // mild singularity at z->1, i.e., at x->1. (jcp 12/02)
      u = (i - 0.5) / (double)npts;
      z = 1 - pow(u, 2);
      wt = 2 * u / (double)npts;
      
      // get starting estimate by asymptotic forms...
      if (z < z0) 
        x = pow((z * (one + c) / c), (one / p));
      else 
        x = one - (one - z) * (one + c) / (one + p * c);
      
      // solve z = (x + c*x**p)/(1 + c) for x by Newton's method...
      for (int iter = 1; iter <= 9; iter++) 
      {
        y = x + c * pow(x, p) - (one + c) * z;
        ypri = one + c * p * pow(x, (p - one));
        x = x - y / ypri;
        if ((iter > 5) && (fabs(y) < 1.e-10)) 
          goto label;
      }
      cout << " " << "metamorph::GetMellinMoment: fatal convergence " << y << " " << z << " " << x << endl;
      exit(3);
    label:
      xvec[i] = x;
      wvec[i] = wt * (one + c) / (one + c * p * pow(x, (p - one)));
    }//for (int i=1...
  }//if (npts != nold...

  // do the integral...
  sum = 0;
  for (int i = 1; i <= npts; i++) 
  {
    x = xvec[i];
    sum = sum + pow(x, MellinPower) * f(x) * wvec[i];
  }
  
  return sum;
} // metamorph::GetMellinMoment-------------------------------------------

double metamorph::GetConditionNumber() 
{
  double abst = 0, absInvt = 0, T2 = 0, InvT2 = 0,
    Tnorm = 0, InvTnorm = 0, CondNum = 0;
  
  for (int i = 0; i < Nm+1; i++)
    for (int j = 0; j < Nm+1; j++) 
    {
      abst = fabs((*T)(i,j));
      absInvt = fabs((*InvT)(i,j));
      T2 += pow(abst, 2);
      InvT2 += pow(absInvt, 2);
    }
  Tnorm = sqrt(T2);
  InvTnorm = sqrt(InvT2);
  CondNum = Tnorm * InvTnorm;
  return CondNum;
}// metamorph::GetConditionNum-------------------------------------

metamorph::~metamorph() 
{
  delete[] P;
  delete[] C;
  delete[] Pdec;
  
  delete M;
  delete InvM;
  delete T;
  delete InvT;
  delete multi;
  
} //~metamorph----------------------------------------------------
