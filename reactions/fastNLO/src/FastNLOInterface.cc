// Modified by: Marco Guzzi and Katerina Lipka
// DESY, 25/06/2014
// New interface for DiffTop is included


////////////////////////////////////////////////////////////////////////
//
//   FastNLOInterface
//                                                                      //
//  The interface through which fortran based xfitter interacts with   //
//  c++ version of FastNLOReader.                                       // 
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <FastNLOxFitter.h>
#include <cmath>
#include <string>


#define INTERP_PTS 4 //(Polynomial 7th order)
#define NGausspoints 25
#define Pi 3.14159265359

using namespace std;
typedef vector<bool> BoolArray; 


double Function_Mu(double s1, double s2);
double polint(double xa[],double ya[],int n,double x);
double interpC(double *A,double xx, double *x, int NGrid);
void gauleg(double x1,double x2,double *x,double *w, int n);


extern "C" {
  int fastnloinit_(const char *s, const int *idataset, const char *thfile, int *I_FIT_ORDER, bool *PublicationUnits , double* murdef, double* murscale, double *mufdef, double* mufscale);
  int fastnlocalc_(const int *idataset, double *xsec);
  int fastnlocalctop_(const int *idataset, double *xsec, double *thbin, double *tot, int *Npt);
  int fastnlopointskip_(const int *idataset, int *point, int *npoints);
  int hf_errlog_(const int* ID, const char* TEXT, long length);
  int hf_stop_();
  double interp_(double *A, double *xx1, double *x, int *NGrid1, double *res);
  int setfastnlotoppar_(const int *idataset);
}






map<int, FastNLOxFitter*> gFastNLO_array;
map<int, BoolArray*>     gUsedPoints_array;
int CreateUsedPointsArray(int idataset, int npoints);

int fastnloinit_(const char *s, const int *idataset, const char *thfile, int *I_FIT_ORDER, bool *PublicationUnits , double* murdef, double* murscale, double *mufdef, double* mufscale) {

  
   map<int, FastNLOxFitter*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   if(FastNLOIterator != gFastNLO_array.end( )) {
     int id = 12032301;
     const char* text = "I: Double initialization of the same fastnlo data set!";
     hf_errlog_(&id, text, (long)strlen(text));
     //hf_stop_();
     return 1;
   }
   
   FastNLOxFitter* fnloreader = new FastNLOxFitter( thfile );  
   
   if(*PublicationUnits)
     fnloreader->SetUnits(fastNLO::kPublicationUnits);
   else 
     fnloreader->SetUnits(fastNLO::kAbsoluteUnits);
     
   if(*murdef>=0.)
     fnloreader->SetMuRFunctionalForm((fastNLO::EScaleFunctionalForm) ((int) (*murdef)));
   if(*mufdef>=0.)
     fnloreader->SetMuFFunctionalForm((fastNLO::EScaleFunctionalForm) ((int) (*mufdef)));

   fnloreader->SetScaleFactorsMuRMuF(  *murscale, *mufscale);
    
   
   int ordercalc = *I_FIT_ORDER;

   if (ordercalc==3)  {
      fnloreader->SetContributionON(fastNLO::kFixedOrder,0,true);
      fnloreader->SetContributionON(fastNLO::kFixedOrder,1,true);
      fnloreader->SetContributionON(fastNLO::kFixedOrder,2,true);
      printf("DiffTop pert. order = NNLO O(alphas^4) ordercalc = %d\n",ordercalc);
   }
   else if (ordercalc==2) {
      fnloreader->SetContributionON(fastNLO::kFixedOrder,0,true);
      fnloreader->SetContributionON(fastNLO::kFixedOrder,1,true);
      fnloreader->SetContributionON(fastNLO::kFixedOrder,2,false);
      printf("DiffTop pert. order = NLO O(alphas^3) ordercalc = %d\n",ordercalc);

   }
   else if (ordercalc==1) {
      fnloreader->SetContributionON(fastNLO::kFixedOrder,0,true);
      fnloreader->SetContributionON(fastNLO::kFixedOrder,1,false);
      fnloreader->SetContributionON(fastNLO::kFixedOrder,2,false);
      printf("DiffTop pert. order = LO O(alphas^2) ordercalc = %d\n",ordercalc);
   }
   else {
      printf("DiffTop pert. order is not defined, ordercalc = %d:\n",ordercalc);
      exit(1);
   }
   
   gFastNLO_array.insert(pair<int, FastNLOxFitter*>(*idataset, fnloreader) );
   return 0;
}

int setfastnlotoppar_(const int *idataset) {
   //!< Dedicated settings for difftop
   map<int, FastNLOxFitter*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   map<int, BoolArray*>::const_iterator UsedPointsIterator = gUsedPoints_array.find(*idataset);
   if(FastNLOIterator == gFastNLO_array.end( )) {
     int id = 12032303;
     char text[256];
     sprintf(text, "S: Can not find FastnloReader for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }
   
   FastNLOxFitter* fnloreader = FastNLOIterator->second;   
   fnloreader->SetExternalFuncForMuF( &Function_Mu );
   fnloreader->SetExternalFuncForMuR( &Function_Mu);
   //fnloreader->SetScaleFactorsMuRMuF(1.0,1.0); //Be reminded that muR and muF scales are hard coded (that's not true!)

   return 0;
}

int fastnlocalc_(const int *idataset, double *xsec) {
  
   map<int, FastNLOxFitter*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   map<int, BoolArray*>::const_iterator UsedPointsIterator = gUsedPoints_array.find(*idataset);
   if(FastNLOIterator == gFastNLO_array.end( )) {
     int id = 12032302;
     char text[256];
     sprintf(text, "S: Can not find FastnloReader for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }
   
   FastNLOxFitter* fnloreader = FastNLOIterator->second;

   if(UsedPointsIterator == gUsedPoints_array.end( )) 
     CreateUsedPointsArray(*idataset, fnloreader->GetNObsBin());
   UsedPointsIterator = gUsedPoints_array.find(*idataset);
   
   if(UsedPointsIterator == gUsedPoints_array.end( )) {
     int id = 12032303;
     char text[256];
     sprintf(text, "S: Can not find proper UsedPointsIterator for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }

   BoolArray*     usedpoints = UsedPointsIterator->second;


   fnloreader->FillAlphasCache();
   fnloreader->FillPDFCache();			// pdf is 'external'! you always have to call FillPDFCache();
   fnloreader->CalcCrossSection();
 

  
   vector < double > xs = fnloreader->GetCrossSection();
 
   int outputidx = 0;
   for ( unsigned i=0;i<xs.size();i++){
     if(usedpoints->at(i)) {
       xsec[outputidx] = xs[i];
       outputidx++;
     }
   }
 

   return 0;
}

//MK14 New function for Difftop calculation: it is called in trunk/src/difftop_fastnlo.f
int fastnlocalctop_(const int *idataset, double *xsec, double *thbin, double *tot, int *Npt){
  
   map<int, FastNLOxFitter*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   map<int, BoolArray*>::const_iterator UsedPointsIterator = gUsedPoints_array.find(*idataset);
   if(FastNLOIterator == gFastNLO_array.end( )) {
     int id = 12032302;
     char text[256];
     sprintf(text, "S: Can not find FastnloReader for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }
   
   FastNLOxFitter* fnloreader = FastNLOIterator->second;

   if(UsedPointsIterator == gUsedPoints_array.end( )) 
     CreateUsedPointsArray(*idataset, fnloreader->GetNObsBin());
   UsedPointsIterator = gUsedPoints_array.find(*idataset);
   
   if(UsedPointsIterator == gUsedPoints_array.end( )) {
     int id = 12032303;
     char text[256];
     sprintf(text, "S: Can not find proper UsedPointsIterator for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }

   BoolArray*     usedpoints = UsedPointsIterator->second;

   //Perturbative order of the Calculation/Fit: LO, NLO, NNLO => 1,2,3

   fnloreader->FillAlphasCache();
   fnloreader->FillPDFCache();	  // pdf is 'external'! you always have to call FillPDFCache();
   fnloreader->CalcCrossSection();
   //fnloreader->PrintCrossSections();


   vector < double > xs = fnloreader->GetCrossSection();


   //Read the number of theory points 
   int Nthpoints;

   Nthpoints = xs.size();
   *Npt = Nthpoints;


   //Assign the values of the theory variable (PT or Y)
   double bin[Nthpoints];
   for ( unsigned int b = 0 ; b < xs.size() ; b++ ) 
     {
       bin[b] = fnloreader->GetObsBinLoBound(b,0);
       thbin[b] = bin[b];
     }

   //Assign theory Xsec values 
   int outputidx = 0;
   for (unsigned i=0; i<xs.size(); i++)
   {
     if(usedpoints->at(i)) 
      {
       xsec[outputidx] = xs[i];
       outputidx++;

      }
   }


   double xg[NGausspoints];
   double wg[NGausspoints];
   double Total;
   int k;

   //MK14 -> hard wired for ttbar FatsNLO Grids: PT spectrum starts at 1 GeV
   if(thbin[0] < 1.0)  thbin[0]=1.0;

   gauleg(thbin[0],thbin[Nthpoints-1],xg,wg,NGausspoints);

   Total=0.0;
   for (k=0; k<NGausspoints; k++)
     Total += interpC(xsec,xg[k],thbin,Nthpoints)*wg[k];

   *tot = Total;



   return 0;
}










int fastnlopointskip_(const int *idataset, int *point, int *npoints) {
  map<int, BoolArray*>::const_iterator UsedPointsIterator = gUsedPoints_array.find(*idataset);
  if(UsedPointsIterator == gUsedPoints_array.end( )) 
    CreateUsedPointsArray(*idataset, *npoints);

  UsedPointsIterator = gUsedPoints_array.find(*idataset);
  if(UsedPointsIterator == gUsedPoints_array.end( )) {
    int id = 12032304;
    char text[256];
    sprintf(text, "S: fastnlopointskip: Can not find proper UsedPointsIterator for DataSet: %d",*idataset);
    hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
  }
  
  BoolArray*     usedpoints = UsedPointsIterator->second;
  usedpoints->at(*point-1) = false;

  return 0;
}

int CreateUsedPointsArray(int idataset, int npoints) {
  cout << "creating new table..."<<endl;
  BoolArray* usedpoints = new BoolArray;
  for (int i=0; i<npoints; i++)
    usedpoints->push_back(true);
  gUsedPoints_array.insert(pair<int, BoolArray*>(idataset, usedpoints) );

  return 0;
}







//MK14 for scale settings: the top mass is assigned already in the FastNLO tabs
double Function_Mu(double s1, double s2) {
  //For dsig/dpT s1 corresponds to pT-> mT=sqrt(s1**2+mu**2)  
  //dsgima/dy still work is progress
   // --- fastNLO user: This is an example function
   //     to demonstrate how you might perform the
   //     definition of the scales using a
   //     'flexible-scale'-table, where a function
   //     of s1 and s2 can be used.
   //     Which variables s1 and s2 stand for are
   //     coded in the fastNLO table.
   double mu = 173.3;

   return mu;
}













void gauleg(double x1,double x2,double *x,double *w, int n)
{
//  (C) Copr. 1986-92 Numerical Recipes Software 

        int m,j,i;
        double eps=3.e-14,z1,z,xm,xl,pp,p3,p2,p1;
        m=(n+1)/2;
        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);

        for(i=1;i<=m;i++)
        {
                z=cos(Pi*((double)i-0.25)/((double)n+0.5));
                do
                {
                        p1=1.;
                        p2=0.;

                        for(j=1;j<=n;j++)
                        {
                                p3=p2;
                                p2=p1;
                                p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j;
                        }

                        pp=n*(z*p1-p2)/(z*z-1.);
                        z1=z;
                        z=z1-p1/pp;
                }
                while(fabs(z-z1)>eps);

                x[i-1]=xm-xl*z;
                x[n-i]=xm+xl*z;
                w[i-1]=2.*xl/((1.-z*z)*pp*pp);
                w[n-i]=w[i-1];
        }
}


double polint(double xa[],double ya[],int n,double x)
{
	int i,m,ns=0;
	double den,dif,dift,ho,hp,w,res;
	double c[n],d[n];

	dif=fabs(x-xa[0]);
	for (i=0;i<=n-1;i++)
	{
		if ((dift=fabs(x-xa[i]))<dif)
		{
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	res=ya[ns--];

	for (m=1;m<n;m++)
	{
		for (i=0;i<=n-m-1;i++)
		{
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];

			if ((den=ho-hp)==0.0)
			{
				printf("Error in routine polint");
				exit(1);
			}

			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}

		res+=((2*ns<(n-1-m) ? c[ns+1] : d[ns--]));
	}

	return res;
}






double interpC(double *A, double xx, double *x, int NGrid)
{
  //	double res,aux;
        double x1[NGrid+1];
	int k,i;
	
	for (i=0;i<NGrid;i++)
	  x1[i]=x[i];

	if (xx < x1[0])
	{
		printf("too low x-value %15.9g\n",x1[0]);
		exit(1);
	}

	for (k=0;xx>=x1[k];k++);

	k-=INTERP_PTS;

	if (k<0) k=0;
	if (k>NGrid-2*INTERP_PTS) k=NGrid-2*INTERP_PTS;

   if(xx==x1[NGrid]) return 0;
	else

	return polint(&x1[k],&A[k],2*INTERP_PTS,xx);

}	




 //-------------------------------------------------------
 //For the Fortran

double interp_(double *A, double *xx1, double *x, int *NGrid1, double *res)
{

  double xx=*xx1;
  int NGrid=*NGrid1;
  double res1;

  double x1[NGrid+1];
  int k,i;
  
  for (i=0;i<NGrid;i++)
    x1[i]=x[i];
  
  if (xx < x1[0])
    {
      printf("too low x-value %15.9g\n",x1[0]);
      exit(1);
    }
  
  for (k=0;xx>=x1[k];k++);
  
  k-=INTERP_PTS;
  
  if (k<0) k=0;
  if (k>NGrid-2*INTERP_PTS) k=NGrid-2*INTERP_PTS;
  
  if(xx==x1[NGrid]) return 0;
  else
    {
      res1 = polint(&x1[k],&A[k],2*INTERP_PTS,xx);
      *res = res1;
      return 0;
    }
}	








