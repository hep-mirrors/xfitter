// ==============================================
//
//      File Name: FantoWrapper.cc
//
//      Description: A C-style wrapper to the Fantomas package
//
// ==============================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>
#include "metamorphCollection.h"
#include "FantoWrapper.h"

//lk25
#define XFITTER

// Declaration of global variables inside fantomas.cc
metamorphCollection *metacol;   
bool SteeringCardExists = false;

void readfantosteer()
// function that reads fantomas input steering card
// readfantosteer() will be called by xFitter in Fantomas_PdfParam::atStart().
// before PDFs are calculated.
{
  metacol=new metamorphCollection();
#ifndef FANTOMAS_XFITTER  
  metacol->SetVerbosity(1);
#endif //FANTOMAS_XFITTER
  
  metacol->ReadCard();

  SteeringCardExists = true;
}

void writefantoout()
// function that writes fantomas output steering card.
// writefantosteer_() will be called by xFitter in maind.f
// after PDFs are calculated to return updated fantomas
// parameters.
{
  if (SteeringCardExists)
  {
    metacol->WriteCard();
#ifdef FANTOMAS_XFITTER
    metacol->MetamorphDiagnostics("output/");
    int ret = system("cp pdfparams/Fantomas/metamorphCollectionPrior.cc output/");
    if (ret != 0) {
        std::cerr << "Failed to copy metamorphCollectionPrior.cc to output/\n";
    }
    metacol->Chi2Diagnostics("output/");
#else    
    metacol->MetamorphDiagnostics("./");
#endif //FANTOMAS_XFITTER   
  }
}


void updatefantopars(int &flavor, double *parsin)
// the array a[] will be passed from xFitter into the MetamorphCollection metacol->
// a[] will be the difference between the initial value and the new updated value
// for all Sc and Sm parameters. The called function will update all metamorph
// objects inside of metacol-> updatefantopars() will be called each time
// xFitter updates the minuit input values inside of Fantomas_PdfParam::atStart()
{
  metacol->UpdateParameters(flavor,parsin);
}

double fantopara(int &flavor, double &x)
// function to be called that returns PDF value at a specified x-value for a given flavor.
// fantomaspara() will be called inside of Fantomas_PdfParam::operator()(double x) when 
// xFitter calculates the PDF for each flavor.
{
  double ftmp;
  ftmp = metacol->f(flavor,x);
  //std::cout << "f(" << x << ") inside fantomas.cc: " << ftmp << std::endl;
  return ftmp;
}

double fantoMellinMoment(int &flavor, int &MellinPower, int npts)
{
  double momenttmp = metacol->MellinMoment(flavor,MellinPower,npts);
  return momenttmp;
}

void getfantochi2(double& fantochi2)
// to be called for in fcn.f when chi2out is calculated.
// function will set argument to chi2 penalty from fantomas constrant
// by calling function in MetamorphCollection.h
{
    fantochi2 = metacol->Chi2prior();  
}//getfantochi2 ->

extern "C"
{
  void readfantosteer_()
  {readfantosteer();}
  
  void writefantoout_()
  {writefantoout();}
  
  void updatefantopars_(int &ifl,double *parsin)
  {updatefantopars(ifl,parsin);}
  
  double fantopara_(int &ifl, double &x)
  {return fantopara(ifl, x);}
  
  double fantomellinmoment_(int &ifl, int &mellinpower, int npts)
  {return fantoMellinMoment(ifl, mellinpower,npts);}
  
  void getfantochi2_(double& fantochi2)
  {getfantochi2(fantochi2);}

}
