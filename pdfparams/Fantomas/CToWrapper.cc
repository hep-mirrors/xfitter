// ==============================================
//
//      File Name: fantomas.cc
//
//      Description: An interface to 
//
//      Date created: 01/28/2022
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
#include "metamorphCollection.h"
#include "CToWrapper.h"

// Declaration of global variables inside fantomas.cc
metamorphCollection *metacol;   
bool xFitterFantomas = false;

extern "C" void readfantosteer_()
// function that reads fantomas input steering card
// readfantosteer() will be called by xFitter in Fantomas_PdfParam::atStart().
// before PDFs are calculated.
{
  metacol=new metamorphCollection();
  metacol->SetVerbosity(0); //set verbosity=0/1 to suppress/print out metamorph diagnostics

  metacol->ReadCard();
  //pn25 Update metamorph modulators using the just read parameters
  metacol->UpdateMetamorphs();

  xFitterFantomas = true;
}

extern "C" void writefantoout_()
// function that writes fantomas output steering card.
// writefantosteer_() will be called by xFitter in maind.f
// after PDFs are calculated to return updated fantomas
// parameters.
// lk24 reduced write functions to one function to write all cards at once when chi2 improves in xFitter
{
  if (xFitterFantomas == true)
  {
    metacol->WriteCard();
  }
}


extern "C" void updatefantopars_(int &flavor, double *deltasin)
// the array a[] will be passed from xFitter into the metamorphCollection metacol->
// a[] will be the difference between the initial value and the new updated value
// for all Sc and Sm parameters. The called function will update all metamorph
// objects inside of metacol-> updatefantopars() will be called each time
// xFitter updates the minuit input values inside of Fantomas_PdfParam::atStart()
{
  metacol->UpdateParameters(flavor,deltasin);
}

extern "C" double fantopara_(int &flavor, double &x)
// function to be called that returns PDF value at a specified x-value for a given flavor.
// fantomaspara() will be called inside of Fantomas_PdfParam::operator()(double x) when 
// xFitter calculates the PDF for each flavor.
{
  double ftmp;
  ftmp = metacol->f(flavor,x);
  //std::cout << "f(" << x << ") inside fantomas.cc: " << ftmp << std::endl;
  return ftmp;
}

extern "C" double fantomellinmoment_(int &flavor, int &MellinPower, int npts)
{
  double momenttmp = metacol->MellinMoment(flavor,MellinPower,npts);
  return momenttmp;
}
