// ==============================================
//
//      File Name: fantomas.cc
//
//      Description: Calculates PDFs using 
//                   metamorph.h to be called
//                   by xfitter
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
#include "fantomas.h"

//lk22
#define XFITTER

// Declaration of global variables inside fantomas.cc
metamorphCollection *metacol;   
bool xFitterFantomas = false;

void readfantosteer()
// function that reads fantomas input steering card
// readfantosteer() will be called by xFitter in Fantomas_PdfParam::atStart().
// before PDFs are calculated.
{
  metacol=new metamorphCollection();

  metacol->ReadCard();

  xFitterFantomas = true;
}

void writefantoout(double& chi2in)
// function that writes fantomas output steering card.
// writefantosteer_() will be called by xFitter in maind.f
// after PDFs are calculated to return updated fantomas
// parameters.
// lk24 reduced write functions to one function to write all cards at once when chi2 improves in xFitter
{
  if (xFitterFantomas == true)
  {
    metacol->WriteCard();
    //metacol->WriteC();
    //metacol->Writechi2log(chi2in);
    //metacol->Writeerrlog();
  }
}

extern "C" void writefantoout_(double& chi2in)
{
  writefantoout(chi2in);
}//writefantosteer_ ->


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
  double fantomoment;
  double momenttmp = metacol->MellinMoment(flavor,MellinPower,npts);
  return momenttmp;
}

//lk23 added function to get extra chi2 from condition of 0<Ci~1 where Ci is bezier coefficients.
//void getfantochi2(double& fantochi2)
// to be called for in fcn.f when chi2out is calculated.
// function will set argument to chi2 penalty from fantomas constrant
// by calling function in MetamorphCollection.h
//{
//  if (xFitterFantomas == true)
//    fantochi2 = metacol->fantomasconstrchi2();  
//}//getfantochi2 ->

//extern "C" void getfantochi2_(double& fantochi2)
//{
//  getfantochi2(fantochi2);
//}//getfantochi2_ ->
