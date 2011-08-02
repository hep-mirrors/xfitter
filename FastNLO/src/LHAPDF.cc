// Author: Krzysztof Nowak
// DESY, 01/08/2011

//  Version 0.1, 

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Dummy LHAPDF interface                                              //
//                                                                      //
//  The interface through which fortran based h1fitter interacts with   //
//  c++ version of FastNLOReader.                                       // 
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <LHAPDF.h>

double LHAPDF::alphasPDF(double d1) { return 0.; }
void LHAPDF::setVerbosity(Verbosity noiselevel) {return; }
void LHAPDF::initPDFSetByName(const std::string& filename) {return;}
int LHAPDF::numberPDF() {return 0;}
void LHAPDF::initPDF(int memset) {return;}
std::vector<double> LHAPDF::xfx(double x, double Q) {std::vector<double> k; return k;}
