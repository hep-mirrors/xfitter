/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _gZDPDF_HDR
#define _gZDPDF_HDR

#include "hdpdf.h"

#define MAX_N_ERR 32
/*
extern hdpdf_t* gZDPDF;
extern int CurGridIndex;
extern bool DeltaMode;
*/

hdpdf_t* CurGrid();
void LoadGrid(const string& Label, int n=0);
void SelectGrid(int n);
void SetDeltaMode(int dm);
void Set_tmin(double t);
void ZeusDpdf3Pom(double xP, double zP, double QQ, double f[7], int xpow=1);
void ZeusDpdf3Reg(double xP, double zP, double QQ, double f[7], int xpow=1);
void ZeusDpdf3(double xP, double zP, double QQ, double f[7], int xpow=1);

#endif
