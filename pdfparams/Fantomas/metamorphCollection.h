/*
Description: class metamorphCollection provides a container of metamorph objects with parametrizations
 for parton distributions or other nonperturbative functions realized using B'ezier curves. This class
 also provides methods to initialize and update the member metamorph functions during a global fit.
*/

#ifndef METAMORPHCOLLECTION_H
#define METAMORPHCOLLECTION_H

using namespace std;
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <array>
#include <string>
#include <math.h>
#include <map>
#include <string>
#include <cstring>
#include <algorithm>
#include "metamorph.h"
#include <unistd.h> // For readlink()
#include <limits.h>  // For PATH_MAX
#include <ctime> // For time and localtime
#include <iomanip>


class metamorphCollection
// metamorphCollection is a container for Nmeta metamorph objects providing parametrizations for
// parton distributions or analogous nonperturbative functions. metamorphCollection also provides methods
// to read, update, and return PDF values for a metamorph object corresponding to specified flavors
// in the provided input card. A metamorph member is created while reading its parameters from a
// steering cards. MetaVector is a vector containing the member metamorph objects. 
// MetaRoster is a map to access members of MetaVector according to the integer PDF flavor iflavor provided
// from the external program. 
{
private:
  
  //Maximal sizes of fixed-size interface arrays; adjust as necessary
  const static int maxSc = 3;                   // maximal number of Sc variables for each flavor
  const static int maxNm=7;               // maximal degree of the Bezier polynomial
  const static int maxctrlpts = maxNm+1;      //maximal number of control points
  const static int maxScm = maxSc + maxctrlpts; // maximal allowed values in array Scm
  const static int maxMet = 10;                 // maximal number of metamorph flavors allowed 
  
  vector<metamorph> MetaVector;    // Vector containiner for all metamorph objects 
  map<int, metamorph*> MetaRoster; // MetaRoster maps the integer physical flavor ifl onto the corresponding metamorph in MetaVector
  map<int, int> PositionRoster;  // Map of physical flavor ifl onto the position iMet of the corresponding metamorph in MetaVector 

  int NMeta = 0;          // Number of metamorph members in MetaCollection
  int paraiMet = 0;      // variable used to select Scm parameters in metamorphCollection::UpdateParameters()

  int iflavor[maxMet] = {0};  
  int Nm[maxMet] = {0};
  int MappingMode[maxMet] = {0};
  double xPower[maxMet] = {0};
  double Xs0[maxMet][maxctrlpts] = {{0}}, fp0[maxMet][maxctrlpts] = {{0}}, fm0[maxMet][maxctrlpts] = {{0}};
  double Xs[maxMet][maxctrlpts] = {{0}}, fp[maxMet][maxctrlpts] = {{0}}, fm[maxMet][maxctrlpts] = {{0}};
  string strScm[maxMet][maxScm] = {{}};
  double Scm0[maxMet][maxScm] = {{0}},  //initial parameters of metamorphs
    Scm[maxMet][maxScm] = {{0}}, //current parameters of metamorphs
    DSwitch[maxMet][maxScm]={{0}}; //array of switches to turn metamorph parameter changes on or off
 
  const  string newflag = "NEW",   // flag signaling to calculate PDF and use as control point for metamorph PDF
         calcflag = "CALC",    // flag signaling to calculate PDF in the output card
         fixflag = "FIX";   // flag signaling to use metamorph::Carrier(x) for control point in PDF calculation
  int k;                         // counts the number of input control points for a given flavor
  int l;                         // counts the number of control points used in a given metamorph
  int m;                         // counts the number of control points for temporary metamorph
  int iPts[maxMet];              // array containing # of entries for each flavor 

  array<array<string, 4>, maxMet> flvcomment{}; // TH2025 revised; array to store headers for each flavor in fantomas card
  unsigned int VerbosityLevel=0; //Set VerbosityLevel=1 to print out diagnostic messages

  void PushMember();
  // MetamorphCollection::PushMember() is called inside metamorphCollection::ReadCard() to push the
  // newly created metamorph into MetaVector and list in two "address books", MetaRoster and PositionRoster.
  // The initial values read from metamorphCollection::ReadCard() are used to create an initial metamorph
  // object inside the vector MetaVector. The boundary conditions are defined as well for the 
  // MetaVector element and then the Modulator function is calculated, allowing for the initial 
  // PDF values to be calculated. The vector element is stored inside of MetaRoster,
  // mapping the flavor and the metamorph object which will be used to easily access the PDF value
  // of any given flavor from the card. Scm is used to calculate the metamorph objects
  // instead of Scm0, since Scm will be updated by metamorphCollection::UpdateParams() when xFitter
  // varies the parameters Sc and Sm.

public:
  metamorphCollection();
  
  void ReadCard(const string& inputcard="steering_fantomas.txt");
  // metamorphCollection::ReadCard reads a fantomas steering card file providing initial values,
  // Sc0, Sm0, etc., for varied parameters Sc, Sm, etc. of the initialized metamorphs. In the course
  // of the fit, the metamorph parameters are updated using changes (deltas) provided by
  // the external fitting program, as described in the header of UpdateParameters.
  // 
  // Input: inputcard = name of the input steering card; default: steering_fantomas.txt

  void WriteCard(const string& outputcard="steering_fantomas_out.txt");
  // Create an output card for Fantomas using the updated parameters of Sc and Sm.
  // Output: outputcard = name of the output steering card; default: steering_fantomas_out.txt
  
  void UpdateParameters(const int ifl, double *deltas);
  // Must be called to update metamorph parameters every time the external fitting program updates
  // the changes (deltas) in these parameters.
  // Inputs: ifl= the physical flavor ID of the metamorph to be updated
  //           deltas=an array of changes in the metamorph parameters, Sc and Sm. 
  // For a metamorph with a Bezier polynomial degree Nm, given a vector deltas[maxSc+Nm+1] from
  // the fitting program, the parameters are recalculated as
  // Sc[i]=Sc0[i]+delta[i], i=0,.., maxSc-1; Sm[j]=Sm0[j]+delta[j+maxSc], j=0,...,Nm,
  // where Sc0[i] and Sm0[j] are constant initial values. The initial value for each fitted delta
  // is set to 0, so that only the change is added. 

  void UpdateMetamorphs();
  //Update modulators for all metamorph members

  double f(const int ifl, const double x);
  // Returns the value of metamorph of flavor ifl and momentum fraction x.

  double MellinMoment(int ifl, double MellinPower, int npts=10000);
  //Returns the Mellin moment <x^(MellinPower+1) f> =\int_0^1 x^(MellinPower+1) f(ifl, x) dx
  //for the metamorph of flavor ifl

  double GetConditionNumber(int ifl);
  //Returns the condition number for metamorph of flavor ifl.

  int GetMetamorphCount();
  //returns the number of metamorph members in MetamorphCollection

  void SetVerbosity(const unsigned int VerbosityLevelIn=0){
  //Sets verbosity level for printing out diagnostics  
    VerbosityLevel=VerbosityLevelIn;
  }

  ~metamorphCollection();
}; // class metamorphCollection

#endif
