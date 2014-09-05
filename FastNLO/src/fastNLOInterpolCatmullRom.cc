// Author: Daniel Britzger
// DESY, 27/06/2013

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "fastNLOInterpolCatmullRom.h"

using namespace std;

//______________________________________________________________________________


fastNLOInterpolCatmullRom::fastNLOInterpolCatmullRom(double min, double max) : fastNLOInterpolBase(min,max,4) {
   debug["fastNLOInterpolCatmullRom"]<<"New fastNLOInterpolCatmullRom instance."<<endl;
}


//______________________________________________________________________________


fastNLOInterpolCatmullRom::~fastNLOInterpolCatmullRom(void) {
}


//______________________________________________________________________________

// vector<pair<int,double> > fastNLOInterpolCatmullRom::GetCopyOfNodeValues(double x) {
//    vector<pair<int,double> > ret = GetNodeValues(x);
//    return ret;
// }

//vector<pair<int,double> > fastNLOInterpolCatmullRom::CalcNodeValues(double x) {
void fastNLOInterpolCatmullRom::CalcNodeValues(vector<pair<int,double> >& nodes, double x) {
   // Performs interpolation of value value on grid 'fgrid'.
   // uses distance measure 'fdm'
   // returns for for all relevant grid points
   //    the integer number of that node
   //    the vale, that this node obtains.
   //

   static const unsigned int nS = 4; // number of nodes that receive contributions from this interpolation

   // --- relative distance delta - in function fdm H(x)
   // deltascale (Interpol(.,.,.delta,.): relative distance of value to node 'nnode'
   double delta = GetDelta(x);
   // --- distances to all nodes
   double  dist[nS] = {1.+delta,0.+delta,1.-delta,2.-delta} ;

   // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
   int nmod = 0;                      // --- variable for final node
   int nnode = FindLargestPossibleNode(x);
   int nmax = fgrid.size()-2;
   if (fLastGridPointWasRemoved) nmax=fgrid.size()-1;


   // nnode: number of the next node to the left of the current value
   // nmax:  number of the last node which could lie to the left of a potential value
   // delta: relative distance of value to node 'nnode'
   // ikern: select interpolation kernel   1:Catmul Rom  2: Lagrange
   // nmod:  modified number of next node to the left (to be used for storage - relevant only at boundaries)
   // kernel: array(4) containing the interpolation kernel

   static vector <double> kern(nS);
   // --- Catmul Rom interpolation kernel
   if (nnode == 0 ) { // --- left boundary
      kern[0] = 1.0 - 7.0/6.0*delta - 1.0/6.0*delta*delta + 1.0/3.0*delta*delta*delta;
      kern[1] = 4.0/3.0*delta + 1.0/3.0*delta*delta - 2.0/3.0*delta*delta*delta;
      kern[2] = -1.0/6.0*delta - 1.0/6.0*delta*delta + 1.0/3.0*delta*delta*delta;
      kern[3] = 0.0;
      nmod = nnode + 1;
   } else if (nnode == nmax) { // --- right boundary
      kern[0] = 0.0;
      kern[1] = -1.0/6.0*dist[2] - 1.0/6.0*dist[2]*dist[2] + 1.0/3.0*dist[2]*dist[2]*dist[2];
      kern[2] = 4.0/3.0*dist[2] + 1.0/3.0*dist[2]*dist[2] - 2.0/3.0*dist[2]*dist[2]*dist[2];
      kern[3] = 1.0 -7.0/6.0*dist[2] -1.0/6.0*dist[2]*dist[2] +1.0/3.0*dist[2]*dist[2]*dist[2];
      nmod = nnode - 1;
   } else { // --- central region
      kern[0] = 2.0 - 4.0*dist[0] + 2.5*dist[0]*dist[0] - 0.5*dist[0]*dist[0]*dist[0];
      kern[1] = 1.0 - 2.5*dist[1]*dist[1] + 1.5*dist[1]*dist[1]*dist[1];
      kern[2] = 1.0 - 2.5*dist[2]*dist[2] + 1.5*dist[2]*dist[2]*dist[2];
      kern[3] = 2.0 - 4.0*dist[3] + 2.5*dist[3]*dist[3] - 0.5*dist[3]*dist[3]*dist[3];
      nmod = nnode;
   }

   // keep value for next time.
   fnmod = nmod;

   //    cout<<"   x="<<x<<"\tdelta="<<delta<<"\tnnode="<<nnode<<"\tk[3]="<<kern[3]<<"\tg-max="<<fgrid.back()<<"\tnx="<<fgrid.size()<<"\tnf="<<fnmod<<"\tx0b="<<nmod-1<<endl;

   // generate return values
   //vector<pair<int,double> >  ret(nS);


   nodes.resize(nS);
   for ( unsigned int i = 0 ; i<nS ; i++ ){
      nodes[i] = make_pair(nmod-1+i,kern[i]);
   }

   if (fLastGridPointWasRemoved ) {
      if ( nodes.back().first==(int)fgrid.size() ) {
         nodes.resize(3);
         if ( nodes.back().first==(int)fgrid.size() ) {
            nodes.resize(2);
            if ( nodes.back().first==(int)fgrid.size() ) {
               nodes.resize(1);
               if ( nodes.back().first==(int)fgrid.size() ) {
                  nodes.resize(0);
               }
            }
         }
      }
   }
}
