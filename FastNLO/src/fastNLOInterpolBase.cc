// Author: Daniel Britzger
// DESY, 27/06/2013

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "fastNLOInterpolBase.h"

using namespace std;

//______________________________________________________________________________


fastNLOInterpolBase::fastNLOInterpolBase(double min, double max, int nMinNodes = -1) :
   PrimalScream("fastNLOInterpol"),fNMinNodes(nMinNodes), fvalmin(min), fvalmax(max) {
   debug["fastNLOInterpolBase"]<<"New fastNLOInterpolBase instance."<<endl;
   fLastVal = -34729.432;
   fLastGridPointWasRemoved=false;
}


//______________________________________________________________________________


fastNLOInterpolBase::~fastNLOInterpolBase(void) {
}


//______________________________________________________________________________
fastNLOGrid::GridType fastNLOInterpolBase::TranslateGridType(string in){
   if ( in == "linear" ) return fastNLOGrid::kLinear;
   else if ( in == "loglog025" ) return fastNLOGrid::kLogLog025;
   else if ( in == "log10" ) return fastNLOGrid::kLog10;
   else if ( in == "sqrtlog10" ) return fastNLOGrid::kSqrtLog10;
   else {
      cout<<"fastNLOInterpolBase::TranslateGridTyp. Error. Cannot identify distance measure. in="<<in<<endl;
      exit(1);
   }
}


//______________________________________________________________________________



vector<pair<int,double> >* fastNLOInterpolBase::GetNodeValuesPtr(double x){
   cout<<"x="<<x<<",\tfLastVal="<<fLastVal<<endl;
   if ( x==fLastVal ) return &fNodes; // nothing todo. I know the nodes already.
   else {
      bool InRange = CheckX(x);
      if ( !InRange ) { // standard return, if there is no
      //         //fRetNodes = vector<pair<int,double> > ()
      //         //return fRetNodes;
      }
      fLastVal = x;
      CalcNodeValues(fNodes,x);
      return &fNodes;
   }
}

const vector<pair<int,double> >& fastNLOInterpolBase::GetNodeValues(double x){
   if ( fHgrid.empty() ) warn["GetNodeValues"]<<"There is no grid."<<endl;
   //cout<<"GetNodeValues for (x="<<x<<"), [fLastVal="<<fLastVal<<")"<<endl;
   if ( x==fLastVal ) return fNodes; // nothing todo. I know the nodes already.
   bool InRange = CheckX(x);
   if ( !InRange ) { // standard return, if there is no
      //fRetNodes = vector<pair<int,double> > ()
      //return fRetNodes;
   }
   fLastVal = x;
   CalcNodeValues(fNodes,x);
   return fNodes;
}

//______________________________________________________________________________


void fastNLOInterpolBase::RemoveLastNode(){
   debug["RemoveLastNode"]<<"Removing last node with highest value, but keep maximum value at fvalmax="<<fvalmax<<endl;
   fgrid.resize(fgrid.size()-1);
   fHgrid.resize(fHgrid.size()-1);
   fLastGridPointWasRemoved = true;
   debug["RemoveLastNode"]<<"last bin removed successful."<<endl;
}

//______________________________________________________________________________


void fastNLOInterpolBase::MakeGridsWithNNodesPerMagnitude(fastNLOGrid::GridType type, int nNodesPerMag){
   if ( fvalmin >= fvalmax ){
      warn["MakeGridsWithNNodesPerMagnitude"]<<"Minimum grid value is smaller/equal maximum value. min="<<fvalmin<<", max="<<fvalmax<<endl;
   }
   int nxtot = (int)(fabs(log10(fvalmax)-log10(fvalmin))*nNodesPerMag);
   if ( nxtot < nNodesPerMag ) nxtot = nNodesPerMag; // at least nxPerMagnitude points
   debug["MakeGridWithNNodesPerMagnitude"]<<"Create "<<nxtot<<" nodes (valmin="<<fvalmin<<",valmax="<<fvalmax<<")."<<endl;
   MakeGrids(type,nxtot+1); // plus 1. We have now a node at 1
}

void fastNLOInterpolBase::MakeGrids(fastNLOGrid::GridType type, int nNodes){
   // Generate the Hgrid first, and then the grid;
   // type: Type of distance measure
   // nNodes. Number of nodes. nNodes must be >= 1.
   // using fvalmin and fvalmin for the grid range

   debug["MakeGrid"]<<"Distance measure = "<<type<<endl;
   fdm = type;

   // check number of nodes
   if ( nNodes == -1 ) {
      error["MakeGrid"]<<"Minimum number of nodes not initialized. It seems that the (inherited) interpolation routine is missing."<<endl;
      exit(1);
   }
   if ( nNodes < fNMinNodes ) {
      error["MakeGrid"]<<"Number of nodes must be larger than "<<fNMinNodes<<" for this interpolation method."<<endl;
      exit(1);
   }
   //    else if ( nNodes == fNMinNodes ) {
   //     info["MakeGrid"]<<"This grid has only the minimum number of required nodes. nNodes="<<nNodes<<endl;
   //    }

   // check min and max values
   if ( fvalmin > fvalmax ){
      error["MakeGrid"]<<"Minimum grid value is smaller/equal maximum value. min="<<fvalmin<<", max="<<fvalmax<<endl;
   }

   MakeGrids(fvalmin,fvalmax,nNodes);

}

int fastNLOInterpolBase::FindLargestPossibleNode(double x){
   // --- find x position in range:  0 <= node1 < nnode1-1
   int node1 = fgrid.size()-2;  // --- initialize with largest possible value
   if ( fLastGridPointWasRemoved ) node1=fgrid.size()-1;
   if ( x < fgrid[0] ) {
      warn["FindLargestPossibleNode"]<<"Value is smaller than smallest node. Using first node. This may bias the result! x="<<x<<endl;
      return 0;
   }
   else if ( x==fgrid[0] ) {
      return 0;
   }
   if ( x > fgrid.back() ) {
      if ( !fLastGridPointWasRemoved )
         warn["FindLargestPossibleNode"]<<"Value is larger than largest node. Using last node. This may bias the result! x="<<x<<endl;
      else if ( x > fvalmax)
         warn["FindLargestPossibleNode"]<<"Value is larger than largest node and than largest grid value. Using last node. Interpolation kernel may lead unreasonable values! x="<<x<<endl;
      return node1;
   }
   //    if ( x > fvalmax ) {
   //       warn["FindLargestPossibleNode"]<<"Value is larger than maximum grid value. Using maximum grid value. Interpolation kernel may lead unreasonable values! x="<<x<<", valmax="<<fvalmax<<endl;
   //       return fvalmax;
   //    }

   //for( unsigned int iNode=1 ; iNode<fgrid.size()-2 ; iNode++ ){
   for( unsigned int iNode=1 ; iNode<fgrid.size() ; iNode++ ){
      if ( x <= fgrid[iNode]) {
         return iNode-1;
         //      node1=iNode-1;
         //      return node1; //done
      }
   }
   //error["FindLargestPossibleNode"]<<"Could not find largest node. x="<<x<<endl;
   return node1;
}


vector<double> fastNLOInterpolBase::MakeLinearGrid(double min, double max, int nNodes){
   vector<double> grid(nNodes);
   double dist = max-min;
   dist /= (nNodes - 1);
   for (int i=0 ; i<nNodes; i++ ){
      grid[i] = min+i*dist;
   }
   return grid;
}

void fastNLOInterpolBase::MakeGrids(double min, double max, int nNodes){
   // first make fhgrid and then make fgrid
   if ( isnan(min) || isnan(max) || nNodes <= 0 ) {
      error["MakeGrids"]<<"Cannot make unreasoanble grid! Requested: nNodes="<<nNodes<<", min="<<min<<", max="<<max<<". Exiting."<<endl;
      exit(1);
   }
   if ( fNMinNodes==1 && nNodes != 1) {
      warn["MakeGrids"]<<"Minimum number of nodes is 1. Number of nodes requested is "<<nNodes<<". Expecting, that this is a one-point grid. Therefore using only one node."<<endl;
      nNodes=1;
   }
   vector<double> hgrid(nNodes);
   double lo=min, hi=max;
   switch (fdm) {
   case fastNLOGrid::kLinear:
      lo = min;
      hi = max;
      break;
   case fastNLOGrid::kLogLog025:
      lo = Function_loglog025(min);
      hi = Function_loglog025(max);
      break;
   case fastNLOGrid::kLog10:
      lo = Function_log10(min);
      hi = Function_log10(max);
      break;
   case fastNLOGrid::kSqrtLog10:
      lo = Function_sqrtlog10(min);
      hi = Function_sqrtlog10(max);
      break;
   default:
      error["MakeGrid"]<<"Unknown grid type."<<endl;
   }
   if ( isnan(lo) || isnan(hi) ) {
      error["MakeGrids"]<<"Cannot convert min and max value to 'H'-space. min="<<min<<", H(min)="<<lo<<", max="<<max<<", H(max)="<<hi<<". Exiting."<<endl;
      exit(1);
   }
   double del = hi-lo;
   if ( nNodes==1 ) {
      hgrid[0] = lo+del/2.;
   }
   else if ( nNodes > 1 ) {
      for(int l=0;l<nNodes;l++){
         hgrid[l]   = lo +  double(l)/double(nNodes-1)*del;
         if ( isnan(hgrid[l]) ) {
            error["MakeGrids"]<<"Grid point could not be calculated. Hnode="<<hgrid[l]<<endl;
         }
      }
   }

   SetHGrid(hgrid);

   vector<double> g = MakeGridFromHGrid(fHgrid);
   SetGrid(g);
   //PrintGrid();
}

vector<double> fastNLOInterpolBase::MakeGridFromHGrid(vector<double> hg){
   if ( fHgrid.empty() ) { error["MakeGridFromHGrid"]<<"There is no HGrid."<<endl; exit(1);}
   vector<double> grid;
   switch (fdm) {
   case fastNLOGrid::kLinear:
      grid = hg;
      break;
   case fastNLOGrid::kLogLog025:
      grid = HGrid_loglog025_inv(hg);
      break;
   case fastNLOGrid::kLog10:
      grid = HGrid_log10_inv(hg);
      break;
   case fastNLOGrid::kSqrtLog10:
      grid = HGrid_sqrtlog10_inv(hg);
      break;
   default:
      error["MakeGridFromHGrid"]<<"Unknown grid type."<<endl;
   }
   return grid;
}


vector<double> fastNLOInterpolBase::HGrid_loglog025_inv(vector<double> grid){
   vector<double> ret = grid;
   for (unsigned int i=0 ; i<grid.size(); i++ ) {
      ret[i] = Function_loglog025_inv(grid[i]);
   }
   return ret;
}
vector<double> fastNLOInterpolBase::HGrid_log10_inv(vector<double> grid){
   vector<double> ret = grid;
   for (unsigned int i=0 ; i<grid.size(); i++ )
      ret[i] = Function_log10_inv(grid[i]);
   return ret;
}
vector<double> fastNLOInterpolBase::HGrid_sqrtlog10_inv(vector<double> grid){
   vector<double> ret = grid;
   for (unsigned int i=0 ; i<grid.size(); i++ )
      ret[i] = Function_sqrtlog10_inv(grid[i]);
   return ret;
}


void fastNLOInterpolBase::SetGrid(vector<double> grid){
   // todo. Add some checks here.
   fgrid = grid;
}
void fastNLOInterpolBase::SetHGrid(vector<double> hgrid){
   // todo. Add some checks here.
   fHgrid = hgrid;
}

bool fastNLOInterpolBase::CheckX(double& x) {
   bool sanity = false;
   if ( fgrid.size() == 1 ) {
      //x = fgrid[0];
      return true;
   }
   if ( x < fgrid[0] ) {
      double xin = x;
      x = fgrid[0];
      if ( x!=fLastVal )
         warn["CheckX"]<<"Value "<<xin<<" is smaller than smallest node (min="<<fgrid[0]<<"). Using this first node."<<endl;
   }
   else if ( x > fgrid.back() ) {
      double xin = x;
      if ( fLastGridPointWasRemoved ) {
         if ( x > fvalmax ) {
            x = fvalmax;
            if ( x!=fLastVal )
               warn["CheckX"]<<"Value "<<xin<<" is larger than largest grid value (max="<<fvalmax<<"). Using this value instead."<<endl;
         }
      }
      else {
         x = fgrid.back();
         if ( fabs(x-fgrid.back())>1.e-4 && x!=fLastVal)
               warn["CheckX"]<<"Value "<<xin<<" is larger than largest node (max="<<fgrid.back()<<"). Using this first node."<<endl;
      }
   }
   else
      sanity = true;

   return sanity;
}

double fastNLOInterpolBase::GetHx(double x ){
   switch(fdm) {
   case fastNLOGrid::kLinear:
      return x;
      break;
   case fastNLOGrid::kLogLog025:
      return Function_loglog025(x);
      break;
   case fastNLOGrid::kLog10:
      return Function_log10(x);
      break;
   case fastNLOGrid::kSqrtLog10:
      return  Function_sqrtlog10(x);
      break;
   default:
      error["GetHx"]<<"Unknown H-function measure."<<endl;
      return 0;
   }
}

double fastNLOInterpolBase::GetDelta(double x ){
   int node1 = FindLargestPossibleNode(x);
   double Hx = GetHx(x);
   if ( node1>= (int)fHgrid.size() ) {
      error["GetDelta"]<<"largest possible node is outside of grid."<<endl;
      exit(1);
   }
   // check if node next to largest node is also in grid
   if ( node1+1>= (int)fHgrid.size() ) {
      if ( fLastGridPointWasRemoved ) {
         //warn["GetDelta"]<<"Last grid point was removed. I assume this point was at x="<<fvalmax<<endl;
         double Hxone = GetHx(fvalmax);
         //cout<<"fHgrid[node1+1="<<fHgrid[node1+1]<<endl;
         double del = Hxone - fHgrid[node1];
         if ( del == 0 ){
            //warn["GetDelta"]<<"Distance between nodes is zero."<<endl;
            return 0;
         }
         return (Hx - fHgrid[node1]) / del;
      }
      else {
         error["GetDelta"]<<" node next to 'largest possible node' is outside of grid."<<endl;
         exit(1);
      }
   }
   double del = fHgrid[node1+1] - fHgrid[node1];
   if ( del == 0 ){
      //warn["GetDelta"]<<"Distance between nodes is zero."<<endl;
      return 0;
   }
   return (Hx - fHgrid[node1]) / del;
}


void fastNLOInterpolBase::PrintGrid() {
   warn["PrintGrid"]<<"\n ---------- printing grid -------------- " <<endl;
   warn>>"n grid nodes: " <<fgrid.size()<<endl;
   for (  unsigned int i = 0 ; i<fgrid.size() ; i++ ) {
      warn>>"i="<<i<<"\tnode="<<fgrid[i]<<endl;
   }
   warn["PrintGrid"]<<"\n ---------- printing Hgrid -------------- " <<endl;
   warn>>"n Hgrid nodes: " <<fHgrid.size()<<endl;
   for (  unsigned int i = 0 ; i<fHgrid.size() ; i++ ) {
      warn>>"i="<<i<<"\tnode="<<fHgrid[i]<<endl;
   }
   warn>>"----------------------------------------------------"<<endl;
}
