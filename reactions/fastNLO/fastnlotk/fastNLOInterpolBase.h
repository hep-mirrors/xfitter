// Author: Daniel Britzger
// DESY, 28/06/2013

#ifndef __fastNLOInterpolBase__
#define __fastNLOInterpolBase__


#include "speaker.h"
#include <string>
#include <vector>
#include <cmath>
#include <utility>

namespace fastNLOGrid {
   enum GridType {
      kLinear           = 0,            // linear grid
      kLog10            = 1,            // log10 grid
      kLogLog025        = 2,            // loglog grid (only valid for mu-grids)
      kLogLog           = 3,            // loglog grid (only valid for mu-grids)
      kSqrtLog10        = 4,            // sqrt(logarithmic) grid  (only valid for x-grids)
      k3rdrtLog10       = 5,            // log(x)^(1/3) (3rd root) (only valid for x-grids)
      k4thrtLog10       = 6             // log(x)^(1/4) (4th root) (only valid for x-grids)
   };
}


class fastNLOInterpolBase : public PrimalScream {

public:

   fastNLOInterpolBase(double min, double max, fastNLOGrid::GridType type, int nMinNodes);
   fastNLOInterpolBase(double density, fastNLOGrid::GridType type, int nMinNodes);
   virtual ~fastNLOInterpolBase(void);

   const std::vector<std::pair<int,double> >& GetNodeValues(double val);

   void MakeGrids(int nNodes,double ReduceXmin=0);
   void MakeGridsWithNNodesPerMagnitude(int nNodes,double ReduceXmin=0);
   void RemoveLastNode();

   void PrintGrid();
   const std::vector<double>& GetGrid() const { return fgrid;}
   const std::vector<double>* GetGridPtr() const { return &fgrid;}
   const std::vector<double>& GetHGrid() const { return fHgrid;}
   double GetDelta(double);
   bool CheckX(double&);

   static fastNLOGrid::GridType TranslateGridType(std::string in);
   std::vector<double> fgrid;

protected:

   void SetGrid(std::vector<double> grid);
   void SetHGrid(std::vector<double> grid);
   void MakeGrids(double min, double max, int nNodes);
   std::vector<double> MakeGridFromHGrid(std::vector<double> g);
   std::vector<double> MakeLinearGrid(double min, double max, int nNodes);

   //virtual std::vector<std::pair<int,double> > CalcNodeValues(double val) = 0;
   virtual void CalcNodeValues(std::vector<std::pair<int,double> >& nodes, double val) = 0;

   int FindLargestPossibleNode(double, bool);

   inline double Function_loglog025( double mu ){
      // function H(mu) = log(log( mu / 0.25 ))
      return log(log(mu/0.25));}
   inline double Function_loglog025_inv( double mu ){
      // inverse of function H(mu) = log(log( mu / 0.25 ))
      return 0.25*exp(exp(mu));}
   inline double Function_loglog( double mu ){
      // function H(mu) = log(log( mu ))
      return log(log(mu));}
   inline double Function_loglog_inv( double mu ){
      // inverse of function H(mu) = log(log( mu ))
      return exp(exp(mu));}
   inline double Function_x( double mu ){
      // function H(mu) = x
      return mu;}
   inline double Function_x_inv( double mu ){
      // inverse of function H(mu) = x;
      return mu;}
   inline double Function_log10( double x ){return log10(x);}
   inline double Function_log10_inv( double x ){return pow(10,x);}
   inline double Function_sqrtlog10( double x ){return -sqrt(-log10(x));}
   inline double Function_sqrtlog10_inv( double x ){return pow(10,-pow(x,2));}
   inline double Function_3rdrtlog10( double mu ){return -pow(fabs(log10(mu)),1./3.);}
   inline double Function_3rdrtlog10_inv( double mu ){return pow(10,-pow(fabs(mu),3));}
   inline double Function_4thrtlog10( double mu ){return -pow(fabs(log10(mu)),0.25);}
   inline double Function_4thrtlog10_inv( double mu ){return pow(10,-pow(mu,4));}

   std::vector<double> HGrid_loglog025_inv(std::vector<double> grid);
   std::vector<double> HGrid_loglog_inv(std::vector<double> grid);
   std::vector<double> HGrid_log10_inv(std::vector<double> grid);
   std::vector<double> HGrid_sqrtlog10_inv(std::vector<double> grid);
   std::vector<double> HGrid_4thrtlog10_inv(std::vector<double> grid);
   std::vector<double> HGrid_3rdrtlog10_inv(std::vector<double> grid);
   int GetNMod() const {return fnmod;}
   double GetHx(double);

protected:
   std::vector<std::pair<int,double> > fNodes;

   int fNMinNodes;
   double fvalmin;
   double fvalmax;
   double fLastVal[5] = {M_PI,M_PI,M_PI,M_PI,M_PI};
   bool fLastGridPointWasRemoved; // odd boolean to agree with original code;
   fastNLOGrid::GridType fdm; // distance measure
   std::vector<double> fHgrid;
   int fnmod ; // variable for final nodes. Has to be filled by inherited algorithm
   bool fExtendLow = false;

};

#endif
