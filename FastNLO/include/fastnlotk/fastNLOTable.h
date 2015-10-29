// Daniel Britzger
// DESY, 08.08.2013
#ifndef __fastNLOTable__
#define __fastNLOTable__
#include <cmath>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "fastNLOBase.h"
#include "fastNLOCoeffBase.h"
#include "fastNLOCoeffAddFix.h"
#include "fastNLOCoeffAddFlex.h"
#include "fastNLOCoeffData.h"
#include "fastNLOCoeffMult.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOTable : public fastNLOBase {

 public:
   fastNLOTable();
   fastNLOTable(string filename);
   ~fastNLOTable();
   fastNLOTable(const fastNLOTable&);

   virtual void ReadTable();
   virtual void WriteTable();
   virtual void WriteTable(string filename);
   bool IsCompatible(const fastNLOTable& other) const;

   // ___________________________________________________________________________________________________
   // Getters for binning structure
   // ___________________________________________________________________________________________________
   // Getters for linear array of "ObsBin"'s running from 0->(NObsBin-1))
   // Returns no. of observable bins
   unsigned int GetNObsBin() const {return NObsBin;}
   // Return lower bin bound for obs. bin iObs in dim. iDim
   double GetObsBinLoBound(unsigned int iObs, unsigned int iDim) const;
   // Return upper bin bound for obs. bin iObs in dim. iDim
   double GetObsBinUpBound(unsigned int iObs, unsigned int iDim) const;
   // Return vector of lower bin bounds in dim. iDim for all obs. bins
   vector < double > GetObsBinsLoBounds(unsigned int iDim) const;
   // Return vector of upper bin bounds in dim. iDim for all obs. bins
   vector < double > GetObsBinsUpBounds(unsigned int iDim) const;
   // Return minimum value of all lower bin bounds for dim. iDim
   double GetObsBinsLoBoundsMin(unsigned int iDim) const;
   // Return maximum value of all upper bin bounds for dim. iDim
   double GetObsBinsUpBoundsMax(unsigned int iDim) const;
   // Return vector of pairs with lower and upper bin bounds in dim. iDim for all obs. bins
   vector < pair < double, double > > GetObsBinsBounds(unsigned int iDim) const;
   // Return observable bin no. for vector of values obs0=var0,obs1=var1,...; -1 if outside range
   int GetObsBinNumber(const vector < double >& vobs) const ;
   // Return observable bin no. for obs0=var0 in 1D binning; -1 if outside range
   int GetObsBinNumber(double var0) const ;
   // Return observable bin no. for obs0=var0,obs1=var1 in 2D binning; -1 if outside range
   int GetObsBinNumber(double var0, double var1) const ;
   // Return observable bin no. for obs0=var0,obs1=var1,obs2=var2 in 3D binning; -1 if outside range
   int GetObsBinNumber(double var0, double var1, double var2) const ;

   // Getters for multidimensional binning, here called Dim<I>Bins
   // Return vector of pairs with unique bin bounds of 1st dim.
   vector < pair < double, double > > GetDim0BinBounds() const;
   // Return vector of pairs with unique bin bounds of 2nd dim. for 'iDim0Bin' of 1st dim.
   vector < pair < double, double > > GetDim1BinBounds(unsigned int iDim0Bin) const;
   // Return vector of pairs with unique bin bounds of 3rd dim. for 'iDim0Bin' and 'iDim1Bin' of 1st two dim.
   vector < pair < double, double > > GetDim2BinBounds(unsigned int iDim0Bin, unsigned int iDim1Bin) const;
   // Return vector of pairs with lower and upper bin bounds for all dimensions for a given obs. bin
   vector < pair < double, double > > GetObsBinDimBounds(unsigned int iObs) const;
   // Return pair with lower and upper bin bounds for given obs. bin and dim. iDim
   pair < double, double > GetObsBinDimBounds(unsigned int iObs, unsigned int iDim) const;
   // Return bin no. in 1st dim. for obs. bin iObs
   unsigned int GetIDim0Bin(unsigned int iObs) const;
   // Return bin no. in 2nd dim. for obs. bin iObs
   unsigned int GetIDim1Bin(unsigned int iObs) const;
   // Return bin no. in 3rd dim. for obs. bin iObs
   unsigned int GetIDim2Bin(unsigned int iObs) const;
   // Return no. of bins in 1st dimension
   unsigned int GetNDim0Bins() const;
   // Return no. of bins in 2nd dimension for given bin in 1st dim.
   unsigned int GetNDim1Bins(unsigned int iDim0Bin) const;
   // Return no. of bins in 3rd dimension for given bins in 1st and 2nd dim.
   unsigned int GetNDim2Bins(unsigned int iDim0Bin, unsigned int iDim1Bin) const;
   // Return bin no. in 1st dim. for obs0=var0; -1 if outside range
   int GetODim0Bin(double var0) const;
   // Return bin no. in 2nd dim. for obs0=var0,obs1=var1; -1 if outside range
   int GetODim1Bin(double var0, double var1) const;
   // Return bin no. in 3rd dim. for obs0=var0,obs1=var1,obs2=var2; -1 if outside range
   int GetODim2Bin(double var0, double var1, double var2) const;
   // DO NOT USE! DOES NOT WORK YET!
   unsigned int GetIDimBin(unsigned int iObs, unsigned int iDim) const;
   vector < pair < double, double > > GetBinBoundaries(int iDim0Bin, int iDim1Bin = -1, int iDim2Bin = -1);


   // ___________________________________________________________________________________________________
   // Some other info getters
   // ___________________________________________________________________________________________________
   vector < double > GetBinSize() const {return BinSize;};                                      //!< Get Binsize = BinSizeDim1 < * BinSizeDim2 >
   double GetBinSize(int bin) const {return BinSize[bin];};                                     //!< Get Binsize = BinSizeDim1 < * BinSizeDim2 >
   void SetNumDiffBin(int iDiff ) { NDim=iDiff; DimLabel.resize(NDim); IDiffBin.resize(NDim);}  //!< Set dimension of calculation. (Single-differential, double-differential, etc...)
   unsigned int GetNumDiffBin() const { return NDim; }                                                   //!< Get dimension of calculation. (Single-differential, double-differential, etc...)

   int GetIDiffBin(int bin) const { return IDiffBin[bin]; }                                     //!< Get if dimension is 'truly differential' or bin-integrated (divided by bin-width or not)

   void SetDimLabel( string label, unsigned int iDim , bool IsDiff = true );
   string GetDimLabel( int iDim  ) const {return DimLabel[iDim];};                              //!< Get label (name) of observable in dimension iDim
   vector<string > GetDimLabels() const {return DimLabel;};                                     //!< Get label (name) of all observables

   void SetIpublunits(int unit){Ipublunits = unit;}
   int GetIpublunits() const {return Ipublunits;}

   double GetEcms() const {return Ecms;}
   void SetEcms(double E) {Ecms = E;}

   int GetLoOrder() const {return ILOord;}
   void SetLoOrder(int LOOrd);

   bool IsNorm() const { return INormFlag == 0 ? false : true;}
   string GetDenomTable() const {return DenomTable;}

   string GetRivetId() const;


   // ___________________________________________________________________________________________________
   // Info print out functionality
   // ___________________________________________________________________________________________________
   //  Print basic info about fastNLO table and its contributions
   void PrintTableInfo(const int iprint = 0) const;
   //  Print (technical) constants of fastNLO table (use iprint) for level of details.
   void PrintFastNLOTableConstants(const int iprint = 0) const;
   void PrintScenario() const;
   virtual void Print() const;


   // ___________________________________________________________________________________________________
   // Other useful functions
   // ___________________________________________________________________________________________________
   // handle coefficient tables
   //int WriteCoeffTable(int no);
   //int WriteCoeffTable(int no,ofstream* outstream );
   //int WriteCoeffTableDividebyN(int no);
   void DeleteAllCoeffTable();
   //int CreateCoeffBase(int no);
   int CreateCoeffTable(int no,fastNLOCoeffBase *newcoeff);
   void AddTable(const fastNLOTable& rhs);
   fastNLOCoeffBase* GetCoeffTable(int no) const;
   fastNLOCoeffData* GetDataTable() const;                                                      //!< returns pointer to data table if available, else returns NULL pointer
   fastNLOCoeffAddBase* GetReferenceTable(ESMOrder eOrder) const;                               //!< returns pointer to reference table if available, else returns NULL pointer


private:
   bool cmp(const double x1, const double x2) const;
   bool cmp(const vector<double>& x1, const vector<double >& x2) const;
   bool cmp(const vector<vector<double> >& x1,const vector<vector<double > >& x2) const;
   bool cmp(const vector<vector<pair<double,double> > >&  x1,const vector<vector<pair<double,double> > >& x2) const;


protected:
   void WriteScenario(ostream& table);
   void ReadScenario(istream& table);
   void ReadCoeffTables(istream& table);
   fastNLOCoeffBase* ReadRestOfCoeffTable(const fastNLOCoeffBase& cB, istream& table);

   vector < fastNLOCoeffBase* > fCoeff;
   //fastNLOCoeffData* fData;

   double Ecms;
   int ILOord;
   int Ipublunits;
   vector <string> ScDescript;

   // CKR: Changed to unsigned int
   unsigned int NObsBin;
   unsigned int NDim;

   vector <string> DimLabel;
   vector <int> IDiffBin;
   vector < vector <pair<double,double> > > Bin; // every bin has a lower and upper bin boundary and belongs to a 'dimension'. If a truely differential measurment, then upper bin boundary is equal lower one
   vector <double> BinSize;

   // contributions for normalization
   int INormFlag;
   string DenomTable;
   vector <int> IDivLoPointer;
   vector <int> IDivUpPointer;

};
#endif
