// Daniel Britzger
// DESY, 08.08.2013
#ifndef __fastNLOTable__
#define __fastNLOTable__
#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include "speaker.h"

#include "fastNLOCoeffBase.h"
#include "fastNLOCoeffAddFix.h"
#include "fastNLOCoeffAddFlex.h"
#include "fastNLOCoeffData.h"
#include "fastNLOCoeffMult.h"
#include "fastNLOConstants.h"



class fastNLOTable {

 public:
   fastNLOTable();
   fastNLOTable(std::string filename);
   //   fastNLOTable(std::string filename, std::string verbosity = "INFO");
   virtual ~fastNLOTable();
   //   fastNLOTable(const fastNLOTable&, std::string verbosity = "INFO");
   fastNLOTable(const fastNLOTable&);

   virtual void ReadTable();
   virtual void WriteTable();
   virtual void WriteTable(std::string filename);
   bool IsCompatible(const fastNLOTable& other) const;
   bool IsCompatibleScenario(const fastNLOTable& other) const;
   bool IsCatenable(const fastNLOTable& other) const;
   bool IsCatenableScenario(const fastNLOTable& other) const;
   bool IsEquivalent(const fastNLOTable& other, double rtol) const;

   // --- function previously included in fastNLOBase
   // header
   void PrintHeader(int iprint) const;                                  //!< Print header variables (BlockA1) to screen
   bool IsCompatibleHeader(const fastNLOTable& other) const;             //!< Compare header with header of another table
   bool IsCatenableHeader(const fastNLOTable& other) const;              //!< Compare header with header of another table

   // getter/setters
   std::string GetFilename() const {return ffilename;}
   void   SetFilename(std::string name){ffilename=name;}

   //   int  GetItabversion() const {return Itabversion;}
   //   void SetItabversion(int version){Itabversion = version;}
   int  GetITabVersionRead() const {return ITabVersionRead;}
   int  GetITabVersionWrite() const {return ITabVersionWrite;}
   void SetITabVersionRead(int version){ITabVersionRead = version;}
   void SetITabVersionWrite(int version){ITabVersionWrite = version;}

   std::string GetScenName() const {return ScenName;}
   void   SetScenName(std::string name){ScenName = name;}

   int  GetNmult() const;
   int  GetNcontrib() const;
   int  GetNdata() const;

   int  GetOutputPrecision() const {return fPrecision;}
   void SetOutputPrecision(int precision) {fPrecision = precision;}

   /// _____________________________________________________________________________________________
   /// Getters for binning structure
   /// _____________________________________________________________________________________________

   /// Get dimensionality of calculation: single-, double-, or triple-differential
   unsigned int GetNumDiffBin() const {return NDim;}

   /// Getters/setters for linear array of observable bins "ObsBin" running from 0->(NObsBin-1)

   /// Return lower bin bound for obs. bin iObs in dim. iDim
   double GetObsBinLoBound(unsigned int iObs, unsigned int iDim) const;
   /// Return upper bin bound for obs. bin iObs in dim. iDim
   double GetObsBinUpBound(unsigned int iObs, unsigned int iDim) const;
   /// Return std::vector of lower bin bounds in dim. iDim for all obs. bins
   std::vector < double > GetObsBinsLoBounds(unsigned int iDim) const;
   /// Set std::vector of lower bin bounds in dim. iDim for all obs. bins
   void SetObsBinsLoBounds(unsigned int iDim, std::vector < double > v);
   /// Return std::vector of upper bin bounds in dim. iDim for all obs. bins
   std::vector < double > GetObsBinsUpBounds(unsigned int iDim) const;
   /// Set std::vector of upper bin bounds in dim. iDim for all obs. bins
   void SetObsBinsUpBounds(unsigned int iDim, std::vector < double > v);
   /// Return minimum value of all lower bin bounds for dim. iDim
   double GetObsBinsLoBoundsMin(unsigned int iDim) const;
   /// Return maximum value of all upper bin bounds for dim. iDim
   double GetObsBinsUpBoundsMax(unsigned int iDim) const;
   /// Return std::vector of pairs with lower and upper bin bounds in dim. iDim for all obs. bins
   std::vector < std::pair < double, double > > GetObsBinsBounds(unsigned int iDim) const;
   /// Return observable bin no. for std::vector of values obs0=var0,obs1=var1,...; -1 if outside range
   int GetObsBinNumber(const std::vector < double >& vobs) const ;
   /// Return observable bin no. for obs0=var0 in 1D binning; -1 if outside range
   int GetObsBinNumber(double var0) const ;
   /// Return observable bin no. for obs0=var0,obs1=var1 in 2D binning; -1 if outside range
   int GetObsBinNumber(double var0, double var1) const ;
   /// Return observable bin no. for obs0=var0,obs1=var1,obs2=var2 in 3D binning; -1 if outside range
   int GetObsBinNumber(double var0, double var1, double var2) const ;

   /// Getters for multidimensional binning, here called Dim<I>Bins

   /// Return std::vector of pairs with unique bin bounds of 1st dim.
   std::vector < std::pair < double, double > > GetDim0BinBounds() const;
   /// Return std::vector of pairs with unique bin bounds of 2nd dim. for 'iDim0Bin' of 1st dim.
   std::vector < std::pair < double, double > > GetDim1BinBounds(unsigned int iDim0Bin) const;
   /// Return std::vector of pairs with unique bin bounds of 3rd dim. for 'iDim0Bin' and 'iDim1Bin' of 1st two dim.
   std::vector < std::pair < double, double > > GetDim2BinBounds(unsigned int iDim0Bin, unsigned int iDim1Bin) const;
   /// Return std::vector of pairs with lower and upper bin bounds for all dimensions for a given obs. bin
   std::vector < std::pair < double, double > > GetObsBinDimBounds(unsigned int iObs) const;
   /// Return pair with lower and upper bin bounds for given obs. bin and dim. iDim
   std::pair < double, double > GetObsBinDimBounds(unsigned int iObs, unsigned int iDim) const;
   /// Return bin no. in 1st dim. for obs. bin iObs
   unsigned int GetIDim0Bin(unsigned int iObs) const;
   /// Return bin no. in 2nd dim. for obs. bin iObs
   unsigned int GetIDim1Bin(unsigned int iObs) const;
   /// Return bin no. in 3rd dim. for obs. bin iObs
   unsigned int GetIDim2Bin(unsigned int iObs) const;
   /// Return no. of bins in 1st dimension
   unsigned int GetNDim0Bins() const;
   /// Return no. of bins in 2nd dimension for given bin in 1st dim.
   unsigned int GetNDim1Bins(unsigned int iDim0Bin) const;
   /// Return no. of bins in 3rd dimension for given bins in 1st and 2nd dim.
   unsigned int GetNDim2Bins(unsigned int iDim0Bin, unsigned int iDim1Bin) const;
   /// Return bin no. in 1st dim. for obs0=var0; -1 if outside range
   int GetODim0Bin(double var0) const;
   /// Return bin no. in 2nd dim. for obs0=var0,obs1=var1; -1 if outside range
   int GetODim1Bin(double var0, double var1) const;
   /// Return bin no. in 3rd dim. for obs0=var0,obs1=var1,obs2=var2; -1 if outside range
   int GetODim2Bin(double var0, double var1, double var2) const;
   // DO NOT USE! DOES NOT WORK!
   //   unsigned int GetIDimBin(unsigned int iObs, unsigned int iDim) const;
   //   std::vector < std::pair < double, double > > GetBinBoundaries(int iDim0Bin, int iDim1Bin = -1, int iDim2Bin = -1);

   /// ___________________________________________________________________________________________________
   /// Some more info getters with respect to observable dimensions
   /// ___________________________________________________________________________________________________

   /// Get if dimension is 'truly differential' or bin-integrated (divided by bin width or not)
   int GetIDiffBin(int bin) const {return IDiffBin[bin];}
   /// Get BinSize for bin = BinSizeDim1 < * BinSizeDim2 >
   double GetBinSize(int bin) const {return BinSize[bin];};
   /// Get vector of dimensions labels
   std::vector < std::string > GetDimLabels() const {return DimLabel;};
   /// Get dimension label for dimension iDim
   std::string GetDimLabel(int iDim) const {return DimLabel[iDim];};

   /// ___________________________________________________________________________________________________
   /// Some info getters with respect to normalization
   /// ___________________________________________________________________________________________________

   /// Get normalization flag:
   ///    def=0        -> no norm.
   ///     1, 2, 3,... -> normalize to slice in NDim of same table
   ///    -1,-2,-3,... -> normalize to slice in NDim of other table
   int GetINormFlag() const {return INormFlag;};
   /// Get normalization logical (def=false)
   bool IsNorm() const {return INormFlag == 0 ? false : true;}
   /// Get filename of normalization table for INormFlag<0
   std::string GetDenomTable() const {return DenomTable;}

   /// ___________________________________________________________________________________________________
   /// Some info getters & setters for table modifications
   /// ___________________________________________________________________________________________________

   /// get/set scenario description
   std::vector <std::string> GetScDescr() const { return ScDescript; }
   void SetScDescr(std::vector <std::string> ScDescr);

   /// get/set cross section units of published results (pb = 12, fb = 15, ...)
   int GetIpublunits() const {return Ipublunits;}
   void SetIpublunits(int unit){Ipublunits = unit;}

   /// get/set center-of-mass energy in units of GeV
   double GetEcms() const {return Ecms;}
   void SetEcms(double E) {Ecms = E;}

   /// get/set power of alpha_s for LO process
   int GetLoOrder() const {return ILOord;}
   void SetLoOrder(int LOOrd);

   /// get/set no. of observable bins
   unsigned int GetNObsBin() const {return NObsBin;}
   void SetNObsBin(int NObs);

   /// get/set Bin vector
   std::vector < std::vector <std::pair<double,double> > > GetBins() const {return Bin;};
   void SetBins(std::vector < std::vector <std::pair<double,double> > >);

   /// get/set BinSize vector
   std::vector < double > GetBinSize() const {return BinSize;};
   void SetBinSize(std::vector < double >);

   // Erase observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to NObsBin
   void EraseBinFromTable(unsigned int iObsIdx);
   template<typename T> void EraseBin(std::vector<T>& v, unsigned int idx);

   // Multiply observable bin; iObsIdx is the C++ array index to be multiplied and
   // not the observable bin no. running from 1 to NObsBin
   void MultiplyBinInTable(unsigned int iObsIdx, double fact);
   void MultiplyBinSize(unsigned int iObsIdx, double fact);
   template<typename T> void MultiplyBin(std::vector<T>& v, unsigned int idx, double fact);
   void MultiplyBinBorders(unsigned int iDim, double fact);

   void CatBinToTable(const fastNLOTable& other, unsigned int iObsIdx, unsigned int table_count);
   void CatBin(const fastNLOTable& other, unsigned int iObsIdx, unsigned int table_count);

   /// ???
   /// Get Rivet ID of analysis
   std::string GetRivetId() const;
   /// Get cross section from analysis description
   std::string GetXSDescr() const;
   void SetDimLabel(std::string label, unsigned int iDim, bool IsDiff = true);
   void SetNumDiffBin(int iDiff) {NDim=iDiff; DimLabel.resize(NDim); IDiffBin.resize(NDim);}

   //void Cat(const fastNLOCoeffBase& other);


   /// ___________________________________________________________________________________________________
   /// Info print out functionality
   /// ___________________________________________________________________________________________________

   ///  Print basic info about fastNLO table and its contributions
   void PrintTableInfo(const int iprint = 0) const; // DEPRECATED, use PrintContributionSummary instead
   void PrintContributionSummary(int iprint) const;
   ///  Print (technical) constants of fastNLO table (use iprint) for level of details.
   void PrintFastNLOTableConstants(const int iprint = 0) const;  // DEPRECATED, use PrintContributionSummary instead
   void PrintScenario(int iprint) const;
   virtual void Print(int iprint) const;

   /// ___________________________________________________________________________________________________
   /// Other useful functions
   /// ___________________________________________________________________________________________________
   void MergeTable(const fastNLOTable& rhs, fastNLO::EMerge option=fastNLO::kMerge ); //!< 'merge'
   void MergeTables(const std::vector<fastNLOTable*>& tables, fastNLO::EMerge option=fastNLO::kMerge, double cutRMS=0 ); //!< 'merge' (also supports 'median' and 'mean')
   void AddTable(const fastNLOTable& rhs, fastNLO::EMerge option=fastNLO::kMerge); //!< 'merge'
   void SetUserWeights(double wgt); //!< Set user weights for subsequent mergeing wgt
   void SetUserWeights(std::vector<double> wgtsObs); //!< Set user weights for subsequent mergeing wgt[obs]
   void SetUserWeights(std::vector<std::vector<double> > wgtsBinProc); //!< Set user weights for subsequent mergeing wgt[proc][obs]

   /// Handle coefficient tables
   //int WriteCoeffTable(int no);
   //int WriteCoeffTable(int no, ofstream* outstream);
   //int WriteCoeffTableDividebyN(int no);
   void DeleteAllCoeffTable();
   //int CreateCoeffBase(int no);
   int CreateCoeffTable(int no, fastNLOCoeffBase *newcoeff);
   void CatenateTable(const fastNLOTable& other);
   fastNLOCoeffBase* GetCoeffTable(int no) const;
   /// Returns pointer to data table if available, else returns NULL pointer
   fastNLOCoeffData* GetDataTable() const;
   /// Returns pointer to reference table if available, else returns NULL pointer
   fastNLOCoeffAddBase* GetReferenceTable(fastNLO::ESMOrder eOrder) const;

private:
   bool cmp(const double x1, const double x2) const;
   bool cmp(const std::vector<double>& x1, const std::vector<double >& x2) const;
   bool cmp(const std::vector<std::vector<double> >& x1,const std::vector<std::vector<double > >& x2) const;
   bool cmp(const std::vector<std::vector<std::pair<double,double> > >&  x1,const std::vector<std::vector<std::pair<double,double> > >& x2) const;

protected:
   // --- functions previously included in fastNLOBase
   void PrintWelcomeMessage();                         //!< Say hello to fastNLO user
   std::ostream* OpenFileWrite(bool compress=false); //!< open std::ofstream for writing tables to ffilename
   std::istream* OpenFileRead();                     //!< open std::ifstream for reading table
   //std::ofstream *OpenFileRewrite();
   void WriteHeader(std::ostream& table);          //!< write (or cout) hader using std::ostream
   int ReadHeader(std::istream& table);           //!< read header of table (BlockA1)
   void CloseFileWrite(std::ostream& table);
   void CloseFileRead(std::istream& table);
   //void CloseStream();

   std::string ffilename;
   int fPrecision;
   //   int Itabversion;
   int ITabVersionRead;
   int ITabVersionWrite = fastNLO::tabversion;
   std::string ScenName;

   PrimalScream logger;
   static bool fWelcomeOnce;
   // ---- fastNLOBase end


   void WriteScenario(std::ostream& table);
   void ReadScenario(std::istream& table);
   void ReadCoeffTables(std::istream& table, int nCoeff);
   fastNLOCoeffBase* ReadRestOfCoeffTable(const fastNLOCoeffBase& cB, std::istream& table, int ITabVersionRead);

   std::vector < fastNLOCoeffBase* > fCoeff;
   //fastNLOCoeffData* fData;

   double Ecms;
   int ILOord;
   int Ipublunits;
   std::vector <std::string> ScDescript;

   // Unsigned int
   unsigned int NObsBin;
   unsigned int NDim;

   std::vector <std::string> DimLabel;
   std::vector <int> IDiffBin;
   // Every bin has a lower and upper bin boundary and belongs to a 'dimension'. In a point-wise differential measurement, the upper bin boundary is equal to the lower one.
   std::vector < std::vector <std::pair<double,double> > > Bin;
   std::vector <double> BinSize;

   // Contributions for normalization
   int INormFlag;
   std::string DenomTable;
   std::vector <int> IDivLoPointer;
   std::vector <int> IDivUpPointer;

};
#endif
