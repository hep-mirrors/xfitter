#ifndef __fastNLOCoeffAddBase__
#define __fastNLOCoeffAddBase__


#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffAddBase : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddBase(int NObsBin);
   fastNLOCoeffAddBase(const fastNLOCoeffBase& base);
   virtual ~fastNLOCoeffAddBase() {;}
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   void Read(istream& table);
   virtual void Write(ostream& table);
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Print() const;
   virtual void Clear();                                                        //!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients();                                        //!< Set number of events to 1 and normalize coefficients accordingly.

   int GetIRef() const {return IRef;}
   double GetNevt() const { return Nevt; }
   double GetNevt(int NObsBin, int NSubproc) const {
      if (Nevt > 0) return Nevt;
      else {cout<<"Todo. Preparation for v2.3."<<endl; return Nevt;}
   }
   void SetNevt(double nevt) { Nevt = nevt;}                                    //!< Set number of events
   int GetNxmax(int Obsbin) const ;
   int GetXIndex(int Obsbin,int x1bin,int x2bin =0) const ;
   int GetNSubproc() const { return NSubproc;}
   int GetIScaleDep() const { return IScaleDep;}
   int GetNPDF() const {return NPDFPDG.size();}
   int GetPDFPDG(int iPDF) const {return NPDFPDG[iPDF];}
   int GetNPDFDim() const {return NPDFDim;}
   int GetIPDFdef1() const { return IPDFdef1; }
   int GetIPDFdef2() const { return IPDFdef2; }
   int GetIPDFdef3() const { return IPDFdef3; }
   int GetNpow() const {return Npow;}
   int GetNScales() const {return NScales;}
   int GetNScaleDim() const {return NScaleDim;}
   //vector<string > GetScaleDescript(int iScale=0) const { return ScaleDescript[iScale]; };
   string GetScaleDescription(int iScale=0) const { return ScaleDescript[0][iScale]; };         // getter for scale description of scale iScale
   vector<vector<string > > GetScaleDescr() const { return ScaleDescript; }
   int GetNxtot1(int iBin) const { return XNode1[iBin].size(); }
   int GetNxtot2(int iBin) const { return XNode2[iBin].size(); }

   double GetXNode1(int iObsBin, int iNode) const { return XNode1[iObsBin][iNode]; }
   double GetXNode2(int iObsBin, int iNode) const { return XNode2[iObsBin][iNode]; }

   bool IsReference() const {return IRef>0;};
   bool IsCompatible(const fastNLOCoeffAddBase& other) const;

   const vector<vector<pair<int,int> > >& GetPDFCoeff() const { return fPDFCoeff;}

protected:
   fastNLOCoeffAddBase();
   void ReadCoeffAddBase(istream& table);
   int GetScaledimfromvar(int scalevar) const;

   int IRef;
   int IScaleDep;
   double Nevt;
   int Npow;
   vector < int > NPDFPDG;
   int NPDFDim;
   vector < int > NFFPDG;
   int NFFDim;
   int NSubproc;
   int IPDFdef1;
   int IPDFdef2;
   int IPDFdef3;
   vector<vector<pair<int,int> > > fPDFCoeff;                                                   //! fPDFCoeff[iSubProc][iPartonPair][pair]
   // Missing: linear PDF combinations for IPDFdef1=0
   vector < double > Hxlim1;
   v2d XNode1;
   vector < double > Hxlim2;
   v2d XNode2;
   vector < int > Nztot;
   vector < double > Hzlim;
   v2d ZNode;
   int NScales;
   int NScaleDim;
   vector < int > Iscale;                                                                       // not used
   vector < vector < string > > ScaleDescript;

};

#endif
