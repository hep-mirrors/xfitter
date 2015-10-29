#include <cmath>
#include <cstdlib>
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOPDFLinearCombinations.h"

using namespace std;
using namespace fastNLO;


fastNLOPDFLinearCombinations::fastNLOPDFLinearCombinations() {
   //: PrimalScream("fastNLOPDFLinearCombinations")  {
}


fastNLOPDFLinearCombinations::~fastNLOPDFLinearCombinations(){
}


//______________________________________________________________________________


vector<double > fastNLOPDFLinearCombinations::CalcPDFLinearCombination(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1, const vector<double>& pdfx2 , bool pdf2IsAntiParticle ) const {
   //
   // return vector of PDF linear combinations which can be used
   // to calculate the cross section.
   //
   // bool pdf2IsAntiParticle specifies, if pdfx2 still has to be 'inverted' or not
   //

   switch ( c->GetNPDF() ) {
   case 0:      // no PDF involved in process; e.g. e+e-
      return vector<double >();
      break;
   case 1:      // one PDF invovled in process: e.g. DIS
      return CalcPDFLCOneHadron(c,pdfx1);
      break;
   case 2:      // two PDFs involved in process: e.g. pp, ppbar
      if ( pdf2IsAntiParticle ) {
         vector<double> Antipdf2 = MakeAntiHadron(pdfx2);
         return CalcPDFLCTwoHadrons(c,pdfx1,Antipdf2);
      }
      else return CalcPDFLCTwoHadrons(c,pdfx1,pdfx2);
      break;
   default:
      //error["CalcPDFLinearCombination"]<<"Unknown number of PDFs involved in process. NPDF="<<c->GetNPDF()<<endl;
      say::error<<"[CalcPDFLinearCombination] Unknown number of PDFs involved in process. NPDF="<<c->GetNPDF()<<endl;
      exit(1);
      return vector<double >();
   }
}


//______________________________________________________________________________


vector<double > fastNLOPDFLinearCombinations::CalcPDFLCOneHadron(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 ) const {
   //
   // Calculate PDF linear combinations for DIS processes
   //

   // ---- check for ep-DIS ---- //
   bool IsONEPDF = ( c->GetNPDF() == 1 );
   bool IsDIS    = ( c->GetIPDFdef1() == 2 );
   bool IsNCDIS  = ( c->GetIPDFdef2() == 1 );
   bool IsProton = ( c->GetPDFPDG(0) == 2212 );
   if ( IsDIS && IsONEPDF && IsNCDIS && IsProton ) return CalcPDFDIS(c,pdfx1);
   // ---- unknown process ---- //
   else {
      //error["CalcPDFLCDIS"]<<"Could not identify process. Printing and exiting"<<endl;
      say::error<<"Error. Could not identify process. Printing and exiting"<<endl;
      c->Print();
      exit(1);
      return vector<double >();
   }
}



//______________________________________________________________________________



vector<double > fastNLOPDFLinearCombinations::CalcPDFLCTwoHadrons(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1, const vector<double>& pdfx2 ) const {
   //
   // Calculate PDF linear combinations for processes with two hadrons
   //
   // -----------------------------------------------
   // implement other processes e.g. like this:
   //    bool IsDrellYan = (...)   // check for the identifier(s) (most likely IPDFDef2)
   //    if ( IsDrellYan ) CalcPDFLinearCombDrellYan(c,pdfx1,pdfx2);
   // and implement a new function called CalcPDFLinearCombDrellYan(...)
   // -----------------------------------------------

   // ---- check for jet-productions  ---- //
   //bool IsTwoPDF = ( c->GetNPDF() == 2 );
   //bool IsTwoIdenticHadrons = (c->GetIPDFdef1() == 3  &&  c->GetPDFPDG(0) == fabs(c->GetPDFPDG(1)) );
   //bool IsHHJets = ( c->GetIPDFdef2() == 1 );
   //bool IsTTBar  = ( c->GetIPDFdef2() == 2 );

   if ( c->GetIPDFdef2()==0 ) // LiCos are stored in table
      return CalcPDFHHCFromTable(c,pdfx1,pdfx2);
   else if ( c->GetIPDFdef2()==1 &&  (c->GetIPDFdef3()==1 || c->GetIPDFdef3()==2))
      return CalcPDFHHC(c,pdfx1,pdfx2);
   else if ( c->GetIPDFdef2()==169 ) // default 169 PDF LiCos
      return CalcDefaultPDFLiCos(c,pdfx1,pdfx2);
   else if ( c->GetIPDFdef2()==121 ) // default 121 PDF LiCos
      return CalcDefaultPDFLiCos(c,pdfx1,pdfx2);
   else if ( c->GetIPDFdef2()==1 &&  c->GetIPDFdef3()==3 )
      return CalcPDFThreshold(c,pdfx1,pdfx2);
   else if ( c->GetIPDFdef2()==2 ) { // IPDFdef3==0 !
      return CalcPDFttbar(c,pdfx1,pdfx2);
   // else if (...)  //space for other processes
   // ---- (yet) unknown process ---- //
   }
   else {
      say::error<<"[CalcPDFLinearCombination] Could not identify process. Printing and exiting..."<<endl;
      say::error<<"PDFFlag1="<<c->GetIPDFdef1()<<endl;
      say::error<<"PDFFlag2="<<c->GetIPDFdef2()<<endl;
      say::error<<"PDFFlag3="<<c->GetIPDFdef3()<<endl;
      c->Print();
      exit(1);
      return vector<double >();
   }
}


//______________________________________________________________________________


vector<double > fastNLOPDFLinearCombinations::MakeAntiHadron(const vector<double>& xfx ) const {
   vector < double > xfxbar(13);
   for (unsigned int p = 0 ; p<13 ; p++) {
      xfxbar[p] = xfx[12-p];
   }
   return xfxbar;
}


//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcDefaultPDFLiCos(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   //! calculate default 121 or 169 PDF linear combinations
   //!
   //! Format is following for 121 subprocesses:
   //! iSubProc  parton1  parton2
   //! 0           -5       -5    #bbarbbar
   //! 1           -5       -4
   //! 2           -5       -3
   //! 3           -5       -2
   //! 4           -5       -1
   //! 5           -5        0
   //! 6           -5        1
   //! 7           -5        2
   //! 8           -5        3
   //! 9           -5        4
   //! 10          -5        5
   //! 11          -4       -5
   //! 12          -4       -4
   //! ...
   //! 60           0        0   #gg
   //! ...
   //! 120          5        5

   //! Format is following for 169 subprocesses:
   //! iSubProc  parton1  parton2
   //! 0        -6     -6
   //! 1        -6     -5
   //! 2        -6     -4
   //! 3        -6     -3
   //! 4        -6     -2
   //! 5        -6     -1
   //! 6        -6      0
   //! 7        -6      1
   //! 8        -6      2
   //! 9        -6      3
   //! 10       -6      4
   //! 11       -6      5
   //! 12       -6      6
   //! 13       -5     -6
   //! 14       -5     -5
   //! 15       -5     -4
   //! 16       -5     -3
   //! 17       -5     -2
   //! 18       -5     -1
   //! 19       -5      0
   //! 20       -5      1
   //!  ...
   //! 84        0      0   #gg
   //! ...
   //! 168       6      6
   int nSubproc = c->GetIPDFdef2();
   vector < double > pdflc(nSubproc);
   int istart = nSubproc==121 ? 1 : 0;
   int iend   = nSubproc==121 ? 12 : 13;
   int n=0;
   for ( int p1 = istart ; p1<iend ; p1++ ) {
      for ( int p2 = istart ; p2<iend ; p2++ ) {
         pdflc[n] = pdfx1[p1] * pdfx2[p2];
         n++;
      }
   }
   return pdflc;
}


//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFHHCFromTable(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   // calculate PDF linear combinations as stored in table
   if ( c->GetNSubproc() != c->GetIPDFdef3() || c->GetIPDFdef3() != (int)c->GetPDFCoeff().size()) {
      say::error["fastNLOPDFLinearCombinations::CalcPDFHHCFromTable"]
         <<"IPDFdef3 must be equal to NSubproc. (IPDFdef3="<<c->GetIPDFdef3()<<", NSubproc="<<c->GetNSubproc()<<"). Exiting."<<endl;
      exit(1);
   }
   const vector<vector<pair<int,int> > >& PDFCoeff = c->GetPDFCoeff();
   vector < double > pdflc(PDFCoeff.size());
   for ( unsigned int k = 0 ; k<PDFCoeff.size() ; k++ ){
      for ( unsigned int i = 0 ; i<PDFCoeff[k].size() ; i++ ){
         pdflc[k] += pdfx1[PDFCoeff[k][i].first+6] * pdfx2[PDFCoeff[k][i].second+6];
      }
   }
   return pdflc;
}


//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFDIS(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1) const {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //

   int NSubproc = c->GetNSubproc();

   vector < double > pdflc(3);
   pdflc[1] = pdfx1[6]; //gluon
   for (int l=0; l<13; l++) {
      double temp = (l==6 ? 0.0 : pdfx1[l]);
      if (!(l&1)) temp *= 4.;
      pdflc[0] += temp; // delta
   }
   pdflc[0] /= 9.;
   if (NSubproc>2) { // only from NLO
      for (int l=0; l<6; l++) {
         pdflc[2] += pdfx1[5-l] + pdfx1[l+7]; // sigma
      }
   }
   return pdflc;
}

//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFHHC(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //

   int NSubproc = c->GetNSubproc();

   double SumQ1  = 0;
   double SumQB1 = 0;
   double SumQ2  = 0;
   double SumQB2 = 0;
   vector <double> Q1(6);
   vector <double> QB1(6);
   vector <double> Q2(6);
   vector <double> QB2(6);
   for (int k = 0 ; k<6 ; k++) {
      Q1[k]  = pdfx1[k+7];  //! read 1st PDF at x1
      QB1[k] = pdfx1[5-k];
      SumQ1  += Q1[k];
      SumQB1 += QB1[k];
      Q2[k]  = pdfx2[k+7];//  ! read 2nd PDF at x2
      QB2[k] = pdfx2[5-k];
      SumQ2  += Q2[k];
      SumQB2 += QB2[k];
   }
   double G1     = pdfx1[6];
   double G2     = pdfx2[6];

   //   - compute S,A
   double S = 0;
   double A = 0;
   for (int k = 0 ; k<6 ; k++) {
      S += (Q1[k]*Q2[k]) + (QB1[k]*QB2[k]);
      A += (Q1[k]*QB2[k]) + (QB1[k]*Q2[k]);
   }

   //c   - compute seven combinations
   vector <double> H(7);
   H[0]  = G1*G2;
   H[1] = SumQ1*SumQ2 + SumQB1*SumQB2 - S;
   H[2] = S;
   H[3] = A;
   H[4] = SumQ1*SumQB2 + SumQB1*SumQ2 - A;
   H[5] = (SumQ1+SumQB1)*G2;
   H[6] = G1*(SumQ2+SumQB2);

   if (NSubproc == 6) {
      H[5] += H[6];
      H.resize(6);
   }

   // KR: For debugging purposes switch off subprocesses by uncommenting the following line(s)
   // H[0] = 0.;
   // H[1] = 0.;
   // H[2] = 0.;
   // H[3] = 0.;
   // H[4] = 0.;
   // H[5] = 0.;
   // H[6] = 0.;
   // KR DEBUG

   return H;

}


//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFThreshold(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //
   /*c---------------------------------------------------------------
     c       subprocess id : is = 1           ! q q' --> q  q'
     c       subprocess id : is = 2           ! q qb --> q' qb'
     c       subprocess id : is = 3           ! q qb --> q  qb
     c       subprocess id : is = 4           ! q q  --> q  q
     c       subprocess id : is = 5           ! q qb --> g  g
     c       subprocess id : is = 6           ! q g  --> q  g
     c       subprocess id : is = 7           ! g q  --> q  g
     c       subprocess id : is = 8           ! g g  --> q  qb
     c       subprocess id : is = 9           ! g g  --> g  g
     c       subprocess id : is = 10          ! q qb'--> q  qb'
     c---------------------------------------------------------------*/

   int NSubproc = c->GetNSubproc(); // MUST BE 10 here !

   double D1  = pdfx1[+1+6];
   double U1  = pdfx1[+2+6];
   double S1  = pdfx1[+3+6];
   double C1  = pdfx1[+4+6];
   double B1  = pdfx1[+5+6];
   //   double T1  = pdfx1[+6+6];
   double G1  = pdfx1[ 0+6];
   double D1b = pdfx1[-1+6];
   double U1b = pdfx1[-2+6];
   double S1b = pdfx1[-3+6];
   double C1b = pdfx1[-4+6];
   double B1b = pdfx1[-5+6];
   //   double T1b = pdfx1[-6+6];

   double D2  = pdfx2[+1+6];
   double U2  = pdfx2[+2+6];
   double S2  = pdfx2[+3+6];
   double C2  = pdfx2[+4+6];
   double B2  = pdfx2[+5+6];
   //   double T2  = pdfx2[+6+6];
   double G2  = pdfx2[ 0+6];
   double D2b = pdfx2[-1+6];
   double U2b = pdfx2[-2+6];
   double S2b = pdfx2[-3+6];
   double C2b = pdfx2[-4+6];
   double B2b = pdfx2[-5+6];
   //   double T2b = pdfx2[-6+6];

   vector <double> H(NSubproc);


   // 'q1q2 --> q1q2'
   //  q q'  + qb qb'
   // fqqp
   H[0] = D1*(U2+S2+C2+B2)
      +U1*(D2+S2+C2+B2)
      +S1*(D2+U2+C2+B2)
      +C1*(D2+U2+S2+B2)
      +B1*(D2+U2+S2+C2)
      +D1b*(U2b+S2b+C2b+B2b)
      +U1b*(D2b+S2b+C2b+B2b)
      +S1b*(D2b+U2b+C2b+B2b)
      +C1b*(D2b+U2b+S2b+B2b)
      +B1b*(D2b+U2b+S2b+C2b);

   //  'q1q1b --> q2q2b'
   // (xlqqb + xlqbq)
   //     q q + qb qb
   double fqq =  D1*D2+D1b*D2b
      +U1*U2+U1b*U2b
      +S1*S2+S1b*S2b
      +C1*C2+C1b*C2b
      +B1*B2+B1b*B2b;

   //     q qb + qb q
   double fqqb= D1*D2b
      +U1*U2b
      +S1*S2b
      +C1*C2b
      +B1*B2b;
   // c     qb q
   double fqbq= D1b*D2
      +U1b*U2
      +S1b*S2
      +C1b*C2
      +B1b*B2;

   ///   g*q pdfs are not working ?!?! no idea why!
   //        G1*G2 is working
   //        D1*D2 is working
   //         ...
   ////c     q g + qb g
   //    double fqg = (D1+U1+S1+C1+B1)*G2
   //       +(D1b+U1b+S1b+C1b+B1b)*G2;
   ////c     g q + g qb
   //    double fgq=  G1*(D2+U2+S2+C2+B2)
   //       +G1*(D2b+U2b+S2b+C2b+B2b);
   double fqg = 0;
   double fgq=  0;


   //c     q qb' + qb q'
   double fqqbp=D1*(U2b+S2b+C2b+B2b)
      +U1*(D2b+S2b+C2b+B2b)
      +S1*(D2b+U2b+C2b+B2b)
      +C1*(D2b+U2b+S2b+B2b)
      +B1*(D2b+U2b+S2b+C2b);

   double fqbpq=D1b*(U2+S2+C2+B2)
      +U1b*(D2+S2+C2+B2)
      +S1b*(D2+U2+C2+B2)
      +C1b*(D2+U2+S2+B2)
      +B1b*(D2+U2+S2+C2);


   H[1] = fqqb+fqbq;
   H[2] = fqqb+fqbq;
   H[3] = fqq;
   H[4] = fqqb+fqbq;

   H[5] = fqg; //geht nicht (s.o)
   H[6] = fgq; //geht nicht (s.o)

   H[7] = G1*G2;
   H[8] = G1*G2;
   H[9] = fqbpq+fqqbp;

   //       if (is .eq. 1) then
   //         sub_process = 'q1q2 --> q1q2'
   //    elseif (is .eq. 2) then
   //         sub_process = 'q1q1b --> q2q2b'
   //    elseif (is .eq. 3) then
   //         sub_process = 'q1q1b --> q1q1b'
   //    elseif (is .eq. 4) then
   //         sub_process = 'q1q1 --> q1q1'
   //    elseif (is .eq. 5) then
   //         sub_process = 'q1q1b --> gg'
   //    elseif (is .eq. 6) then
   //         sub_process = 'qg --> qg'
   //    elseif (is .eq. 7) then
   //         sub_process = 'gq --> qg'
   //    elseif (is .eq. 8) then
   //         sub_process = 'gg --> qqb'
   //    elseif (is .eq. 9) then
   //         sub_process = 'gg --> gg'
   //    elseif (is .eq. 10) then
   //         sub_process = 'q1q2b --> q1q2b'

   return H;

}


//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFttbar(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   //
   // Calculate pdf-lcs for ttbar cross sections in pp/ppbar
   // Used by NNLO generator of Marco Guzzi
   //
   // pdf[0] = gg
   // pdf[1] = qq
   //
   // DB 24.08.13
   //

   //int NSubproc = c->GetNSubproc();
   int IPDF3 = c->GetIPDFdef3();
   if ( IPDF3==0 ) { //NSubproc==2
      vector <double> pdflc(2);
      pdflc[0] += pdfx1[6]*pdfx2[6]; // gg
      // qq
      for (int k = 0 ; k<6 ; k++) {
         pdflc[1] += pdfx1[k]*pdfx2[12-k];
         pdflc[1] += pdfx1[12-k]*pdfx2[k];
      }
      return pdflc;
   }
   else if ( IPDF3==1 ) { //NSubproc==4
      // 0,1: gg, qq (inelastic)
      // 2,3: gg, qq (elastic)
      vector <double> pdflc(4);
      pdflc[0] += pdfx1[6]*pdfx2[6]; // gg
      pdflc[2]=pdflc[0];
      // qq
      for (int k = 0 ; k<6 ; k++) {
         pdflc[1] += pdfx1[k]*pdfx2[12-k];
         pdflc[1] += pdfx1[12-k]*pdfx2[k];
      }
      pdflc[3]=pdflc[1];
      return pdflc;
   }

   return vector<double>();

}
