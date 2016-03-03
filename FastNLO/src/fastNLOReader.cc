// Author: Daniel Britzger
// DESY, 08/2013
// Updated: Klaus Rabbertz
// KIT, 01/2014

//______________________________________________________________________________
/**
    --- fastNLO user: Hello!
        If you use fastNLO for the first time, please read through the
        documentation and comments carefully in order to calculate
        a reasonable cross section.
        All comments that start with '--- fastNLO user:' are intended as a
        short documentation for various options, that can be changed by you.

        In fastNLO version 2, there are two different types of tables.
        Although internally they are implemented slightly differently, both are called
        v2 for their larger flexiblity compared to version 1.4.
        The simpler ones, v2.0, are extended versions of this previous format
        v1.4 from which a conversion into v2.0 is possible, but without profiting
        of the improvements, of course.
        The second type of tables, v2.1, are called 'flexible-scale' tables
        which store the matrix elements in a scale independent way and
        also the scale variables are stored more generally.
        These tables give you the possibility to change in addition to the renormalization
        also the factorization scale by arbitrary factors and have the possiblity to
        change the formula according to which the scale is derived.

        Please chek (see point 2 below), which type of table you are using and
        then refer to the comments and functions suitable for this fastNLO table.


         0.  This Introduction
         1.  Instantiation of fastNLO classes
              a. basic FastNLOUser class
              b. FastNLOLHAPDF
              c. fastNLOAlphas
              d. FastNLOCRunDec
         2.  Print table information
         3.  Calculate cross sections
         4.  Modify PDF settings for LHAPDF-interfaces
         5.  Change alphas values
         6.  Units of the calculation
         7.  Contributions and order of calculation
         8.  Scale settings
         9.  Flexible-scale concept
        10.  Access cross section and k-factor
        11.  Print-out cross section
        12.  Verbosity level
        13.  FastNLO for jet-production diffractive DIS
              a. Introduction
              b. The FastNLOUser.h class
              c. Calculate diffractive cross sections.
              d. Diffractive DIS example code
        14.  Example code
        15.  Example analysis


    1a.
    ------- Initialize table for fastNLOReader ------- //
    --- fastNLO user:
        In addition to a fastNLO table two additional ingredients are required:
        - the PDF set and
        - the alpha_s evolution to be used
        These can be freely defined by the user by making an instance of your class
        that derives from the FastNLOReader class and passing the name of the
        fastNLO table as an argument, e.g.:
           FastNLOUser* fnlo = new FastNLOUser( tablename );

        To facilitate using fastNLOReader a number of predefined user classes
        of FastNLOUser exist, interfacing to
        LHAPDF (PDF and alpha_s, see M. Whalley, D. Bourilkov, R. Group, hep-ph/0508110),
        GRV Alphas (default alpha_s evolution used in fastNLO for crosschecks,
                    based on M. Glueck, E. Reya, A. Vogt, Eur.Phys.J.C5:461-470,1998, hep-ph/9806404;
                    PDF from LHAPDF),
        CRunDec (alpha_s evolution up to 4 loops, see B. Schmidt, M. Steinhauser,
                 Comput.Phys.Commun. 183 (2012) 1845-1848, arXiv:1201.6149;
                 PDF from LHAPDF).

        Their use is explained in the following.


    1b.
        Initialize with PDF from LHAPDF and corresponding alphas value and
        evolution for this PDF set. A change of the alpha_s value is only
        possible through the choice of the PDF file/set and member, e.g. CT10as.LHgrid
            FastNLOLHAPDF fnlolhapdf( tablename , PDFFile , PDFMember );

        Print information from LHAPDF
            fnlolhapdf.PrintPDFInformation();
            int npdf = fnlolhapdf.GetNPDFMembers();
            int imaxpdf = fnlolhapdf.GetNPDFMaxMember(); // imaxpdf = npdf - 1

        ( Please note that because of a feature in gfortran the output via your LHAPDF
          installation may be asynchronous to the C++ output. Usually, the gfortran
          output comes at the end after all C++ output, but this depends on your actual system.
          You can try to set the environment variable GFORTRAN_UNBUFFERED_ALL to yes
          in your shell to get it synchronized. Keep your fingers crossed. )


    1c.
        Initialize with PDF from LHAPDF and GRV alphas evolution (default)
            fastNLOAlphas fnlo( tablename , PDFFile , PDFMember );

        Change the alpha_s value through
            fnlo.SetAlphasMz(0.1179);
        Change values of the alpha_s evolution code through:
            Alphas::SetNf(5);
            Alphas::SetMz(91.70);
        For all options see Alphas.h


    1d.
        Initialize with PDF from LHAPDF and RunDec alpha_s evolution.
            FastNLOCRunDec fnlo( tablename , PDFFile , PDFMember );
        Change the alpha_s value for all instances, by:
            fnlo.SetAlphasMz(0.1179);
        Change values of the alpha_s evolution code through:
            fnlo.SetNf(5);
            fnlo.SetNloop(4);
            fnlo.SetMz(91.70);

        (Note: CTEQ6M:   M_Z = 91.70,   alpha_s(M_Z) = 0.1179;
               PDG 2012: M_Z = 91.1876, alpha_s(M_Z) = 0.1184)



    2.
    ---- Table information ---- //
    --- fastNLO user: For a comprehensive insight into the fastNLO variables
        you can use:
                fnlo.PrintFastNLOTableConstants(0);



    3.
    ---- (Re-)calculate cross sections ---- //
    --- fastNLO user: Before you can access the fastNLO computed
        cross sections, you always have to call CalcCrossSection()!
        So, before accessing the cross sections, please call:
                fnlo.CalcCrossSection();


    4.
    ------- Select another PDF set and member ------- //
    --- fastNLO user: You can select another PDF set and member here.
        With LHAPDF, you can set the PDF set and member using e.g.:
              fnlo.SetLHAPDFFilename( string PDFFile );
              fnlo.SetLHAPDFMember( int PDFMember );



    5.
    ------- Changing the alpha_s(M_Z) value and/or evolution ------- //
    --- fastNLO user:
        The alpha_s evolution is provided by the code of the chosen
        interface, e.g. GRV alpha_s for the fnlo instance here.
        The value of alpha_s(M_Z) can be changed from its default PDG 2012 values
        like this:

               fnlo.SetAlphasMz(0.1179);

        (Note: CTEQ6M:   M_Z = 91.70,   alpha_s(M_Z) = 0.1179;
               PDG 2012: M_Z = 91.1876, alpha_s(M_Z) = 0.1184)

        To use a different alpha_s evolution code one has to interface it.
        Here, for example, we use the above-mentioned CRunDec code:

     FastNLOCRunDec fnlocrundec( tablename , PDFFile , 0 );
     fnlocrundec.SetMz(91.1876);
     fnlocrundec.SetAlphasMz(0.1184);
     fnlocrundec.CalcCrossSection();


    6.
    ------- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ------- //
    --- fastNLO user: You can choose the units in which you want
        to access (or print) your cross-section results.
        There are two possibilites:
          - The default option is 'publication units', i.e. divided by
            bin widths if done so in the relevant publication
               fnlo.SetUnits(fastNLO::kPublicationUnits);
          - The other option is 'absolute' units in barn, but still in
            the same magnitude as in the publication (e.g. pb, fb, nb, etc.)

          fnlo.SetUnits(kAbsoluteUnits); // in namespace fastNLO
        or
          fnlo.SetUnits(kPublicationUnits); // in namespace fastNLO


    7.
    ------- Set the calculation order (if available) ------- //
    --- fastNLO user: Each fastNLO table comes typically with
        various contributions.
        Currently, five different types of contributions have been tested.
        Three can be combined to give a scale, PDF and alpha_s dependent
        cross-section, one is a fixed multiplicative correction and, at last,
        also data points with uncertainties might be included in a table.
        For calculating a cross section, by default only the LO & NLO contributions
        are used. However, each contribution can be swiched on or off separately.
        Please make sure to avoid combinations that do not make sense,
        e.g. 2-loop threshold corrections with LO pQCD.

        For switching a contribution on/off, its type must be known:
          - kFixedOrder                  -> Fixed order calculation (in alpha_s)
          - kThresholdCorrection         -> Threshold corrections
          - kElectroWeakCorrection       -> Electroweak corrections (not derived yet)
          - kNonPerturbativeCorrections  -> Non-perturbative corrections|Hadronisation corrections
        plus one must know the 'Id' of this contribution, which can be printed e.g.
        by calling
           fnlo.PrintTableInfo();

        To switch a contribution on/off please use:
               bool SetOn = fnlo.SetContributionON( contrib, Id, on/off )
        and in particular for switching on check on the return value SetOn that it actually worked.
        Here, 'contrib' is not the contribution number, but the type
        as given above: kFixedOrder, ...
        Within each type the contributions are counted separately starting with Id=0.
        The total number of contributions then counts all contributions of all types.


    8.
    ------- Selecting the scale treatment ------- //
    --- fastNLO user: The simplest way to modify the predefined renormalization and
        factorization scales is to provide a scale factor by which the default scale
        is multiplied. These factors must be positive and not too small (> 1.e-6).
        Otherwise they can in principal (within reason) be set arbitrarily for
        flexible-scale tables. For the normal v2 tables the choice of factors for the
        factorization scale is limited to some fixed values, usually 0.5, 1.0, and 2.0
        plus sometimes also 0.25, see the respective table information.
        Note: If threshold corrections are available and switched on for evaluation,
        the scale factors for the renormalization and factorization scale must be identical.

        The function call to set the scale factors is:
            bool SetScales = fnlo.SetScaleFactorsMuRMuF(xmur, xmuf);
        where xmur, xmuf are the scale factors. Check the return value in order to verify
        that the selected scale factors could actually be activated.

        The return value of this function call is boolean and returns false, if the
        the requested scale factors can not be chosen. In this case, the last legal
        values remain unchanged.


    9.
    ----- Additional possibilities for scales in 'flexible-scale' tables (v2.1) ----- //
        First check, if your table is a flexible-scale table or not
             bool IsFlex = fnlo.GetIsFlexibleScaleTable()
        You can choose a function to define how
        to compute the renormalization and factorization scale.
        Each 'flexible-scale' table comes with two variables that can be used
        for calculating the scales. They are called scale1 and scale2 and
        at least one needs to have a dimension in "GeV".
        DIS tables have typically stored scale1 = Q and scale2 = pt, while
        hadron-hadron tables might have for example scale1 = pt and scale2 = y.
        Other settings are imaginable. Please check, which obervables exactly
        are stored as scale variables!

        There are two possibilities, how you can define your scale now:

          - use predefined functions using e.g.
               fnlo.SetMuRFunctionalForm(fastNLO::EScaleFunctionalForm);
            for changing the calculation of the renormalizatoin scale.
            Please refer to FastNLOReader.h for all options of EScaleFunctionalForm.

          - or you can pass a function pointer to FastNLOReader using
               fnlo.SetExternalFuncForMuR( double (*Func)(double,double) );
            to pass any function using scale1 and scale2 to fastNLO.

        WARNING: Some choice had to be made for the default settings. Please think
        carefully about the choice of the scales ...
        Default setting for DIS tables:
          - mu_r:  kQuadraticMean      -> mu_r = sqrt( (Q^2 + scale2^2)/2. ) // because scale1=Q!
          - mu_f:  kScale1             -> mu_f = Q
        Default setting for pp and ppbar tables:
          - mu_r:  kScale1             -> mu_r = scale1
          - mu_f:  kScale1             -> mu_f = scale1

        Valid calls are e.g.:
        fnlo.SetMuRFunctionalForm(fastNLO::kScale1);        // set function how to calculate mu_r from scale1 and scale2
        fnlo.SetMuFFunctionalForm(fastNLO::kScale1);        // set function how to calculate mu_f from scale1 and scale2
        fnlo.SetMuRFunctionalForm(fastNLO::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
        fnlo.SetMuFFunctionalForm(fastNLO::kScale1);        // set function how to calculate mu_f from scale1 and scale2
        fnlo.SetExternalFuncForMuR( &Function_Mu );         // set external function to calculate mu_r from scale1 and scale2
        fnlo.SetMuRFunctionalForm(fastNLO::kExpProd2);      // set function how to calculate mu_f from scale1 and scale2
        fnlo.SetMuFFunctionalForm(fastNLO::kExpProd2);      // set function how to calculate mu_f from scale1 and scale2

    INFO: All above-mentioned scale changing functions automatically perform a refilling of the
          fastNLO internal PDF cache. To switch it off you can use a boolean, like:
          fnlo.SetMuFFunctionalForm(fastNLO::kScale1 , false );


    10.
    ---- Access cross sections ---- //
    --- fastNLO user: To access the cross section from fastNLO
        you should use:
              vector < double > xs = fnlo.GetCrossSection();
        If you want to have a pointer to an array of numbers you might use
              vector < double > xs = fnlo.GetCrossSection();
              double* cs = &xs[0];

        Further you can access the "k-factor", which is calculated with all
        'contributions' that are switched on (e.g. non-perturbative corrections)
        against the LO fixed-order contribution.
        Remark:
             - the proverbial k-factor is NLO vs. LO
             - 1-loop threshold corrections are vs. LO
             - 2-loop threshold corrections are vs. NLO
             - non-perturbative corrections usually are vs. NLO

              vector < double > kFactors = fnlo.GetKFactors();


    11.
    ---- Printing ---- //
    --- fastNLO user: For an easy overview of your cross section calculation
        you might use the following print methods:
                fnlo.PrintCrossSections();

        Or print it like the Fortran reader code:
                fnlo.PrintCrossSectionsDefault();


    12.
    ------- Set fastNLOReader verbosity ------- //
    --- fastNLO user:
        The following line sets the verbosity level of fastNLOReader
        Six different levels are implemented, the default is INFO:
        DEBUG, MANUAL, INFO, WARNING, ERROR, SILENT
            SetGlobalVerbosity(WARNING);
        Alternatively, a specific verbosity level can be set
        to any instance:
            fnlo.SetVerbosity(level);


    13.
    ------- FastNLO for jets in diffractive DIS ------- //
    13a.
     FastNLO is also applicable to jets in diffractive DIS.
     The calculation of jet cross sections in diffractive
     DIS is performed by adapting the slicing method,
     where the xpom integration is performed during the evaluation
     of the fastNLO table. The differential cross section
     in xpom is calcualted by a rescaling of the center-of-mass
     energy of the incident hadron.
     The boundaries of the integration interval are automatically
     smoothed out.
     More details on the applied method can be found on the
     website, i.e.
     http://fastnlo.hepforge.org/docs/talks/20120912_fastNLOv2_DBritzger_DiffractiveFastNLO.pdf

    13b.
    --- fastNLO user:
     In order to calculate diffractive DIS processes, the user
     has to provide a diffractive PDF, as well as an alpha_s
     evolution code. Both pieces have to be implemented in the
     FastNLODiffUser.h file, where the functions
        double FastNLODiffUser::EvolveAlphas(double Q)
        bool FastNLODiffUser::InitPDF()
        vector<double> FastNLODiffUser::GetDiffXFX(double xpom, double zpom, double muf)
     have to be implemented in a reasonable way.
     Some examples and more help on this, can provide the authors.
     The implementation of the alpha_s evolution code can also be
     adapted e.g. from fastNLOAlphas.h or FastNLOCRunDec.h.

    13c.
     The calculation of diffractive cross sections performs
     an integration of xpom. This is done by a simple Riemann integration.
     Four possibilities to define the slicing are implemented.
     1. Use a logarithmic xpom slicing
        Set the number of slices, the xpom_min and xpom_max range, e.g.:
          fnlodiff->SetXPomLogSlicing( 12, pow(10.,-2.3), pow(10.,-1) );
     2. Use a linear xpom slicing
          fnlodiff->SetXPomLinSlicing( 12, 0.0, 0.1 );
     3. Use an exponential xpom slicing
     4. Set your individual xpom slicing. This basically also allows
        to implement a MC integration.
            nStep:     number of slices
            xpom[]:    central value of each slice
            dxpom[]:   width of each slice
        fnlodiff->SetXPomSlicing(int nStep, double* xpom, double* dxpom);

     To calculate and access the cross sections use:
           vector<double> xs = fnlodiff->GetDiffCrossSection();

     If you want to calculate cross sections as fucntion of xpom,
     you have to calculate each xpom bin by setting the 'xpomslicing', and
     summing all bins by yourself.
     WARNING:
     In this case, one always have to call SetUnits(fastNLO::kAbsoluteUnits) !

     Tipp 1: Some brief studies showed, that already with ca. 10 slices, the
     cross section converges sufficiently fast. The linear slicing is
     preferred over the logarithmic slicing.
     Tipp 2:
     Choosing Q2 (or pT) as factorization scale increases the speed significantly.

    13d.
     In the following example code the class description FastNLODiffUser may be
     replaced by a specific interface class to a diffractive PDF (see 13b). This class
     has to be added to the include statements above, e.g.:
        #include "fastnlo/FastNLODiffUser.h"
     Some example code could look like (uncomment the following lines,
     comment out the other examples under 14. and 15., and recompile):

       // ---- Example code for jet-cross sections in diffractive DIS ---- //
       //  we setup an instance of the FastNLODiffUser class
       FastNLODiffUser fnlodiff( tablename );

       //  If you want to receive your cross section in
       //   pb/GeV or in pb. Here we choose pb/GeV
       fnlodiff.SetUnits(fastNLO::kPublicationUnits);

       // Set the xpom integration interval and method
       fnlodiff.SetXPomLinSlicing( 12, 0.0, 0.1 );

       // Optional:
       // make your scale definition (see above)
       fnlodiff.SetMuFFunctionalForm(kQuadraticSum);
       fnlodiff.SetMuRFunctionalForm(kQuadraticSum);
       fnlodiff.SetScaleFactorsMuRMuF(1.0,1.0);

       // calculate and access the cross section
       vector<double>  xs = fnlodiff.GetDiffCrossSection();
       // Print it
       fnlodiff.PrintCrossSections();
       // ------------------------------------------------------------------ //


    14.
    ---- Example of a cross section calculation with some nice standardized output
    fastNLOLHAPDF fnlo(tablename,"cteq6m.LHpdf",0);
    fnlo.CalcCrossSection();
    fnlo.PrintCrossSectionsDefault();

*/
//______________________________________________________________________________

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <string>
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOReader.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOCoeffAddFix.h"
#include "fastnlotk/fastNLOCoeffAddFlex.h"
#ifdef FNLO_HOPPET
#include "fastnlotk/HoppetInterface.h"
#endif

using namespace std;
using namespace fastNLO;
using namespace say;

//______________________________________________________________________________
fastNLOReader::fastNLOReader() : fastNLOTable() {
   logger.SetClassName("fastNLOReader");
   fUnits               = fastNLO::kPublicationUnits;
   fMuRFunc             = fastNLO::kScale1;
   fMuFFunc             = fastNLO::kScale1;
   fPDFSuccess          = false;
   fAlphasCached        = 0.;
   fPDFCached           = 0.;
   fUseHoppet            = false;
}


//______________________________________________________________________________

fastNLOReader::fastNLOReader(string filename) : fastNLOTable(filename) {
   //SetGlobalVerbosity(DEBUG); // Temporary for debugging
   logger.SetClassName("fastNLOReader");
   debug["fastNLOReader"]<<"New fastNLOReader reading filename="<<filename<<endl;
   //fCoeffData           = NULL;
   //    Coeff_LO_Ref         = NULL;
   //    Coeff_NLO_Ref        = NULL;
   fUnits               = fastNLO::kPublicationUnits;
   fMuRFunc             = fastNLO::kScale1;
   fMuFFunc             = fastNLO::kScale1;
   fPDFSuccess          = false;
   fAlphasCached        = 0.;
   fPDFCached           = 0.;
   fUseHoppet            = false;
   SetFilename(filename);
}


//______________________________________________________________________________
fastNLOReader::~fastNLOReader(void) {
}



//______________________________________________________________________________
fastNLOReader::fastNLOReader(const fastNLOReader& other) :
   fastNLOTable(other),
   ffilename(other.ffilename), fScalevar(other.fScalevar), fScaleFacMuR(other.fScaleFacMuR),
   fUnits(other.fUnits), fPDFSuccess(other.fPDFSuccess), fPDFCached(other.fPDFCached),
   fAlphasCached(other.fAlphasCached), Fct_MuR(other.Fct_MuR), Fct_MuF(other.Fct_MuF),
   bUseSMCalc(other.bUseSMCalc),
   XSection_LO(other.XSection_LO), XSection(other.XSection), kFactor(other.kFactor),
   QScale_LO(other.QScale_LO), QScale(other.QScale), XSectionRef(other.XSectionRef),
   XSectionRefMixed(other.XSectionRefMixed), XSectionRef_s1(other.XSectionRef_s1),
   XSectionRef_s2(other.XSectionRef_s2)
{
   //! Copy constructor
   OrderCoefficients(); // initialize pointers to fCoeff's
}



//______________________________________________________________________________
void fastNLOReader::SetFilename(string filename) {
   debug["SetFilename"]<<"New filename="<<filename<<endl;
   ffilename    = filename;
   OrderCoefficients();
   SetCoefficientUsageDefault();
   InitScalevariation();
}


//______________________________________________________________________________
void fastNLOReader::OrderCoefficients() {
   debug["OrderCoefficients"]<<endl;

   // Initialize Coeff's
   fastNLOCoeffBase* Coeff_LO   = NULL;
   fastNLOCoeffBase* Coeff_NLO  = NULL;
   fastNLOCoeffBase* Coeff_NNLO = NULL;
   fastNLOCoeffBase* Coeff_THC1 = NULL;
   fastNLOCoeffBase* Coeff_THC2 = NULL;
   fastNLOCoeffBase* Coeff_NPC1 = NULL;

   // run over all coefficient tables, identify and sort contributions.
   for (unsigned int i= 0; i<fCoeff.size() ; i++ ){
      fastNLOCoeffBase* c = GetCoeffTable(i);
      // give contribution a reasonable name
      //char nbuf[400];
      //       sprintf(nbuf,"Coeff. %s %s %s",
      //              _ContrName[c->GetIContrFlag1()-1].c_str(),_OrdName[c->GetIContrFlag1()-1][c->GetIContrFlag2()-1].c_str(),_fNSDep[c->GetNScaleDep()].c_str());
      //       c->SetName(nbuf);

      // data
      if ( fastNLOCoeffData::CheckCoeffConstants(c,true) ) {
         debug["OrderCoefficients"]<<"Found data table."<<endl;
      }
      // additive contributions
      else if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
         // Reference table
         if ( ((fastNLOCoeffAddBase*)c)->IsReference() ) {
            debug["OrderCoefficients"]<<"Found reference table."<<endl;
         }
         // Additive fixed order (perturbative) contribution
         else if ( c->GetIContrFlag1() == 1 ) {
            if ( c->IsLO() )            Coeff_LO  = c;
            else if ( c->IsNLO() )      Coeff_NLO = c;
            else if ( c->IsNNLO() )     Coeff_NNLO = c;
         }
         // Threshold corrections
         else if ( c->GetIContrFlag1() == 2 ) {
            if ( c->GetIContrFlag2() == 1 )             Coeff_THC1 = c;
            else if ( c->GetIContrFlag2() == 2 )        Coeff_THC2 = c;
            else {
               error["OrderCoefficients"]<<"Threshold correction implemented only up to 2-loops, exiting!\n";
               exit(1);
            }
         }
      }
      // multiplicative corrections
      else if ( fastNLOCoeffMult::CheckCoeffConstants(c,true) ) {
         // Non-perturbative corrections
         if ( c->GetIContrFlag1()==4 )  Coeff_NPC1 = c;
         else {
            error["ReadTable"]<<"Further multiplicative corrections not yet implemented, stopped!\n";
            exit(1);
         }
      }
   }

   // Delete and re-initialize lists for BlockB's
   const int defsize = 10;
   BBlocksSMCalc.clear();
   BBlocksSMCalc.resize(defsize);

   // Assign non-perturbative corrections, switch off by default
   if (Coeff_NPC1) {
      BBlocksSMCalc[kNonPerturbativeCorrection].push_back(Coeff_NPC1);
   }

   // Assign threshold corrections, switch off by default
   if (Coeff_THC1) {
      BBlocksSMCalc[kThresholdCorrection].push_back(Coeff_THC1);
   }
   if (Coeff_THC2) {
      BBlocksSMCalc[kThresholdCorrection].push_back(Coeff_THC2);
   }

   // Assign fixed order calculations (LO must be [0]), switch on by default
   if (Coeff_LO)  {
      BBlocksSMCalc[kFixedOrder].push_back(Coeff_LO);
   } else {
      error["OrderCoefficients"]<<"Could not find any LO Calculation. Exiting!"<<endl;
      exit(1);
   }
   if (Coeff_NLO) {
      BBlocksSMCalc[kFixedOrder].push_back(Coeff_NLO);
   } else {
      info["OrderCoefficients"]<<"Could not find any NLO calculation."<<endl;
   }
   if (Coeff_NNLO) {
      BBlocksSMCalc[kFixedOrder].push_back(Coeff_NNLO);
   } else {
      info["OrderCoefficients"]<<"Could not find any NNLO calculation."<<endl;
   }

   //int iprint = 2;
   //PrintFastNLOTableConstants(iprint);
}


//______________________________________________________________________________
void fastNLOReader::SetCoefficientUsageDefault() {
   //! Switch on LO and NLO contribution.
   //! Deactivate all other contributions
   bUseSMCalc.clear();
   bUseSMCalc.resize(BBlocksSMCalc.size());

   // Switch all off
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
         bUseSMCalc[j].push_back(false);
      }
   }
   // Activate LO and NLO
   bUseSMCalc[kFixedOrder][kLeading] = true; //LO
   bUseSMCalc[kFixedOrder][kNextToLeading] = true;//NLO
}


//______________________________________________________________________________
void fastNLOReader::InitScalevariation() {
   //! Initialize to scale factors of (MuR,MuF) = (1,1)
   debug["InitScalevariation"]<<"Try to initialize scale factors MuR and MuF to (1,1)."<<endl;

   if (!GetIsFlexibleScaleTable()) {
      bool SetScales = SetScaleFactorsMuRMuF(1.,1.);
      if (!SetScales){
         error["InitScalevariation"]<<"Could not find scale variation with scale factor 1.0, stopped!"<<endl;
         exit(1);
      }
   // if (!GetIsFlexibleScaleTable()) {
   //    if ( B_NLO() ) {
   //       fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();
   //       for (int iscls=0; iscls< GetNScaleVariations(); iscls++) {
   //          const double muFac = cNLO->GetScaleFactor(iscls);
   //          if (fabs(muFac-1.0) < 1.e-7) {
   //             SetScaleVariation(iscls,true);
   //             break;
   //          }
   //       }
   //       if (fScalevar == -1) {
   //          error["InitScalevariation"]<<"Could not found scale variation with scale factor 1.0. Exiting.\n";
   //          exit(1);
   //       }
   //    } else {
   //       // no NLO table -> no scale dependent contribrutions
   //    }
   } else {
      // this is a MuVar table. You can vary mu_f and mu_r independently by any factor
      // and you can choose the functional form of mu_f and mu_r as functions of
      // scale1 and scale1 (called partly scaleQ2 and scalePt).
      fScaleFacMuR = 1.;
      fScaleFacMuF = 1.;
      fastNLOCoeffAddFlex* cNLO = (fastNLOCoeffAddFlex*)B_NLO();
      if ( !cNLO ) cNLO = (fastNLOCoeffAddFlex*)B_LO();
      // KR: Commented out since anyway always false (warning from clang compiler)!
      //     cNLO->GetScaleDescr()[0].size() is of type unsigned int
      //     The functional form for scale variations is set below.
      // if (cNLO->GetScaleDescr()[0].size() <0) { // ???
      //    warn["InitScalevariation"]<<"No scaledescription available.\n";
      //    SetFunctionalForm(kScale1 , kMuR);
      //    SetFunctionalForm(kScale1 , kMuF);
      //    return;
      // }

      // ---- DIS ---- //
      if (cNLO->GetNPDF() == 1) {
         SetFunctionalForm(kQuadraticMean , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
      }
      // ---- HHC --- //
      else if (cNLO->GetNPDF() == 2) {
         SetFunctionalForm(kScale1 , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
      } else {
         error<<"Unknown process.\n";
         exit(1);
      }
   }
}


//______________________________________________________________________________
bool fastNLOReader::SetScaleVariation(int scalevar) {
   //! ------------------------------------------------
   //!   NEVER call this setter directly, only via
   //!   the method SetScaleFactorsMuRMuF!
   //!
   //!   Set the scale variation table to correspond
   //!   to the selected MuF factor if possible.
   //!   Usually, v2.0 tables are stored for multiple
   //!   MuF settings like factors of 0.5, 1.0 and 2.0
   //!   times the nominal scale, e.g.
   //!     scalevar -> scalefactor
   //!        '0'   ->   1.0
   //!        '1'   ->   0.5
   //!        '2'   ->   2.0
   //!   If tables for multiple MuF factors are present,
   //!   then they MUST correspond to exactly the same
   //!   factors in the SAME order for all such contrbutions,
   //!   e.g. NLO plus 2-loop threshold corrections!
   //!
   //!   This method returns true if the chosen
   //!   'scalevar' table exists for all selected
   //!   contributions with extra scale tables.
   //! ------------------------------------------------
   debug["SetScaleVariation"]<<"Setting to scalevar table "<<scalevar<<endl;

   if (GetIsFlexibleScaleTable()) {
      warn["SetScaleVariation"]<<"WARNING! This is a flexible-scale table. MuF scale variation tables are not necessary!"<<endl;
      warn["SetScaleVariation"]<<"You should not have called this method for the active table. Nothing changed!"<<endl;
      return false;
   }

   // Check for maximal scale variation of all active SM calcs
   int scalevarmax = GetNScaleVariations();
   if (scalevar >= scalevarmax) {
      error["SetScaleVariation"]<<"This table has only "<<scalevarmax<<" scale variation(s) stored for all active contributions!"<<endl;
      error["SetScaleVariation"]<<"You wanted to access the non-existing number "<<scalevar<<", stopped!"<<endl;
      exit(1);
   }

   fScalevar     = scalevar;

   // TBD
   // The following is only reasonable if called from SetScaleFactorsMuRMuF
   // Is it necessary here ?
   fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();
   if ( !cNLO ) {
      info["SetScaleVariation"]<<"No NLO calculation available."<<endl;
      return true;
   }

   double fScaleFacMuF = cNLO->GetScaleFactor(fScalevar);
   info["SetScaleVariation"]
      <<"Selecting MuF table according to a multiplicative scale factor of the factorization scale of "
      <<fScaleFacMuF<<" times the nominal scale."<<endl;

   // check for threshold corrections.
   if (!BBlocksSMCalc[kThresholdCorrection].empty()) {
      bool lkthc = false;
      for (unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++) {
         if (bUseSMCalc[kThresholdCorrection][i]) {
            lkthc = true;
         }
      }

      if (lkthc) {
         if (abs(fScaleFacMuR-fScaleFacMuF) > DBL_MIN) {
            error["SetScaleVariation."]<<"Threshold corrections only allow for symmetric variations of the renormalization and factorization scales,"<<endl;
            error["SetScaleVariation."]<<"but fScaleFacMuR = "<<fScaleFacMuR<<" is different from fScaleFacMuF = "<<fScaleFacMuF<<", stopped!"<<endl;
            exit(1);
         }
         fastNLOCoeffAddFix* cThC = (fastNLOCoeffAddFix*)B_ThC();
         double fScaleFacMuF2 = cThC->GetScaleFactor(fScalevar);
         if (abs(fScaleFacMuF2-fScaleFacMuF) > DBL_MIN) {
            error["SetScaleVariation."]<<"Scale variations different for NLO and ThC contributions. This should never happen!"<<endl;
            error["SetScaleVariation."]<<"Please do not use this method directly but only via SetScaleFactorsMuRMuF and check the return code!"<<endl;
            exit(1);
         }
      }
   }
   return true;
}


void fastNLOReader::UseHoppetScaleVariations(bool useHoppet){

   if ( FNLO_HOPPET[0] == '\0' ) {
      error["UseHoppetScaleVariation."] << "Hoppet support was not compiled with fastNLO. "
         << "Therefore you can't use Hoppet to calculate the scale variations." <<endl;
      exit(1);
   } else {
      if (useHoppet) {
         if (GetIsFlexibleScaleTable()) {
            info["UseHoppetScaleVariations"]<<"This is a 'flexible-scale' table, therefore you can already choose all desired scale variations without Hoppet."<<endl;
            fUseHoppet = false;
            return;
         }
         fastNLOCoeffAddBase * c = (fastNLOCoeffAddBase*)B_LO();
         if (c->GetIPDFdef1() == 2) {
            error["UseHoppetScaleVariations"] << "Hoppet scale variations not yet implemented for DIS." << std::endl;
            exit(1);
         }

         info["UseHoppetScaleVariations"] << "Hoppet will be used to calculate scale variations." << std::endl;
         fUseHoppet = true;
         HoppetInterface::InitHoppet(*this);
         FillPDFCache(1.);
      }
      else {
         info["UseHoppetScaleVariations"] << "Hoppet will NOT be used to calculate scale variations." << std::endl;
         fUseHoppet = false;
      }
   }
}

//______________________________________________________________________________
bool fastNLOReader::SetContributionON(ESMCalculation eCalc , unsigned int Id , bool SetOn) {
   //! Enable or disable a contribution to be considered in the cross section calculation
   //!  - Use SetOn=true, to switch contribution ON,
   //!  - Use SetOn=false, to switch a contribution off
   //!
   //! Each contribution is identified by an ESMCalculation and by a universal Id.
   //! For all available contributions in your table, call PrintTableInfo().
   //!
   //! The LO contribution can be e.g. addressed by (eCalc=fastNLO::kFixedOrder, Id=0);
   //! The NLO contribution can be e.g. addressed by (eCalc=fastNLO::kFixedOrder, Id=1);
   //!
   //! If an additional additive contribution is switched on, then the PDFCache and AlphasCache
   //! are refilled.

   info["SetContributionON"]<<(SetOn?"Activating":"Deactivating")<<" contribution "<<_ContrName[eCalc]<<" with Id = "<<Id<<endl;

   // sanity checks 1
   // existence of contribution
   if (bUseSMCalc[eCalc].size()<=Id || BBlocksSMCalc[eCalc].size() <=Id) {
      warn["SetContributionON"]
         <<"Contribution "<<_ContrName[eCalc]<<" does not exist in this table, cannot switch it On/Off! Ignoring call."<<endl;
      // set to backed up original value
      return false;
   }

   // backup original value
   bool SetOld = bUseSMCalc[eCalc][Id];
   // set the new value immediately, otherwise GetNScaleVariations(), which is used in FillAlphasCache, will give wrong result.
   bUseSMCalc[eCalc][Id] = SetOn;

   // existence of scale variation for additive contributions (otherwise cache filling will fail!)
   fastNLOCoeffAddBase* c = (fastNLOCoeffAddBase*)BBlocksSMCalc[eCalc][Id];
   if (!GetIsFlexibleScaleTable(c) && !c->GetIAddMultFlag()) {
      int scaleVar = c->GetNpow() == ILOord ? 0 : fScalevar;
      // check that scaleVar is in allowed range, can otherwise lead to segfaults!
      if (scaleVar >= GetNScaleVariations()) {
         warn["SetContributionON"]
            <<"Scale variation "<<scaleVar<<" of contribution "<<_ContrName[eCalc]<<" , Id = "<<Id<<", is > number of available scale variations "<<GetNScaleVariations()<<"! Ignoring call."<<endl;
         // set to backed up original value
         bUseSMCalc[eCalc][Id] = SetOld;
         return false;
      }
   }

   if (!SetOld && SetOn && !fastNLOTools::IsEmptyVector(XSection_LO) ) {
      if (!c->GetIAddMultFlag()) { // if 'new' additive contribution, then refill PDF and alpha_s cache.
         // Fill alpha_s cache
         debug["SetContributionON"]<<"Call FillAlphasCache for contribution eCalc="<<eCalc<<"\tId="<<Id<<endl;
         fAlphasCached = 0;
         FillAlphasCache();
         // Fill PDF cache
         debug["SetContributionON"]<<"Call FillPDFCache for contribution eCalc="<<eCalc<<"\tId="<<Id<<endl;
         fPDFCached = 0;
         FillPDFCache(0.);
      }
   }
   // Needs to be done in the beginning!
   // set the new value
   //   bUseSMCalc[eCalc][Id] = SetOn;
   return true;
}


//______________________________________________________________________________
int fastNLOReader::GetNScaleVariations() const {
   if (GetIsFlexibleScaleTable()) {
      info["GetNScaleVariations"]<<"This is a 'flexible-scale' table, therefore you can choose all desired scale variations."<<endl;
      return 0;
   }

   // Check whether only contributions without extra scale tables (LO, multiplicative) are present
   bool NoExtra = true;
   // Check for maximal scale variation of all active SM calcs with extra scale tables
   // Assume a maximum of 10!
   unsigned int scalevarmax = 10;
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
         fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)BBlocksSMCalc[j][i];
         // Check on contributions with extra scale tables (NLO, threshold corrections)
         int kType  = c->GetIContrFlag1()-1;
         int kOrder = c->GetIContrFlag2()-1;
         debug["GetNScaleVariations"]<<"Contribution type is = "<<kType<<", contribution order is = "<<kOrder<<", contribution switch is = " <<bUseSMCalc[j][i]<<endl;
         // Do not check pQCD LO or multiplicative corrections
         if (bUseSMCalc[j][i] && !c->GetIAddMultFlag() &&
             !(kType == kFixedOrder && kOrder == kLeading)) {
            NoExtra = false;
            if (c->GetNScalevar() < (int)scalevarmax) {
               scalevarmax = c->GetNScalevar();
            }
         }
      }
   }
   if (NoExtra) {scalevarmax = 1;}
   debug["GetNScaleVariations"]<<"Found "<<scalevarmax<<" scale variations."<<endl;
   return scalevarmax;
}


//______________________________________________________________________________
vector < double > fastNLOReader::GetScaleFactors() const {
   if (GetIsFlexibleScaleTable()) {
      info["GetScaleFactors"]<<"This is a 'flexible scale table', therefore you can choose all desired scale variations."<<endl;
      return vector<double>();
   }
   else
      return ((fastNLOCoeffAddFix*)BBlocksSMCalc[kFixedOrder][kNextToLeading])->GetAvailableScaleFactors();
}

//______________________________________________________________________________
string fastNLOReader::GetScaleDescription(const ESMOrder eOrder, int iScale) const {
   //! Get label of scale iScale for order eOrder of the fixed order calculation.
   fastNLOCoeffAddBase* coeff = NULL;
   if (eOrder < (int)BBlocksSMCalc[kFixedOrder].size())
      coeff = (fastNLOCoeffAddBase*) BBlocksSMCalc[kFixedOrder][eOrder];
   else {
      error["GetScaleDescription"]<<"Requested contribution not found." << endl;
      exit(1);
   }
   return coeff->GetScaleDescription(iScale);
}

//______________________________________________________________________________
double fastNLOReader::GetNevt(const ESMOrder eOrder) const {
   //! Get label of scale iScale for order eOrder of the fixed order calculation.
   fastNLOCoeffAddBase* coeff = NULL;
   if (eOrder < (int)BBlocksSMCalc[kFixedOrder].size())
      coeff = (fastNLOCoeffAddBase*) BBlocksSMCalc[kFixedOrder][eOrder];
   else {
      error["GetScaleDescription"]<<"Requested contribution not found." << endl;
      exit(1);
   }
   return coeff->GetNevt();
}
//______________________________________________________________________________
vector < double > fastNLOReader::GetCrossSection() {
   // Get fast calculated NLO cross section
   if (XSection.empty()) CalcCrossSection();
   return XSection;
}

//______________________________________________________________________________
vector< vector < double > > fastNLOReader::GetCrossSection2Dim() {
   //! Get cross section as 2-dimensional vector according to defined binning
   if (GetNumDiffBin() != 2)
      error["GetCrossSection2Dim"]<<"This function is only valid for NDiffBin=2"<<endl;
   // Get fast calculated NLO cross section
   if (XSection.empty()) CalcCrossSection();
   vector< vector < double > > XSection2Dim;
   int k = 0;
   for (unsigned int i = 0; i < GetNDim0Bins(); i++) {
      XSection2Dim.push_back(vector < double >());
      int  NBinDim  = GetNDim1Bins(i);
      for (int j = 0; j < NBinDim; j++) {
         XSection2Dim[i].push_back(XSection[k]);
         k++;
      }
   }
   return XSection2Dim;
}
//______________________________________________________________________________
vector < double > fastNLOReader::GetKFactors() {
   // Get ratio of fast calculated NLO to LO cross section
   if (XSection.empty()) CalcCrossSection();
   return kFactor;
}


//______________________________________________________________________________
vector < double > fastNLOReader::GetQScales(int irelord) {
   // Get XSection weighted Q scale in bin
   if (XSection.empty()) CalcCrossSection();
   if (irelord == 0) {
      return QScale_LO;
   } else {
      return QScale;
   }
}


//______________________________________________________________________________
vector < double > fastNLOReader::GetReferenceCrossSection() {
   // Get reference cross section from direct nlojet++ calculation
   if (XSectionRef.empty() && XSectionRef_s1.empty()) {
      CalcReferenceCrossSection();
   }
   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc == kScale1 && fMuRFunc == kScale1)                   return XSectionRef_s1;
      else if (fMuFFunc == kScale2 && fMuRFunc == kScale2)              return XSectionRef_s2;
      else if (fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean)return XSectionRefMixed;
      else return XSectionRefMixed;
   } else return XSectionRef; // XSectionRef from BlockB-Ref
   return XSectionRef;
}


//______________________________________________________________________________
void fastNLOReader::CalcReferenceCrossSection() {
   //!
   //!  Initialize the internal arrays for the reference cross
   //!  sections with the information from the FastNLO file
   //!

   XSectionRef.clear();
   XSectionRef.resize(NObsBin);
   XSectionRefMixed.clear();
   XSectionRef_s1.clear();
   XSectionRef_s2.clear();
   XSectionRefMixed.resize(NObsBin);
   XSectionRef_s1.resize(NObsBin);
   XSectionRef_s2.resize(NObsBin);

   if (!GetIsFlexibleScaleTable()) {
      fastNLOCoeffAddBase* Coeff_LO_Ref = GetReferenceTable(kLeading);
      fastNLOCoeffAddBase* Coeff_NLO_Ref = GetReferenceTable(kNextToLeading);
      fastNLOCoeffAddBase* Coeff_NNLO_Ref = GetReferenceTable(kNextToNextToLeading);
      if ( Coeff_LO_Ref && Coeff_NLO_Ref && Coeff_NNLO_Ref )
         warn["CalcReferenceCrossSection"]<<"Found NNLO reference cross section. Returning reference of LO+NLO+NNLO.\n";
      if (Coeff_LO_Ref && Coeff_NLO_Ref) {
         for (unsigned int i=0; i<NObsBin; i++) {
            for (int l=0; l<Coeff_LO_Ref->GetNSubproc(); l++) {
               fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)Coeff_LO_Ref;
               int xUnits = c->GetIXsectUnits();
               double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
               XSectionRef[i] +=  c->GetSigmaTilde(i,0,0,0,l) * unit / c->GetNevt(i,l) ; // no scalevariations in LO tables
            }
            for (int l=0; l<Coeff_NLO_Ref->GetNSubproc(); l++) {
               fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)Coeff_NLO_Ref;
               int xUnits = c->GetIXsectUnits();
               double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
               XSectionRef[i] +=  c->GetSigmaTilde(i,fScalevar,0,0,l) * unit / c->GetNevt(i,l);
            }
            if ( Coeff_NNLO_Ref ) {
               for (int l=0; l<Coeff_NNLO_Ref->GetNSubproc(); l++) {
                  fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)Coeff_NNLO_Ref;
                  int xUnits = c->GetIXsectUnits();
                  double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
                  XSectionRef[i] +=  c->GetSigmaTilde(i,fScalevar,0,0,l) * unit / c->GetNevt(i,l);
               }
            }
         }
      }
      else
         warn["CalcReferenceCrossSection"]<<"No reference cross sections for LO and NLO available.\n";
   }
   else {
      for (unsigned int i=0; i<NObsBin; i++) {
         fastNLOCoeffAddFlex* cLO = (fastNLOCoeffAddFlex*)BBlocksSMCalc[kFixedOrder][kLeading];
         int xUnits = cLO->GetIXsectUnits();
         double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
         for (int n=0; n<cLO->GetNSubproc(); n++) {
            XSectionRefMixed[i]             += cLO->SigmaRefMixed[i][n] * unit / cLO->GetNevt(i,n);
            XSectionRef_s1[i]               += cLO->SigmaRef_s1[i][n] * unit / cLO->GetNevt(i,n);
            XSectionRef_s2[i]               += cLO->SigmaRef_s2[i][n] * unit / cLO->GetNevt(i,n);
         }
         fastNLOCoeffAddFlex* cNLO = (fastNLOCoeffAddFlex*)BBlocksSMCalc[kFixedOrder][kNextToLeading];
         xUnits = cNLO->GetIXsectUnits();
         unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
         for (int n=0; n<cNLO->GetNSubproc(); n++) {
            XSectionRefMixed[i]             += cNLO->SigmaRefMixed[i][n] * unit / cNLO->GetNevt(i,n);
            XSectionRef_s1[i]               += cNLO->SigmaRef_s1[i][n] * unit / cNLO->GetNevt(i,n);
            XSectionRef_s2[i]               += cNLO->SigmaRef_s2[i][n] * unit / cNLO->GetNevt(i,n);
         }
         // todo: nnlo reference cross section
      }
   }
}


//______________________________________________________________________________
bool fastNLOReader::PrepareCache() {
   // check pdf cache
   const double PDFcks = CalcNewPDFChecksum();
   if (fPDFCached==0. || (fPDFCached!=0. && fabs(PDFcks/fPDFCached -1.) > 1.e-7)) {
      debug["PrepareCache"]<<"Need to refill PDFCache, since PDFCecksum="<<PDFcks<<" and fPDFCached="<<fPDFCached<<endl;
      FillPDFCache(PDFcks);
   } else  debug["PrepareCache"]<<"No need to refill PDFCache."<<endl;

   // check pdf cache
   if (!fPDFSuccess) {
      error["PrepareCache"]<<"Cannot calculate cross sections. PDF has not been initalized successfully."<<endl;
      return false;
   }

   // check alpha_s cache
   const double asref = CalcReferenceAlphas();
   if (fAlphasCached == 0. || fAlphasCached != asref) {
      debug["PrepareCache"]<<"Need to refill AlphasCache, since fAlphasCached="<<fAlphasCached<<endl;
      FillAlphasCache();
   }
   // do we now have an alphas?
   if (fAlphasCached==0. || fAlphasCached != asref) {
      error["PrepareCache"]<<"Filling of alpha_s cache failed. fAlphasCached="<<fAlphasCached<<"\tasref="<<asref<<endl;
      return false;
   }
   return true;
}


//______________________________________________________________________________
void fastNLOReader::CalcCrossSection() {
   //!
   //!  Initialize the internal arrays with the NLO cross sections
   //!  with the information from the FastNLO file, the pdf and
   //!  the defined alpha_s
   //!
   // xs = (sum(all active add. (perturbative) contr.) + sum(all other active add. contr.) * prod(all active multipl. contr.)
   debug["CalcCrossSection"]<<endl;

   XSection_LO.clear();
   XSection.clear();
   XSection_LO.resize(NObsBin);
   XSection.resize(NObsBin);
   kFactor.clear();
   kFactor.resize(NObsBin);
   QScale_LO.clear();
   QScale.clear();
   QScale_LO.resize(NObsBin);
   QScale.resize(NObsBin);

   // handle alpha_s and PDF Cache
   bool CacheOK = PrepareCache();
   if (!CacheOK) {
      error["CalcCrossSection"]<<"Caching failed. Cannot calculate cross sections."<<endl;
      return;
   }

   // Perturbative (additive) contributions
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
         if ( bUseSMCalc[j][i] ) {
            if ( fastNLOCoeffAddFlex::CheckCoeffConstants(BBlocksSMCalc[j][i],true) )
               CalcCrossSectionv21((fastNLOCoeffAddFlex*)BBlocksSMCalc[j][i]);
            else if ( fastNLOCoeffAddFix::CheckCoeffConstants(BBlocksSMCalc[j][i],true) )
               CalcCrossSectionv20((fastNLOCoeffAddFix*)BBlocksSMCalc[j][i]);
         }
      }
   }

   // Check whether pQCD contributions beyond LO exist and are activated
   bool lknlo = false;
   if (!BBlocksSMCalc[kFixedOrder].empty()) {
      for (unsigned int i = 0 ; i <BBlocksSMCalc[kFixedOrder].size() ; i++) {
         int kOrder = BBlocksSMCalc[kFixedOrder][i]->GetIContrFlag2()-1;
         if (bUseSMCalc[kFixedOrder][i] && kOrder > 0) {
            lknlo = true;
            break;
         }
      }
   }

   // Contributions from the a-posteriori scale variation
   if (!GetIsFlexibleScaleTable() && lknlo) {
      fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();
      if ( fabs(fScaleFacMuF - cNLO->GetScaleFactor(fScalevar)) > DBL_MIN ) {
         if (!fUseHoppet){
            error["CalcCrossSection"] << "Inconsistent choice of chosen factorization scale table and fScaleFacMuF." << endl;
            exit(1);
         }
         CalcAposterioriScaleVariationMuF();
      }
      if ( fabs(fScaleFacMuR - cNLO->GetScaleFactor(fScalevar)) > DBL_MIN ) {
         CalcAposterioriScaleVariationMuR();
      }
   }

   // calculate LO cross sections
   if (GetIsFlexibleScaleTable())
      CalcCrossSectionv21((fastNLOCoeffAddFlex*)B_LO(),true);
   else
      CalcCrossSectionv20((fastNLOCoeffAddFix*)B_LO(),true);


   // non-perturbative corrections (multiplicative corrections)
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
         if (bUseSMCalc[j][i] ) {
            if ( fastNLOCoeffMult::CheckCoeffConstants( BBlocksSMCalc[j][i] , true ) ) {
               fastNLOCoeffMult* cMult = (fastNLOCoeffMult*) BBlocksSMCalc[j][i];
               if ( cMult->GetIContrFlag1() == 4 && cMult->GetIContrFlag2() == 1) {
                  debug["CalcCrossSection"]<<"Multiplying with non-perturbative correction."<<endl;
                  for (unsigned int iB=0; iB<NObsBin; iB++) {
                     XSection[iB] *= cMult->GetMultFactor(iB);
                     //            XSection_LO[iB]     *= BBlocksSMCalc[j][i]->fact[iB];
                  }
               } else {
                  error["CalcCrossSection"]<<"Found unknown multiplicative correction. Printing coeff table and exiting..."<<endl;
                  cMult->Print();
                  exit(1);
               }
            }
         }
      }
   }

   // ---- k-factor calculation ---- //
   debug["CalcCrossSection"]<<"Calculate k-factors: xs/xs_LO"<<endl;
   for (unsigned int i=0; i<NObsBin; i++) {
      kFactor[i] = XSection[i] / XSection_LO[i];
   }

   // ---- Q-scale calculation ---- //
   debug["CalcCrossSection"]<<"Calculate Q-scales: xsQ/xs"<<endl;
   for (unsigned int i=0; i<NObsBin; i++) {
      QScale_LO[i] = QScale_LO[i]/XSection_LO[i];
      QScale[i]    = QScale[i]/XSection[i];
   }

}


//______________________________________________________________________________
void fastNLOReader::CalcAposterioriScaleVariationMuR() {

   fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();
   int scaleVar          = cNLO->GetNpow() == ILOord ? 0 : fScalevar;
   double scalefac       = fScaleFacMuR / cNLO->GetScaleFactor(scaleVar);

   logger.debug["CalcAposterioriScaleVariationMuR"]<<"scalefac="<<scalefac<<endl;
   if ( GetIsFlexibleScaleTable() ) { logger.error["CalcAposterioriScaleVariationMuR"]<<"This function is applicable only to non-flexible scale tables."<<endl; exit(1);}
   fastNLOCoeffAddFix* cLO  = (fastNLOCoeffAddFix*) B_LO();
   vector<double>* XS    = &XSection;
   vector<double>* QS    = &QScale;
   int xUnits = cLO->GetIXsectUnits();
   const double n     = cLO->GetNpow();
   const double L     = log(scalefac);
   //TBD: 5 must be replaced by Nf here!
   const double beta0 = (11.*3.-2.*5)/3.;
   for (unsigned int i=0; i<NObsBin; i++) {
      double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
      int nxmax = cLO->GetNxmax(i);
      for (int j=0; j<cLO->GetTotalScalenodes(); j++) {
         double asnp1 = pow(cLO->AlphasTwoPi_v20[i][j],(n+1)/n);//as^n+1
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<cLO->GetNSubproc(); l++) {
               double clo  = cLO->GetSigmaTilde(i,0,j,k,l) *  cLO->PdfLc[i][j][k][l] * unit / cLO->GetNevt(i,l);
               double xsci = asnp1 * clo * n * L * beta0;
               double mur  = fScaleFacMuR * cLO->GetScaleNode(i,0,j);
               XS->at(i) +=  xsci;
               QS->at(i) +=  xsci*mur;
            }
         }
      }
   }
}

//______________________________________________________________________________
void fastNLOReader::CalcAposterioriScaleVariationMuF() {

   fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();
   int scaleVar          = cNLO->GetNpow() == ILOord ? 0 : fScalevar;
   double scalefac       = fScaleFacMuF / cNLO->GetScaleFactor(scaleVar);

   debug["CalcAposterioriScaleVariationMuF"]<<"scalefac="<<scalefac<<endl;
   if ( GetIsFlexibleScaleTable() ) { error["CalcAposterioriScaleVariationMuF"]<<"This function is only reasonable for non-flexible scale tables."<<endl; exit(1);}
   fastNLOCoeffAddFix* cLO  = (fastNLOCoeffAddFix*) B_LO();
   vector<double>* XS    = &XSection;
   int xUnits = cLO->GetIXsectUnits();
   const double n     = cLO->GetNpow();
   debug["CalcAposterioriScaleVariationMuF"] << "Npow=" << n <<endl;
   for (unsigned int i=0; i<NObsBin; i++) {
      double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
      int nxmax = cLO->GetNxmax(i);
      for (int j=0; j<cLO->GetTotalScalenodes(); j++) {
         double asnp1 = pow(cLO->AlphasTwoPi_v20[i][j],(n+1)/n);//as^n+1
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<cLO->GetNSubproc(); l++) {
               // TODO: Not implemented correctly. Need to fix DIS case.
               double clo  = cLO->GetSigmaTilde(i,0,j,k,l) *(cLO->PdfSplLc1[i][j][k][l] + cLO->PdfSplLc2[i][j][k][l]) * unit / cLO->GetNevt(i,l);
               double xsci = asnp1 * n * log(scalefac) * clo;
               //double xsci = asnp1 * n * log(scalefac) * clo;
               XS->at(i) -= xsci;
            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::CalcCrossSectionv21(fastNLOCoeffAddFlex* c , bool IsLO) {
   //!
   //!  Cross section calculation for DIS and HHC tables in v2.1 format
   //!
   debug["CalcCrossSectionv21"]<<"Npow="<<c->GetNpow()<<"\tIsLO="<<IsLO<<endl;

   vector<double>* XS = IsLO ? &XSection_LO : &XSection;
   vector<double>* QS = IsLO ? &QScale_LO : &QScale;
   // KR: Having different IXsectUnits in different contributions only works when
   //     everything always scaled to Ipublunits (unique per table)
   // Get x section units of each contribution
   int xUnits = c->GetIXsectUnits();
   debug["CalcCrossSectionv21"]<<"Ipublunits = " << Ipublunits << ", xUnits = " << xUnits << endl;

   for (unsigned int i=0; i<NObsBin; i++) {
      double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
      int nxmax = c->GetNxmax(i);
      for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
         for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
            double Q2           = pow(c->GetScaleNode1(i,jS1),2);
            double mur          = CalcMu(kMuR , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuR);
            double muf          = CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
            double mur2         = pow(mur,2);
            double muf2         = pow(muf,2);
            for (int x=0; x<nxmax; x++) {
               for (int n=0; n<c->GetNSubproc(); n++) {
                  double as             = c->AlphasTwoPi[i][jS1][kS2];
                  double pdflc          = c->PdfLcMuVar[i][x][jS1][kS2][n];
                  if (pdflc == 0.) continue;
                  double fac  = as * pdflc * unit;
                  double xsci = c->SigmaTildeMuIndep[i][x][jS1][kS2][n] * fac / c->GetNevt(i,n);
                  if ( c->GetNScaleDep() >= 5 ) {
                     xsci             += c->SigmaTildeMuFDep [i][x][jS1][kS2][n] * log(muf2) * fac / c->GetNevt(i,n);
                     xsci             += c->SigmaTildeMuRDep [i][x][jS1][kS2][n] * log(mur2) * fac / c->GetNevt(i,n);
                     if ( c->GetIPDFdef1() == 2 ) {   // DIS tables use log(mu/Q2) instead of log(mu)
                        xsci -= c->SigmaTildeMuFDep [i][x][jS1][kS2][n] * log(Q2) * fac / c->GetNevt(i,n);
                        xsci -= c->SigmaTildeMuRDep [i][x][jS1][kS2][n] * log(Q2) * fac / c->GetNevt(i,n);
                     }
                     if ( c->GetNScaleDep() >= 6 ) {
                        xsci             += c->SigmaTildeMuRRDep [i][x][jS1][kS2][n] * pow(log(mur2),2) * fac / c->GetNevt(i,n);
                     }
                     if ( c->GetNScaleDep() >= 7 ) {
                        xsci             += c->SigmaTildeMuFFDep [i][x][jS1][kS2][n] * pow(log(muf2),2) * fac / c->GetNevt(i,n);
                        xsci             += c->SigmaTildeMuRFDep [i][x][jS1][kS2][n] * log(mur2) * log(muf2) * fac / c->GetNevt(i,n);
                     }
                  }
                  XS->at(i)   += xsci;
                  QS->at(i)   += xsci*mur;
               }
            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::CalcCrossSectionv20(fastNLOCoeffAddFix* c , bool IsLO) {
   //!
   //!  Cross section calculation in v2.0 format
   //!
   debug["CalcCrossSectionv20"]<<"Npow="<<c->GetNpow()<<"\tIsLO="<<IsLO<<endl;

   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar;
   vector<double>* XS    = IsLO ? &XSection_LO : &XSection;
   vector<double>* QS    = IsLO ? &QScale_LO : &QScale;
   // KR: Having different IXsectUnits in different contributions only works when
   //     everything always scaled to Ipublunits (unique per table)
   // Get x section units of each contribution
   int xUnits = c->GetIXsectUnits();
   debug["CalcCrossSectionv20"]<<"Ipublunits = " << Ipublunits << ", xUnits = " << xUnits << endl;

   for (unsigned int i=0; i<NObsBin; i++) {
      double unit = RescaleCrossSectionUnits(BinSize[i], xUnits);
      int nxmax = c->GetNxmax(i);
      for (int j=0; j<c->GetTotalScalenodes(); j++) {
         double scalefac = fScaleFacMuR/c->GetScaleFactor(scaleVar);
         double mur      = scalefac * c->GetScaleNode(i,scaleVar,j);
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<c->GetNSubproc(); l++) {
               double xsci     = c->GetSigmaTilde(i,scaleVar,j,k,l) *  c->AlphasTwoPi_v20[i][j]  * c->PdfLc[i][j][k][l] * unit / c->GetNevt(i,l);
               XS->at(i)      +=  xsci;
               QS->at(i)      +=  xsci*mur;
            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::SetUnits(EUnits Unit) {
   if (fUnits != Unit) {
      fUnits  = Unit;
      //CalcCrossSection();
   } else {
      // nothing todo
   }
}


//______________________________________________________________________________
void fastNLOReader::FillAlphasCache() {
   debug["FillAlphasCache"]<<endl;
   //!
   //!  Fill the internal alpha_s cache.
   //!  This is usally called automatically. Only if you
   //!  make use of ReFillCache==false options, you have
   //!  to take care of this filling by yourself.
   //!

   // check if the alpha_s value is somehow reasonable
   debug["FillAlphasCache"]<<"Sanity check!"<<endl;
   TestAlphas();

   // is there a need for a recalculation?
   const double asNew = CalcReferenceAlphas();
   if (asNew == fAlphasCached) {
      debug["FillAlphasCache"]<<"No need for a refilling of AlphasCache. asNew==fAlphasCached="<<asNew<<endl;
   } else {
      fAlphasCached = asNew;
      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            // Check that this contribution type j and no. i should actually be used
            // Otherwise deactivation of e.g. threshold corr. is not respected here
            if (bUseSMCalc[j][i] ) {
               fastNLOCoeffBase* c = BBlocksSMCalc[j][i];
               if ( fastNLOCoeffAddFlex::CheckCoeffConstants(c,true) )
                  FillAlphasCacheInBlockBv21((fastNLOCoeffAddFlex*)c);
               else if ( fastNLOCoeffAddFix::CheckCoeffConstants(c,true) )
                  FillAlphasCacheInBlockBv20((fastNLOCoeffAddFix*)c);
               else if ( fastNLOCoeffMult::CheckCoeffConstants(c,true) )
                  info["FillAlphasCache"]<<"Nothing to be done for multiplicative contribution."<<endl;
               else {
                  error["FillAlphasCache"]<<"Could not identify contribution. Printing."<<endl;
                  c->Print();
               }
            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::FillAlphasCacheInBlockBv20(fastNLOCoeffAddFix* c) {
   //!
   //!  Internal method for filling alpha_s cache
   //!

   // todo: the flag IScaleDep should also indicate whether this contribution may contain scale variations
   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar;

   // Sanity check that scaleVar is in allowed range
   // For thresh. corr. can otherwise lead to inf and then segfault!
   if (scaleVar >= GetNScaleVariations()) {
      error<<"Trying to refresh cache for non-existing scale variation no. "<<scaleVar<<" while only "<<GetNScaleVariations()<<" exist in total. Exiting."<<endl;
      exit(1);
   }
   double scalefac       = fScaleFacMuR/c->GetScaleFactor(scaleVar);
   debug["FillAlphasCacheInBlockBv20"]<<"scalefac="<<scalefac<<"\tscaleVar="<<scaleVar<<endl;

   for (unsigned int i=0; i<NObsBin; i++) {
      for (int j=0; j<c->GetTotalScalenodes(); j++) {
         double mur        = scalefac * c->GetScaleNode(i,scaleVar,j);
         double as         = CalcAlphas(mur);
         c->AlphasTwoPi_v20[i][j] = pow(as/TWOPI , c->GetNpow());
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::FillAlphasCacheInBlockBv21(fastNLOCoeffAddFlex* c) {
   //!
   //!  Internal method for filling alpha_s cache
   //!

   for (unsigned int i=0; i<NObsBin; i++) {
      for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
         for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
            double mur              = CalcMu(kMuR , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuR);
            double as               = CalcAlphas(mur);
            double alphastwopi      = pow(as/TWOPI, c->GetNpow() );
            c->AlphasTwoPi[i][jS1][kS2] = alphastwopi;
         }
      }
   }
}


//______________________________________________________________________________
double fastNLOReader::CalcAlphas(double Q) {
   //!
   //!  Internal method for calculating the alpha_s(mu)
   //!
   return EvolveAlphas(Q);
}


//______________________________________________________________________________
double fastNLOReader::CalcReferenceAlphas() {
   double mu = 0;
   if (GetIsFlexibleScaleTable()) {
      if (fMuRFunc==kExtern) mu = (*Fct_MuR)(91.,1.)*(fScaleFacMuR+0.1);
      else mu = 91.1876111111+(fMuRFunc*0.1)+(fScaleFacMuR);
   } else mu = 91.187611111115*(fScaleFacMuR+0.1)+fScalevar*0.1;
   double as = CalcAlphas(mu);
   if (std::isnan(as)) {
      error["CalcReferenceAlphas"]<<"Reference alphas is a 'nan' for scale mu="<<mu<<endl;
      //exit(1);
   }
   return as;
}


//______________________________________________________________________________
double fastNLOReader::CalcNewPDFChecksum() {
   //! calculate a PDF checksum to
   //! decide, whether PDF cache has to be refilled

   // init PDF and check success
   debug["CalcNewPDFChecksum"]<<"Call InitPDF() in user module."<<endl;
   fPDFSuccess = InitPDF();
   debug["CalcNewPDFChecksum"]<<"Return value InitPDF() = "<<fPDFSuccess<<endl;
   if (!fPDFSuccess) {
      warn["CalcPDFChecksum"]<<"PDF initialization failed. Please check PDF interface in your FastNLO user module."<<endl;
      return 0.;
   }

   // calculate checksum for some scales and flavors
   double muf = 0;
   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc==kExtern) muf = (*Fct_MuF)(91.,10.)/91.*(fScaleFacMuF+0.5) ;
      else muf = (91.1+0.1*fMuFFunc)/91.+fScaleFacMuF;
   } else muf=(fScaleFacMuF+0.1)+fScalevar*0.1;
   double cks = CalcChecksum(muf);
   return cks;
}


//______________________________________________________________________________
double fastNLOReader::CalcChecksum(double mufac) {
   //! caculate a checksum from the PDF in order to check
   //! if the PDF has changed. This is mandatory
   //! since the old LHAPDF code is written in fortran
   //! and may PDFs may change without any notice.
   debug["CalcChecksum"]<<"Calculate checksum of 13 flavors, 3 mu_f values, and 3 x-values, for scalefac="<<mufac<<endl;
   double cks = 0;
   vector<double> xfx(13);
   const double mf[3] = { 3,10,91.18};
   const double x[3] = {1.e-1,1.e-2,1.e-3};
   for (int jf = 0 ; jf<3 ; jf++) {
      double mu = mf[jf]* mufac;//(fScaleFacMuF+0.1)+fScalevar*0.1;
      for (int ix = 0 ; ix<3 ; ix++) {
         xfx = GetXFX(x[ix],mu);
         for (unsigned int fl = 0 ; fl<xfx.size() ; fl++) {
            cks+=xfx[fl];
         }
      }
   }
   debug["CalcChecksum"]<<"Calculated checksum = "<<cks<<endl;
   return cks;
}


//______________________________________________________________________________
bool fastNLOReader::TestAlphas() {
   //! Test if the alpha_s evolution provided by the user
   //! yields realistic results.
   const double as = CalcAlphas(91.18);
   if (as < 0.01 || as > 0.5) {
      warn["TestAlphas"]<<"The alphas value, returned by the user class seems to be unreasonably small/large."<<endl;
      warn["TestAlphas"]<<"The evolution code calculated alphas(Mz~91.18GeV) = "<<as<<endl;
      return false;
   }
   debug["TestAlphas"]<<"Sanity check of alpha_s(MZ=91.18) = "<<as<<endl;
   return true;
}


//______________________________________________________________________________
bool fastNLOReader::TestXFX() {
   vector<double> pdftest = GetXFX(1.e-2,10);
   if (pdftest.size() <= 13) {
      error["TestXFX"]<<"The pdf array must have the size of 13 flavors.\n";
      return false;
   }
   // if ( pdftest[6] == 0. )printf("fastNLOReader. Warning. There seems to be no gluon in the pdf.\n");
   // double sum = 0;
   // for ( int i = 0 ; i<13 ; i++ ) sum+=fabs(pdftest[i]);
   // if ( sum== 0. ) printf("fastNLOReader. Error. All 13 pdf probabilities are 0. There might be sth. wrong in the pdf interface. Please check FastNLOUser::GetXFX().\n");
   for (int i = 0 ; i<13 ; i++) {
      if (pdftest[i] > 1.e10 || (pdftest[i] < 1.e-10 && pdftest[i] > 1.e-15)) {
         warn["TestXFX"]<<"The pdf probability of the "<<i<<"'s flavor seeems to be unreasonably large/small (pdf="<<pdftest[i]<<").\\n";
      }
   }
   return true;
}



//______________________________________________________________________________
void fastNLOReader::FillPDFCache(double chksum) {
   debug["FillPDFCache"]<<"Passed chksum="<<chksum<<". Do not recalculate checksum (which calls InitPDF()) if chksum!=0."<<endl;
   //!
   //!  Fill the internal pdf cache.
   //!  This function has to be called by the user, since the
   //!  pdf parameters and evolutions are calculated externally.
   //!

   // reset checknum
   // check if the alpha_s value is somehow reasonable
   double PDFnew = chksum;
   if (chksum == 0.) {
      debug["FillPDFCache"]<<"Calculate Checksum!"<<endl;
      PDFnew = CalcNewPDFChecksum();
      if (PDFnew==0.) {
         warn["FillPDFCache"]<<"PDF Checksum is zero."<<endl;
      }
      debug["FillPDFCache"]<<"PDF Checksum = "<<PDFnew<<endl;
   }

   // is there a need for a recalculation?
   if (fPDFCached != 0. && fabs(PDFnew/fPDFCached - 1.) < 1.e-7) {
      debug["FillPDFCache"]<<"No need for a refilling of PDFCache. fPDFCached=RefreshPDFChecksum()"<<PDFnew<<endl;
   } else {
      debug["FillPDFCache"]<<"Refilling PDF cache"<<endl;
      fPDFCached = PDFnew;

      // check (or not) if the pdf is somehow reasonable
      TestXFX();
      if ( FNLO_HOPPET[0] != '\0' ) {
         if (fUseHoppet){
            //Also refill Hoppet cache and assign new PDF
            HoppetInterface::InitHoppet(*this);
         }
      }

      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
         for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
            // Check that this contribution type j and no. i should actually be used
            // Otherwise deactivation of e.g. threshold corr. is not respected here
            if (bUseSMCalc[j][i] ) {
               fastNLOCoeffBase* c = BBlocksSMCalc[j][i];
               if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
                  fastNLOCoeffAddBase* c = (fastNLOCoeffAddBase*)BBlocksSMCalc[j][i];
                  // linear: DIS-case
                  // ---- DIS ---- //
                  if (c->GetIPDFdef1() == 2) {
                     if (c->GetNPDFDim() == 0) {
                        if (!GetIsFlexibleScaleTable(c))
                           FillBlockBPDFLCsDISv20((fastNLOCoeffAddFix*)c);
                        else
                           FillBlockBPDFLCsDISv21((fastNLOCoeffAddFlex*)c);
                     }
                  }
                  // ---- pp ---- //
                  else if (c->GetIPDFdef1() == 3) {
                     if (!GetIsFlexibleScaleTable(c)) FillBlockBPDFLCsHHCv20((fastNLOCoeffAddFix*)c);
                     else FillBlockBPDFLCsHHCv21((fastNLOCoeffAddFlex*)c);
                  } else {
                     error["FillPDFCache"]<<"IPDFdef of tables must be 1 or 2.\n";
                  }
               }
               else if ( fastNLOCoeffMult::CheckCoeffConstants(c,true) ) {
                  info["FillPDFCache"]<<"Nothing to be done for multiplicative contribution."<<endl;
               }
               else {
                  error["FillPDFCache"]<<"Could not identify contribution. Printing."<<endl;
                  c->Print();
               }
            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::FillBlockBPDFLCsDISv20(fastNLOCoeffAddFix* c) {
   //! Fill member variables in fastNLOCoeffAddFix with PDFCache
   debug["FillBlockBPDFLCsDISv20"]<<endl;
   // todo: flag IScaleDep should indicate whether scale variations may exist or not.
   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar;
   double scalefac       = (c->GetScaleFactor(scaleVar) == fScaleFacMuF) ? 1. : fScaleFacMuF;
   vector<double> xfx(13); // PDFs of all partons
   vector<double> xfxspl(13); // PDFs splitting functions of all partons
   if (!GetIsFlexibleScaleTable(c)) {
      for (unsigned int i=0; i<NObsBin; i++) {
         int nxmax = c->GetNxmax(i);
         for (int j=0; j<c->GetNScaleNode(); j++) {
            for (int k=0; k<nxmax; k++) {
               double xp     = c->GetXNode1(i,k);
               double muf    = scalefac * c->GetScaleNode(i,scaleVar,j);
               xfx = GetXFX(xp,muf);

               if ( FNLO_HOPPET[0] != '\0' ) {
                  if (fUseHoppet)
                     xfxspl        = HoppetInterface::GetSpl(xp,muf);
               }
               c->PdfLc[i][j][k] = CalcPDFLinearCombination(c,xfx);
               if (fUseHoppet){
                  c->PdfSplLc1[i][j][k] = CalcPDFLinearCombination(c, xfxspl);
               }
               //                vector < double > buffer = CalcPDFLinearCombDIS(xfx , c->GetNSubproc());
               //                for (int l=0; l<c->GetNSubproc(); l++) {
               //                   c->PdfLc[i][j][k][l] = buffer[l];
               //                }
            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::FillBlockBPDFLCsDISv21(fastNLOCoeffAddFlex* c) {
   //! Fill member variables in fastNLOCoeffAddFlex with PDFCache
   debug["FillBlockBPDFLCsDISv21"]<<endl;//<<"CoeffTable = "<<endl;

   if (c->PdfLcMuVar.empty()) {
      error<< "PdfLcMuVar is empty in CoeffTable. Printing and exiting."<<endl;
      c->Print();
      exit(1);
   }

   for (unsigned int i=0; i<NObsBin; i++) {
      // speed up! if mu_f is only dependent on one variable, we can safe the loop over the other one
      for (int x=0; x<c->GetNxmax(i); x++) {
         //double xp = c->GetXNode1(i,x);
         double xp = c->GetXNode1(i,x);
         if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2) {   // that't the standard case!
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
                  double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
                  c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,GetXFX(xp,muf));
                  //c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombDIS(GetXFX(xp,muf) , c->GetNSubproc() );
               }
            }
         }
         else if (fMuFFunc == kScale2) { // speed up
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               double muf = CalcMu(kMuF , 0 ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
               //vector < double > buffer = CalcPDFLinearCombDIS(GetXFX(xp,muf) , c->GetNSubproc() );
               vector<double > buffer = CalcPDFLinearCombination(c,GetXFX(xp,muf));
               for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
                  c->PdfLcMuVar[i][x][jS1][kS2] = buffer;
               }
            }
         }
         else if (fMuFFunc == kScale1) { // speed up
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) , 0 , fScaleFacMuF);
               //vector < double > buffer = CalcPDFLinearCombDIS(GetXFX(xp,muf) , c->GetNSubproc() );
               vector<double > buffer = CalcPDFLinearCombination(c,GetXFX(xp,muf));
                for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
                   c->PdfLcMuVar[i][x][jS1][kS2] = buffer;
               }
            }
         }
      }
   }
   debug["FillBlockBPDFLCsDISv21"]<<"done." <<endl;

}


//______________________________________________________________________________
void fastNLOReader::FillBlockBPDFLCsHHCv20(fastNLOCoeffAddFix* c) {
   //! Fill member variables in fastNLOCoeffAddFix with PDFCache
   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar; // Use IScaleDep
   double scalefac       = fScaleFacMuF/c->GetScaleFactor(scaleVar);
   debug["FillBlockBPDFLCsHHCv20"]<<"scalefac="<<scalefac<<endl;

   bool IsPPBar;
   // ----- if ppbar ---- //
   if (c->NPDFPDG[0] == -c->NPDFPDG[1] && abs(c->NPDFPDG[0]) == 2212) {
      IsPPBar = true;
   }
   // ----- if pp ---- //
   else if (c->NPDFPDG[0] == c->NPDFPDG[1] && (c->NPDFPDG[0] == 2212 || c->NPDFPDG[1] == 2212)) {
      IsPPBar = false;
   }
   // ----- anything else ---- //
   else {
      error<<"Found beam particles to have PDG codes " << c->NPDFPDG[0] << " and " << c->NPDFPDG[1] << ",\n";
      error<<"but cannot deal with tables other than pp or ppbar, aborting! \n";
      exit(1);
   }


   // half matrix notation
   if ( c->GetNPDFDim() == 1 ) {
      vector < vector < double > > xfx; // PDFs of all partons
      vector < vector < double > > xfxspl; // PDFs splitting functions of all partons
      for (unsigned int i=0; i<NObsBin; i++) {
         int nxmax = c->GetNxmax(i);
         int nxbins1 = c->GetNxtot1(i); // number of columns in half matrix
         xfx.resize(nxbins1);

         if (fUseHoppet){
            xfxspl.resize(nxbins1);
         }
         for (int j=0; j<c->GetNScaleNode(); j++) {
            // determine all pdfs of hadron1
            for (int k=0; k<nxbins1; k++) {
               double xp     = c->GetXNode1(i,k);
               double muf    = scalefac * c->GetScaleNode(i,scaleVar,j);
               xfx[k]        = GetXFX(xp,muf);
               if ( FNLO_HOPPET[0] != '\0' ) {
                  if (fUseHoppet)
                     xfxspl[k]        = HoppetInterface::GetSpl(xp,muf);
               }
            }
            int x1bin = 0;
            int x2bin = 0;
            // half-matrix notation
            for (int k=0; k<nxmax; k++) {
               // Original code calling (x2, x1) cancelling the inverted naming below --> OK in original version
               //               c->PdfLc[i][j][k] = CalcPDFLinearCombination(c,xfx[x2bin],xfx[x1bin], IsPPBar);
               // Fixed code: With correct naming of x1bin, x2bin below (x1, x2) has to be called
               c->PdfLc[i][j][k] = CalcPDFLinearCombination(c,xfx[x1bin],xfx[x2bin], IsPPBar);
               // TODO: Georg was using (x1, x2) that was wrong with original x1bin, x2bin naming
               // TODO: Check with Georg what is correct now after the fix
               if (fUseHoppet){
                  c->PdfSplLc1[i][j][k] = CalcPDFLinearCombination(c, xfx[x1bin], xfxspl[x2bin], IsPPBar);
                  c->PdfSplLc2[i][j][k] = CalcPDFLinearCombination(c, xfxspl[x1bin], xfx[x2bin], IsPPBar);
               }
               // Original code: But x1bin, x2bin are exchanged with respect to GetXIndex for filling.
               //                This is wrong!
               // The 2-dim. x1, x2 bins are mapped onto ix like this:
               //    [0,0] --> [0]
               //    [1,0] --> [1]
               //    [1,1] --> [2]
               //    [2,0] --> [3]
               //    etc.
               // x1bin++;
               // if (x1bin>x2bin) {
               //    x1bin = 0;
               //    x2bin++;
               // Invert naming of x1bin and x2bin
               x2bin++;
               if (x2bin>x1bin) {
                  x2bin = 0;
                  x1bin++;
               }
            }
         }
      }
   }

   // full matrix notation
   else if ( c->GetNPDFDim() == 2 ){
      vector < vector < double > > xfx1; // PDFs of all partons
      vector < vector < double > > xfx2; // PDFs of all partons
      vector < vector < double > > xfxspl1; // PDFs splitting functions of all partons
      vector < vector < double > > xfxspl2; // PDFs splitting functions of all partons
      for (unsigned int i=0; i<NObsBin; i++) {
         int nxmax = c->GetNxmax(i);
         int nxbins1 = c->GetNxtot1(i); // number of xnodes ( == nxmax / Nxtot2[i] )
         int nxbins2 = c->GetNxtot2(i); // number of xnodes ( == nxmax / Nxtot1[i] )
         xfx1.resize(nxbins1);
         xfx2.resize(nxbins2);
         if (fUseHoppet){
            xfxspl1.resize(nxbins1);
            xfxspl2.resize(nxbins1);
         }
         for (int j=0; j<c->GetNScaleNode(); j++) {
            // determine all pdfs of hadron1
            double muf    = scalefac * c->GetScaleNode(i,scaleVar,j);
            for (int k=0; k<nxbins1; k++) {
               double xp     = c->GetXNode1(i,k);
               xfx1[k]        = GetXFX(xp,muf);
               if ( FNLO_HOPPET[0] != '\0' ) {
                  if (fUseHoppet)
                     xfxspl1[k]        = HoppetInterface::GetSpl(xp,muf);
               }

            }
            // determine all pdfs of hadron2
            for (int k=0; k<nxbins2; k++) {
               double xp     = c->GetXNode2(i,k);
               xfx2[k]       = GetXFX(xp,muf);
               if ( FNLO_HOPPET[0] != '\0' ) {
                  if (fUseHoppet)
                     xfxspl2[k]        = HoppetInterface::GetSpl(xp,muf);
               }

            }
            // full matrix notation
            for (int k=0; k<nxmax; k++) {
               int x1bin = k % c->GetNxtot1(i);
               int x2bin = k / c->GetNxtot1(i);
               c->PdfLc[i][j][k] = CalcPDFLinearCombination(c,xfx1[x1bin],xfx2[x2bin], IsPPBar);
               if (fUseHoppet){
                  c->PdfSplLc1[i][j][k] = CalcPDFLinearCombination(c, xfx1[x1bin], xfxspl2[x2bin], IsPPBar);
                  c->PdfSplLc2[i][j][k] = CalcPDFLinearCombination(c, xfxspl1[x1bin], xfx2[x2bin], IsPPBar);
               }

            }
         }
      }
   }
}


//______________________________________________________________________________
void fastNLOReader::FillBlockBPDFLCsHHCv21(fastNLOCoeffAddFlex* c) {
   //! Fill member variables in fastNLOCoeffAddFlex with PDFCache
   //! The calculation is improved, if the factorization scale is
   //! calculated from only one scale variable (i.e. kScale1 or kScale2)
   debug["FillBlockBPDFLCsHHCv21"]<<endl;
   if (c->PdfLcMuVar.empty()) {
      cout<< "PdfLcMuVar in CoeffTable is not accessible (resized)."<<endl;
      exit(1);
   }

   bool IsPPBar;
   // ----- if ppbar ---- //
   if (c->NPDFPDG[0] == -c->NPDFPDG[1] && abs(c->NPDFPDG[0]) == 2212) {
      IsPPBar = true;
   }
   // ----- if pp ---- //
   else if (c->NPDFPDG[0] == c->NPDFPDG[1] && (c->NPDFPDG[0] == 2212 || c->NPDFPDG[1] == 2212)) {
      IsPPBar = false;
   }
   // ----- anything else ---- //
   else {
      error<<"Found beam particles to have PDG codes " << c->NPDFPDG[0] << " and " << c->NPDFPDG[1] << ",\n";
      error<<"but cannot deal with tables other than pp or ppbar, aborting! \n";
      exit(1);
   }


   // half-matrix notation
   if ( c->GetNPDFDim() == 1 ) {
      vector < vector < double > > xfx; // PDFs of all partons
      for (unsigned int i=0; i<NObsBin; i++) {
         int nxmax = c->GetNxmax(i);
         int nxbins1 = c->GetNxtot1(i); // number of columns in half matrix
         xfx.resize(nxbins1);
         if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2)  {   // that't the standard case!
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
                  // determine all pdfs of hadron1
                  for (int k=0; k<nxbins1; k++) {
                     double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
                     double xp   = c->GetXNode1(i,k);
                     xfx[k] = GetXFX(xp,muf);
                  }
                  int x1bin = 0;
                  int x2bin = 0;
                  for (int x=0; x<nxmax; x++) {
                     // CalcPDFLinearCombination calculats Anti-proton from proton
                     c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx[x1bin],xfx[x2bin], IsPPBar);
                     x2bin++;
                     if (x2bin>x1bin) {
                        x2bin = 0;
                        x1bin++;
                     }
                  }
               }
            }
         }
         else if (fMuFFunc == kScale2) {   // speed up
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               // determine all pdfs of hadron1
               for (int k=0; k<nxbins1; k++) {
                  double muf = CalcMu(kMuF , 0 ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
                  double xp     = c->GetXNode1(i,k);
                  xfx[k] = GetXFX(xp,muf);
               }
               for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
                  int x1bin = 0;
                  int x2bin = 0;
                  for (int x=0; x<nxmax; x++) {
                     c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx[x1bin],xfx[x2bin], IsPPBar);
                     x2bin++;
                     if (x2bin>x1bin) {
                        x2bin = 0;
                        x1bin++;
                     }
                  }
               }
            }
         }
         else if (fMuFFunc == kScale1) {   // speed up
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               // determine all pdfs of hadron1
               for (int k=0; k<nxbins1; k++) {
                  double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) , 0 , fScaleFacMuF);
                  double xp     = c->GetXNode1(i,k);
                  xfx[k] = GetXFX(xp,muf);
               }
               for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
                  int x1bin = 0;
                  int x2bin = 0;
                  for (int x=0; x<nxmax; x++) {
                     c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx[x1bin],xfx[x2bin], IsPPBar);
                     x2bin++;
                     if (x2bin>x1bin) {
                        x2bin = 0;
                        x1bin++;
                     }
                  }
               }
            }
         }
      }
   }

   // full-matrix notation
   else if ( c->GetNPDFDim() == 2 ) {
      vector < vector < double > > xfx1; // hadron1
      vector < vector < double > > xfx2; // hadron2
      for (unsigned int i=0; i<NObsBin; i++) {
         int nxmax = c->GetNxmax(i);
         int nxbins1 = c->GetNxtot1(i); // number of xnodes ( == nxmax / Nxtot2[i] )
         int nxbins2 = c->GetNxtot2(i); // number of xnodes ( == nxmax / Nxtot1[i] )
         xfx1.resize(nxbins1);
         xfx2.resize(nxbins2);
         if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2)  {   // that't the standard case!
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
                  // determine all pdfs of hadron1
                  double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
                  for (int k=0; k<nxbins1; k++) {
                     double xp   = c->GetXNode1(i,k);
                     xfx1[k] = GetXFX(xp,muf);
                  }
                  // determine all pdfs of hadron2
                  for (int k=0; k<nxbins2; k++) {
                     double xp   = c->GetXNode2(i,k);
                     xfx2[k] = GetXFX(xp,muf);
                  }
                  for (int x=0; x<nxmax; x++) {
                     // CalcPDFLinearCombination calculats Anti-proton from proton
                     int x1bin = x % c->GetNxtot1(i);
                     int x2bin = x / c->GetNxtot1(i);
                     c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx1[x1bin],xfx2[x2bin], IsPPBar);
                  }
               }
            }
         }
         else if (fMuFFunc == kScale2) {   // speed up
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               double muf = CalcMu(kMuF , 0 ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
               // determine all pdfs of hadron1
               for (int k=0; k<nxbins1; k++) {
                  double xp  = c->GetXNode1(i,k);
                  xfx1[k] = GetXFX(xp,muf);
               }
               // determine all pdfs of hadron2
               for (int k=0; k<nxbins2; k++) {
                  double xp  = c->GetXNode2(i,k);
                  xfx2[k] = GetXFX(xp,muf);
               }
               for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
                  for (int x=0; x<nxmax; x++) {
                     int x1bin = x % c->GetNxtot1(i);
                     int x2bin = x / c->GetNxtot1(i);
                     c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx1[x1bin],xfx2[x2bin]);
                  }
               }
            }
         }
         else if (fMuFFunc == kScale1) {   // speed up
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               // determine all pdfs of hadron1
               double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) , 0 , fScaleFacMuF);
               for (int k=0; k<nxbins1; k++) {
                  double xp   = c->GetXNode1(i,k);
                  xfx1[k] = GetXFX(xp,muf);
               }
               // determine all pdfs of hadron2
               for (int k=0; k<nxbins2; k++) {
                  double xp   = c->GetXNode2(i,k);
                  xfx2[k] = GetXFX(xp,muf);
               }
               for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
                  for (int x=0; x<nxmax; x++) {
                     int x1bin = x % c->GetNxtot1(i);
                     int x2bin = x / c->GetNxtot1(i);
                     c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx1[x1bin],xfx2[x2bin]);
                  }
               }
            }
         }
      }
   }

}


//______________________________________________________________________________
void fastNLOReader::SetExternalFuncForMuR(double(*Func)(double,double)) {
   if (!GetIsFlexibleScaleTable()) {
      warn["SetExternalFuncForMuR"]<<"This is not a flexible-scale table and SetExternalFuncForMuR has no impact.\n";
      man<<"Please use a flexible-scale table, if you want to change your scale definition.\n";
      return;
   }

   Fct_MuR = Func;
   SetFunctionalForm(kExtern , kMuR);
   info["SetExternalFuncForMuR"]<<"Testing external function:"<<endl;
   info<<"Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = "<<(*Fct_MuR)(1,1)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = "<<(*Fct_MuR)(91.1876,91.1876)<<endl;
   info<<"Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = "<<(*Fct_MuR)(1,91.1876)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = "<<(*Fct_MuR)(91.1876,1)<<endl;
}


//______________________________________________________________________________
void fastNLOReader::SetExternalFuncForMuF(double(*Func)(double,double)) {
   if (!GetIsFlexibleScaleTable()) {
      warn["SetExternalFuncForMuF"]<<"This is not a flexible-scale table and SetExternalFuncForMuF has no impact.\n";
      man<<"Please use a flexible-scale table, if you want to change your scale definition.\n";
      return;
   }

   Fct_MuF = Func;
   SetFunctionalForm(kExtern , kMuF);
   info["SetExternalFuncForMuF"]<<"Testing external function:"<<endl;
   info<<"Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = "<<(*Fct_MuF)(1,1)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = "<<(*Fct_MuF)(91.1876,91.1876)<<endl;
   info<<"Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = "<<(*Fct_MuF)(1,91.1876)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = "<<(*Fct_MuF)(91.1876,1)<<endl;
}


//______________________________________________________________________________
void fastNLOReader::SetFunctionalForm(EScaleFunctionalForm func , fastNLO::EMuX MuX) {
   //!
   //!  For MuVar tables this method sets the functional form of
   //!  the renormalization or the factorization scale.
   //!     func:  Choose a pre-defined function
   //!     kMuX:  is it for mu_r or for mu_f ?
   //!

   if (!GetIsFlexibleScaleTable()) {
      warn<<"This is not a flexible-scale table. SetFunctionalForm cannot be used.\n";
      return;
   }

   // ---- setting scale ---- //
   if (MuX == kMuR) fMuRFunc = func;
   else fMuFFunc = func;


   // ---- cross check ---- //
   if (func == kScale2 || func == kQuadraticSum ||  func == kQuadraticMean || func == kQuadraticSumOver4 ||
       func == kLinearMean || func == kLinearSum  ||  func == kScaleMax || func == kScaleMin ||
       func == kProd || func == kExpProd2 ) {

      fastNLOCoeffAddFlex* cNLO = (fastNLOCoeffAddFlex*)B_NLO();
      if ( !cNLO ) cNLO = (fastNLOCoeffAddFlex*)B_LO(); //crash safe
      int nnode = cNLO->GetNScaleNode2(0);
      if (nnode <= 3) {
         error<<"There is no second scale variable available in this table. Using fastNLO::kScale1 only.\n";
         SetFunctionalForm(kScale1,MuX);
      }
      for (unsigned int i=0; i<NObsBin; i++) {
         nnode = cNLO->GetNScaleNode2(i);
         if (nnode < 4) {
            warn<<"Scale2 has only very little nodes (n="<<nnode<<") in bin "<<i<<".\n";
         }
      }
   }
   //PrintScaleSettings(MuX);//not yet ported to v2.2
}


//______________________________________________________________________________
void fastNLOReader::SetMuRFunctionalForm(EScaleFunctionalForm func) {
   SetFunctionalForm(func,kMuR);
}


//______________________________________________________________________________
void fastNLOReader::SetMuFFunctionalForm(EScaleFunctionalForm func) {
   SetFunctionalForm(func,kMuF);
}


//______________________________________________________________________________



bool fastNLOReader::SetScaleFactorsMuRMuF(double xmur, double xmuf) {
   debug["SetScaleFactorsMuRMuF"]<<"Setting to scale factors xmur = "<<xmur<<" and xmuf = "<<xmuf<<endl;
   /**
   // Set renormalization and factorization scale factors simultaneously for scale variations in all v2 tables.
   // You have to ReFill your cache!
   // This is done automatically, but if you want to do it by yourself set ReFillCache = false.
   //
   // The function aborts the whole program if non-sensical scale factors < 1.E-6 are requested.
   // The function returns true if the requested scale factors can be used with the available table:
   //
   //   If it is NOT a flexibleScaleTable and there is no NLO scalevar table for xmuf,
   //   xmur and xmuf are unchanged, a warning is printed and the function returns false!
   //   If threshold corrections are selected, then
   //   - only symmetric scale variations, i.e. xmur / xmuf = 1., are allowed,
   //   - the scale variations for xmuf must be stored in IDENTICAL order
   //     for the NLO and the threshold corrections (there is only one fScalevar!)
   //   If either is not the case, xmur and xmuf are unchanged,
   //   a warning is printed and the function returns false!
   */

   // Check whether xmur and xmuf are positive and at least larger than 1.E-6
   if (xmur < 1.E-6 || xmuf < 1.E-6) {
      error<<"Selected scale factors too small ( < 1.E-6 )! Ignoring call."<<endl;
      return false;
   }

   // Check whether pQCD contributions beyond LO exist and are activated
   bool lknlo = false;
   if (!BBlocksSMCalc[kFixedOrder].empty()) {
      for (unsigned int i = 0 ; i <BBlocksSMCalc[kFixedOrder].size() ; i++) {
         int kOrder = BBlocksSMCalc[kFixedOrder][i]->GetIContrFlag2()-1;
         if (bUseSMCalc[kFixedOrder][i] && kOrder > 0) {
            lknlo = true;
            break;
         }
      }
   }

   // Check whether threshold corrections exist and are activated
   bool lkthc = false;
   if (!BBlocksSMCalc[kThresholdCorrection].empty()) {
      for (unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++) {
         if (bUseSMCalc[kThresholdCorrection][i]) {
            lkthc = true;
            break;
         }
      }
   }
   if (lkthc && fabs(xmur-xmuf) > DBL_MIN) {
      warn["SetScaleFactorsMuRMuF"]
         <<"Threshold corrections do not allow different scale factors for MuR and MuF, nothing changed!\n";
      warn["SetScaleFactorsMuRMuF"]
         <<"The method returns 'false', please check the return code and act appropriately.\n";
      man<<"Please do only symmetric scale variations, i.e. xmur = xmuf, with threshold corrections switched on\n";
      man<<"or deactivate threshold corrections first using\n";
      man<<"FastNLOReader::SetContributionON(kTresholdCorrections,Id,false).\n";
      return false;
   }

   // Deal with factorization scale first
   // Check whether corresponding xmuf variation exists in case of v2.0 table
   if (!GetIsFlexibleScaleTable() && !fUseHoppet) {
      const int ns = GetNScaleVariations();
      debug["SetScaleFactorsMuRMuF"]<<"Found "<<ns<<" scale variations for contributions switched ON."<<endl;
      int sfnlo = -1;
      if (lknlo) {
         for (int is = 0 ; is<ns ; is++) {
            if (fabs(((fastNLOCoeffAddFix*)B_NLO())->GetScaleFactor(is)-xmuf) < DBL_MIN) {
               sfnlo = is;
               break;
            }
         }
      }
      int sfthc = -1;
      if (lkthc) {
         for (int is = 0 ; is<ns ; is++) {
            if (fabs(((fastNLOCoeffAddFix*)B_ThC())->GetScaleFactor(is)-xmuf) < DBL_MIN) {
               sfthc = is;
               break;
            }
         }
      }
      if (lknlo && sfnlo == -1) {
         warn["SetScaleFactorsMuRMuF"]<<"Could not find NLO table with given mu_f scale factor of "<<xmuf<<", nothing changed!"<<endl;
         warn["SetScaleFactorsMuRMuF"]
            <<"The method returns 'false', please check the return code and act appropriately.\n";
         return false;
      }
      if (lkthc && sfthc == -1) {
         warn["SetScaleFactorsMuRMuF"]<<"Could not find ThC table with given mu_f scale factor of "<<xmuf<<", nothing changed!"<<endl;
         warn["SetScaleFactorsMuRMuF"]
            <<"The method returns 'false', please check the return code and act appropriately.\n";
         return false;
      }
      if (lkthc && lknlo && sfnlo != sfthc) {
         warn["SetScaleFactorsMuRMuF"]<<"Order of scale variation tables different in NLO and ThC tables, "<<sfnlo<<" != "<<sfthc<<" !"<<endl;
         warn["SetScaleFactorsMuRMuF"]<<"This is currently not supported, nothing changed!"<<endl;
         warn["SetScaleFactorsMuRMuF"]
            <<"The method returns 'false', please check the return code and act appropriately.\n";
         return false;
      }

      // Finally change renormalization scale first. Otherwise safety check in SetScaleVariationfails!
      fScaleFacMuR = xmur;
      // Now set factorization scale
      fScaleFacMuF = xmuf;
      bool bSetScales = false;
      if (lknlo) {
         bSetScales = SetScaleVariation(sfnlo);
         if (!bSetScales) {
            error["SetScaleFactorsMuRMuF"]<<"NLO scale variation table "<<sfnlo<<" could not be selected, stopped!"<<endl;
            exit(1);
         }
      } else { // LO only
         bSetScales = SetScaleVariation(0);
         if (!bSetScales) {
            error["SetScaleFactorsMuRMuF"]<<"LO scale variation table "<< 0 <<" could not be selected, stopped!"<<endl;
            exit(1);
         }
      }
      PrintScaleSettings();
   }
   else if (fUseHoppet) {

      debug["SetScaleFactorsMuRMuF"]<< "UseHoppet==true: Default factorization scale variation table (muF=1.0) will be enabled and "
         << "Hoppet will be used to calculate scale variation contributions on the fly." << endl;

      const int ns = GetNScaleVariations();
      debug["SetScaleFactorsMuRMuF"]<<"Found "<<ns<<" scale variations for contributions switched ON."<<endl;
      int sfnlo = -1;
      if (lknlo) {
         for (int is = 0 ; is<ns ; is++) {
            if (fabs(((fastNLOCoeffAddFix*)B_NLO())->GetScaleFactor(is) - 1.0) < DBL_MIN) {
               sfnlo = is;
               break;
            }
         }
      }

      bool bSetScales = false;
      if (lknlo) {
         bSetScales = SetScaleVariation(sfnlo);
         if (!bSetScales) {
            error["SetScaleFactorsMuRMuF"]<<"NLO scale variation table "<<sfnlo<<" could not be selected, stopped!"<<endl;
            exit(1);
         }
      }
      fScaleFacMuR = xmur;
      fScaleFacMuF = xmuf;
      PrintScaleSettings();
   }
   else {
      fScaleFacMuR = xmur;
      fScaleFacMuF = xmuf;
      PrintScaleSettings(kMuR);
      PrintScaleSettings(kMuF);
   }
   return true;
}


//______________________________________________________________________________
void fastNLOReader::PrintScaleSettings(fastNLO::EMuX MuX) {
   if (!GetIsFlexibleScaleTable()) {
      info["PrintScaleSettings"]<<"Renormalization scale chosen to be mu_r = "<<fScaleFacMuR<<" * "<<B_LO()->GetScaleDescription()<<endl;
      info["PrintScaleSettings"]<<"Factorization scale chosen to be   mu_f = "<<fScaleFacMuF<<" * "<<B_LO()->GetScaleDescription()<<endl;
   } else {
      // ---- prepare printout ---- //
      static const string sname[2] = {"Renormalization","Factorization"};
      static const string smu[2] = {"mu_r","  mu_f"};
      const int isc = MuX==kMuR?0:1;
      const double sfac = MuX==kMuR?fScaleFacMuR:fScaleFacMuF;
      EScaleFunctionalForm func = MuX==kMuR?fMuRFunc:fMuFFunc;
      char fname[100];
      switch (func) {
      case kScale1:
         sprintf(fname,"%s^2",B_LO()->GetScaleDescription(0).c_str());
         break;
      case kScale2:
         sprintf(fname,"%s^2",B_LO()->GetScaleDescription(1).c_str());
         break;
      case kQuadraticSum:
         sprintf(fname,"(%s^2 + %s^2)",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kQuadraticMean:
         sprintf(fname,"(%s^2 + %s^2)/2",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kQuadraticSumOver4:
         sprintf(fname,"(%s^2 + %s^2)/4",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kLinearMean:
         sprintf(fname,"((%s+%s)/2)^2",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kLinearSum:
         sprintf(fname,"(%s+%s)^2",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kScaleMax:
         sprintf(fname,"max(%s^2,%s^2)",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kScaleMin:
         sprintf(fname,"min(%s^2,%s^2)",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kProd:
         sprintf(fname,"(%s*%s)^2)",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kExpProd2:
         sprintf(fname,"(%s*exp(0.3*%s)^2)",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      case kExtern:
         sprintf(fname,"f_ext(%s,%s)",B_LO()->GetScaleDescription(0).c_str(),B_LO()->GetScaleDescription(1).c_str());
         break;
      default:
         error<<"unknown scale choice.\n";

      }
      info["PrintScaleSettings"]<<sname[isc]<<" scale chosen to be "<<smu[isc]<<"^2 = "
          <<sfac<<"^2 * "<<fname<<endl;
   }
}


//_____________________________________________________________________________
double fastNLOReader::CalcMu(fastNLO::EMuX kMuX , double scale1, double scale2, double scalefac) {
   //!
   //!  Calculate the scales with the defined function and the
   //!  corresponding prefactor.
   //!
   if (kMuX == kMuR && fScaleFacMuR != scalefac) error<<"Sth. went wrong with the scales.\n";
   if (kMuX == kMuF && fScaleFacMuF != scalefac) error<<"Sth. went wrong with the scales.\n";

   EScaleFunctionalForm Func = (kMuX == kMuR) ? fMuRFunc : fMuFFunc;
   double mu = 0;
   if (Func == fastNLO::kScale1)            mu      = scale1;
   else if (Func == fastNLO::kScale2)            mu      = scale2;
   else if (Func == fastNLO::kQuadraticSum)      mu      = FuncMixedOver1(scale1,scale2);
   else if (Func == fastNLO::kQuadraticMean)     mu      = FuncMixedOver2(scale1,scale2);
   else if (Func == fastNLO::kQuadraticSumOver4) mu      = FuncMixedOver4(scale1,scale2);
   else if (Func == fastNLO::kLinearMean)        mu      = FuncLinearMean(scale1,scale2);
   else if (Func == fastNLO::kLinearSum)         mu      = FuncLinearSum(scale1,scale2);
   else if (Func == fastNLO::kScaleMax)          mu      = FuncMax(scale1,scale2);
   else if (Func == fastNLO::kScaleMin)          mu      = FuncMin(scale1,scale2);
   else if (Func == fastNLO::kProd)              mu      = FuncProd(scale1,scale2);
   else if (Func == fastNLO::kExpProd2)          mu      = FuncExpProd2(scale1,scale2);
   else if (Func == fastNLO::kExtern)            mu      = (kMuX==kMuR) ? (*Fct_MuR)(scale1,scale2) : (*Fct_MuF)(scale1,scale2);
   else error["CalcMu"]<<"Could not identify functional form for scales calculation.\n";

   return scalefac * mu;
}


//______________________________________________________________________________
double fastNLOReader::FuncMixedOver1(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 1.));
}


//______________________________________________________________________________
double fastNLOReader::FuncMixedOver2(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 2.));
}


//______________________________________________________________________________
double fastNLOReader::FuncMixedOver4(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 4.));
}


//______________________________________________________________________________
double fastNLOReader::FuncLinearMean(double scale1 , double scale2) {
   return (scale1 + scale2) / 2.;
}


//______________________________________________________________________________
double fastNLOReader::FuncLinearSum(double scale1 , double scale2) {
   return scale1 + scale2;
}


//______________________________________________________________________________
double fastNLOReader::FuncMax(double scale1 , double scale2) {
   if (scale1 > scale2) return scale1;
   else return scale2;
}


//______________________________________________________________________________
double fastNLOReader::FuncMin(double scale1 , double scale2) {
   if (scale1 < scale2) return scale1;
   else return scale2;
}


//______________________________________________________________________________
double fastNLOReader::FuncProd(double scale1 , double scale2) {
   return (scale1 * scale2);
}


//______________________________________________________________________________
double fastNLOReader::FuncExpProd2(double scale1 , double scale2) {
   return (scale1 * exp(0.3*scale2));
}


//______________________________________________________________________________
int fastNLOReader::ContrId(const ESMCalculation eCalc, const ESMOrder eOrder) const {
   int Id = -1;
   if (BBlocksSMCalc.empty() || bUseSMCalc[eCalc].empty()) {
      return Id;
   }

   // Requested order
   string requested = _OrdName[eCalc][eOrder];
   // Loop over all available orders of contribution type eCalc
   for (unsigned int i=0; i<BBlocksSMCalc[eCalc].size(); i++) {
      // Found order
      int iFlag1 = BBlocksSMCalc[eCalc][i]->GetIContrFlag1();
      int iFlag2 = BBlocksSMCalc[eCalc][i]->GetIContrFlag2();
      string available = _OrdName[iFlag1-1][iFlag2-1];
      if (available == requested) {
         Id = i;
      }
   }
   return Id;
}



//______________________________________________________________________________
//
//                      Print outs
//______________________________________________________________________________



//______________________________________________________________________________
void fastNLOReader::PrintTableInfo(const int iprint) const {
   //! this function is inherited from fastNLOTable.
   fastNLOTable::PrintTableInfo(iprint);
}


//______________________________________________________________________________
void fastNLOReader::PrintFastNLOTableConstants(const int iprint) const {
   //! this function is inherited from fastNLOTable.
   this->fastNLOTable::PrintFastNLOTableConstants(iprint);
}


//______________________________________________________________________________
void fastNLOReader::PrintCrossSections() const {
   //!
   //! Print Cross sections in NLO, k-factors and Reference table cross sections
   //!

   //   if ( XSection.empty() )    CalcCrossSection();
   //   if ( XSectionRef.empty() && XSectionRef_s1.empty() )    CalcReferenceCrossSection();}

   vector < double > xs = XSection;

   printf(" #  \n");
   printf(" #  FastNLO Cross sections for\n");
   for (unsigned int i = 0 ; i < ScDescript.size() ; i++) {
      printf(" #     %s\n",ScDescript[i].c_str());
   }
   printf(" #  at sqrt(s) = %8.2f GeV\n", Ecms);
   printf(" #  \n");
   printf(" #  This is a %s-differential table in %s", ((NDim==1)?"single":"double"),GetDimLabel(0).c_str());
   if (NDim==2) printf(" and in %s",GetDimLabel(1).c_str());
   printf(".\n");
   printf(" #\n");

   string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
   string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
   string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;


   if (NDim == 2) {
      double lobindim2 = -42;
      printf(" #  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |\n",GetDimLabel(1).c_str(),unit[Ipublunits].c_str());
      printf(" #  --------------------------------------------------------------------\n");
      for (unsigned int i=0; i<xs.size(); i++) {
         if (GetObsBinLoBound(i,0) != lobindim2) {
            printf(" #                  ---->  from %9.3f to %9.3f in %s  <----\n",GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),GetDimLabel(0).c_str());
            lobindim2 = GetObsBinLoBound(i,0);
         }
         printf(" #   %4.0f   | %9.3f - %9.3f       % 9.4e           % 5.2f      |\n",i*1.,GetObsBinLoBound(i,1),GetObsBinUpBound(i,1),xs[i],kFactor[i]);
      }
   }

   else {
      printf("   ---  %5s  ---        - Bin -       -- XS-FNLO --  \n",GetDimLabel(NDim-1).c_str());
      for (unsigned int i=0; i<xs.size(); i++) {
         printf("  %9.3f - %9.3f   %3.0f         % 9.4e\n",GetObsBinLoBound(i,NDim-1),GetObsBinUpBound(i,NDim-1),i*1.,xs[i]);
      }
   }
   printf(" #  --------------------------------------------------------------------\n");
}


//______________________________________________________________________________
void fastNLOReader::PrintCrossSectionsWithReference() {
   //!
   //!  Print Cross sections in NLO, k-factors and Reference table cross sections
   //!
   //!  Please mention, that the reference cross section can be easily deviating
   //!  more than 20% (scales, pdfs, alpha_s, etc...). This does not mean that
   //!  the table is wrong!
   //!

   vector < double > xs = XSection;
   vector < double > xsref;
   if (XSection.empty())      CalcCrossSection();
   if (XSectionRef.empty() && XSectionRef_s1.empty())      CalcReferenceCrossSection();

   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc == kScale1 && fMuRFunc == kScale1)   {
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's1'\n");
         xsref = XSectionRef_s1;
      } else if (fMuFFunc == kScale2 && fMuRFunc == kScale2) {
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's2'\n");
         xsref = XSectionRef_s2;
      } else if (fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean) {
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
         xsref = XSectionRefMixed;
      } else {
         xsref = XSectionRefMixed;
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
      }
   } else xsref = XSectionRef;


   printf(" #  \n");
   printf(" #  FastNLO Cross sections for\n");
   for (unsigned int i = 0 ; i < ScDescript.size() ; i++) {
      printf(" #     %s\n",ScDescript[i].c_str());
   }
   printf(" #  at sqrt(s) = %8.2f GeV\n", Ecms);
   printf(" #  \n");
   printf(" #  This is a %s-differential table in %s", ((NDim==1)?"single":"double"),GetDimLabel(0).c_str());
   if (NDim==2) printf(" and %s",GetDimLabel(1).c_str());
   printf(" #  \n");
   printf(" #  Please mention, that the reference cross section can easily deviating up to more\n *  than 20%% due to different scale choices, alhpa_s value/evolution, PDFs, etc.");
   printf(" #  This does not mean, that this FastNLO table is wrong!\n\n");
   printf(" #  There are three reference cross sections stored for different scale choices.\n");
   printf(" #  If you have choosen mu_r=mu_f=%s, or mu_r=mu_f=%s or mu_r=mu_f=sqrt((%s^2+%s^2)/2), then you access automatically the corresponding reference cross section.\n",
          B_NLO()->GetScaleDescription(0).c_str(),B_NLO()->GetScaleDescription(1).c_str(),
          B_NLO()->GetScaleDescription(0).c_str(),B_NLO()->GetScaleDescription(1).c_str());
   printf(" #  In any other case your reference cross section is calculated using mu_r=mu_f=sqrt((%s^2+%s^2)/2).\n",
          B_NLO()->GetScaleDescription(0).c_str(),B_NLO()->GetScaleDescription(1).c_str());
   printf(" #  To be fully consistent with the nlojet++ reference cross section, you also have to adjust alpha_s and the alpha_s evolution accordingly.\n\n");

   printf("\n");
   printf(" #\n");

   string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
   string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
   string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;

   if (NDim == 2) {
      double lobindim2 = -321312;
      printf(" #  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |  -- XS-ref (NLOJET++) --    Diff [%%]\n",GetDimLabel(NDim-1).c_str(),unit[Ipublunits].c_str());
      printf(" #  -----------------------------------------------------------------------------------------------------------\n");
      for (unsigned int i=0; i<xs.size(); i++) {
         if (GetObsBinLoBound(i,0) != lobindim2) {
            printf(" #                    ---->  from %9.3f to %9.3f in %s  <----\n",GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),GetDimLabel(0).c_str());
            lobindim2 = GetObsBinLoBound(i,0);
         }
         printf(" #   %4.0f   | %9.3f - %9.3f      % 9.4e           % 5.3f      |     % 9.4e            % 5.4f\n",
                i*1.,GetObsBinLoBound(i,1),GetObsBinUpBound(i,1),xs[i],kFactor[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
      }
   }

   else {
      printf("FastNLOReader::PrintCrossSections( ). Info. Single differential printing of cross sections not yet nicely implemented.\n");
      printf("   ---  %s  ---        - Bin -    -- XS-FNLO  --       -- XS-ref (NLOJET++) --    Diff [%%]\n",GetDimLabel(NDim-1).c_str());
      for (unsigned int i=0; i<xs.size(); i++) {
         printf("  %9.3f - %9.3f   %3.0f         % 9.4e           % 9.4e          % 5.4f\n",GetObsBinLoBound(i,NDim-1),GetObsBinUpBound(i,NDim-1),i*1.,xs[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
      }
   }
   printf(" #  ------------------------------------------------------------------------------------------------------------\n");
}


//______________________________________________________________________________
void fastNLOReader::PrintCrossSectionsDefault(const vector <double> kthc) const {
   //!
   //! Print observable binnning and cross sections at
   //! LO, NLO and K factors like in Fortran Reader for comparison
   //!
   //! This function can only print double-differential crosss section tables.
   //!

   // Check on existence of 2-loop threshold corrections
   //const int ithc2 = kthc.empty() ? -1 : ContrId( fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
   const int ithc2 = kthc.empty() ? -1 : ContrId(kThresholdCorrection,kNextToLeading);

   cout << _DSEPLC << endl;
   printf(" Cross Sections\n");
   if (!GetIsFlexibleScaleTable())
      printf(" The scale chosen here are: mu_f = % #6.3f * %s, and mu_r = % #6.3f * %s \n",
             fScaleFacMuF, B_LO()->GetScaleDescription().c_str(), fScaleFacMuR, B_LO()->GetScaleDescription().c_str());
   cout << _SSEPLC << endl;

   if (NDim == 1) {
      // non-perturbative corrections (just first np correction)
      const int inpc1 = ContrId(kNonPerturbativeCorrection,kLeading);
      const vector < double > knpc = inpc1>-1 ? ((fastNLOCoeffMult*)BBlocksSMCalc[kNonPerturbativeCorrection][kLeading])->GetMultFactor() : vector<double>(NObsBin);

      string header0 = "  IObs  Bin Size IODimO   ";
      string header2 = "   LO cross section   NLO cross section   K NLO";
      if (ithc2>-1)header2 += "     K THC";
      if (inpc1>-1)header2 += "     K NPC";
      unsigned int NDimBins[NDim];
      printf("%s [ %-12s ] %s\n",
             header0.c_str(),DimLabel[0].c_str(),header2.c_str());
      cout << _SSEPLC << endl;
      for (unsigned int i=0; i<NObsBin; i++) {
         NDimBins[0] = 1;
         if (ithc2<0 && inpc1<0) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %#18.11E %#18.11E %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),
                   XSection_LO[i],XSection[i],kFactor[i]);
         } else if (inpc1<0 && ithc2 != -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),
                   XSection_LO[i],XSection[i],kFactor[i],kthc[i]);
         } else if (inpc1>-1 && ithc2 == -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),
                   XSection_LO[i],XSection[i],kFactor[i],knpc[i]);
         } else {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),
                   XSection_LO[i],XSection[i],kFactor[i],kthc[i],knpc[i]);
         }
         printf("\n");
      }
   } else if (NDim == 2) {
      // non-perturbative corrections (just first np correction)
      const int inpc1 = ContrId(kNonPerturbativeCorrection,kLeading);
      const vector < double > knpc = inpc1>-1 ? ((fastNLOCoeffMult*)BBlocksSMCalc[kNonPerturbativeCorrection][kLeading])->GetMultFactor() : vector<double>(NObsBin);

      string header0 = "  IObs  Bin Size IODimO ";
      string header1 = "   IODimI ";
      string header2 = " LO cross section   NLO cross section   K NLO";
      if (ithc2>-1)header2 += "     K THC";
      if (inpc1>-1)header2 += "     K NPC";
      unsigned int NDimBins[NDim];
      printf("%s [ %-12s ] %s [  %-12s  ] %s\n",
             header0.c_str(),DimLabel[1].c_str(),header1.c_str(),DimLabel[0].c_str(),header2.c_str());
      cout << _SSEPLC << endl;
      for (unsigned int i=0; i<NObsBin; i++) {
         for (unsigned int j=0; j<NDim; j++) {
            if (i==0)                                  NDimBins[j] = 1;
            else if (GetObsBinLoBound(i-1,j) < GetObsBinLoBound(i,j))       NDimBins[j]++;
            else if (GetObsBinLoBound(i,j) < GetObsBinLoBound(i-1,j))       NDimBins[j] = 1;
         }
         if (ithc2<0 && inpc1<0) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,1),GetObsBinUpBound(i,1),
                   NDimBins[1],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),XSection_LO[i],XSection[i],kFactor[i]);
         } else if (inpc1<0 && ithc2 != -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,1),GetObsBinUpBound(i,1),
                   NDimBins[1],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),XSection_LO[i],XSection[i],kFactor[i],kthc[i]);
         } else if (inpc1>-1 && ithc2 == -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,1),GetObsBinUpBound(i,1),
                   NDimBins[1],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),XSection_LO[i],XSection[i],kFactor[i],knpc[i]);
         } else {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetObsBinLoBound(i,1),GetObsBinUpBound(i,1),
                   NDimBins[1],GetObsBinLoBound(i,0),GetObsBinUpBound(i,0),XSection_LO[i],XSection[i],kFactor[i],kthc[i],knpc[i]);
         }
         printf("\n");
      }
   } else {
      warn["PrintCrossSectionsDefault"]<<"Print out optimized for two dimensions. No output for "<<NDim<<" dimensions."<<endl;
   }

}


//______________________________________________________________________________
void fastNLOReader::RunFastNLODemo() {
   //!
   //! This method prints out cross sections for different scale
   //! variation tables. Though it also changes the currently stored
   //! settings of this instance!
   //!
   //! RunFastNLODemo is changing settings (like scale choices) of this reader.
   //!

   info["RunFastNLODemo"]<<"RunFastNLODemo is changing settings (like scale choices) of this reader."<<endl;

   // If flexible-scale table, set MuR and MuF functional forms
   if (GetIsFlexibleScaleTable()) {
      SetMuRFunctionalForm(fastNLO::kScale1);
      SetMuFFunctionalForm(fastNLO::kScale1);
      //SetMuRFunctionalForm(fastNLO::kExpProd2);
      //SetMuRFunctionalForm(fastNLO::kExpProd2);
   }

   // Check on existence of LO and NLO (Id = -1 if not existing)
   int ilo   = ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
   int inlo  = ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
   if (ilo < 0 || inlo < 0) {
      error["PrintFastNLODemo"]<<"LO and/or NLO not found, nothing to be done!\n";
      return;
   }
   // Check on existence of 2-loop threshold corrections
   int ithc2 = ContrId(fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);

   // Pre-define desired order of scale variations
   const int nxmu = 4;
   double xmu[nxmu] = {1.0, 0.25, 0.5, 2.0};
   int   ixmu[nxmu] = { -1,   -1,  -1,  -1};
   // Get number of available scale variations and check on available scale factors,
   // in particular for MuF; set pointers
   int nscls = GetNScaleVariations();
   // With threshold corrections, allow only default scale (0)
   if (ithc2 > -1) {
      nscls = 1;
   }
   for (int iscls=0; iscls<nscls; iscls++) {
      SetScaleVariation(iscls);
      double fxmu = fScaleFacMuF;
      for (int i=0; i<nxmu; i++) {
         if (fabs(xmu[i]-fxmu) < 0.000001) {
            ixmu[i] = iscls;
         }
      }
   }

   // Loop over scales
   for (int iscls=0; iscls<nxmu; iscls++) {
      // First result is with NLO, LO result via division by K factor
      if (ixmu[iscls] > -1) {
         SetScaleVariation(ixmu[iscls]);
         CalcCrossSection();

         // Second result: Include threshold corrections for NLO if available
         vector < double > kthc;
         if (ithc2 > -1) {
            vector < double > stdk = kFactor;
            bool SetOn = SetContributionON(fastNLO::kThresholdCorrection, ithc2, true);
            if (!SetOn) {
               warn["RunFastNLODemo"]<<"Contribution "<<_ContrName[fastNLO::kThresholdCorrection]<<" could not be switched on for scale variation number "<<iscls<<", skipped!"<<endl;
               continue;
            }
            CalcCrossSection();
            kthc = kFactor;
            // Threshold K factor is NLO including 2-loop vs. NLO
            for (unsigned int i=0; i<kthc.size(); i++) {
               if (fabs(kFactor[i]) > DBL_MIN) {
                  kthc[i] = kFactor[i]/stdk[i];
               } else {
                  kthc[i] = -1.;
               }
            }
            SetContributionON(fastNLO::kThresholdCorrection, ithc2, false);
         }

         PrintCrossSectionsDefault(kthc);
      }
   }
}


//______________________________________________________________________________
double fastNLOReader::RescaleCrossSectionUnits(double binsize, int xunits) {
   //!
   //! This method rescales the stored cross section units according to
   //! the chosen Ipublunits and settings for [kAbsoluteUnits | kPublicationUNits].
   //!
   double unit = 1.;
   // For kAbsoluteUnits remove division by BinSize
   if (fUnits == kAbsoluteUnits) {
      unit *= binsize;
   }
   // Rescale SigmaTilde to Ipublunits (accounts for potentially different settings in XSectUnits)
   if ( xunits != Ipublunits ) {
      unit /= pow(10.,xunits-Ipublunits);
   }
   return unit;
}


//______________________________________________________________________________
fastNLOReader::XsUncertainty fastNLOReader::GetScaleUncertainty(const EScaleUncertaintyStyle eScaleUnc) {
   // Get 2- or 6-point scale uncertainty
   const double xmurs[7] = {1.0, 0.5, 2.0, 0.5, 1.0, 1.0, 2.0};
   const double xmufs[7] = {1.0, 0.5, 2.0, 1.0, 0.5, 2.0, 1.0};
   fastNLOReader::XsUncertainty XsUnc;

   unsigned int NObsBin = GetNObsBin();
   unsigned int npoint = 0;
   if ( eScaleUnc == kSymmetricTwoPoint ) {
      npoint = 2;
   } else if ( eScaleUnc == kAsymmetricSixPoint ) {
      npoint = 6;
   }

   debug["GetScaleUncertainty"]<<"npoint = "<<npoint<<endl;
   if ( npoint == 0 ) {
      info["GetScaleUncertainty"]<<"Only default scale selected, uncertainties will be zero."<<endl;
   } else if ( npoint == 2 ) {
      info["GetScaleUncertainty"]<<"Symmetric 2-point scale variations selected,"<<endl;
   } else if ( npoint == 6 ) {
      info["GetScaleUncertainty"]<<"Asymmetric 6-point scale variations selected,"<<endl;
   } else {
      error["GetScaleUncertainty"]<<"ERROR! No usual scale variation scheme selected, exiting."<<endl;
      error["GetScaleUncertainty"]<<"npoint = "<<npoint<<endl;
      exit(1);
   }

   //! Cross section and absolute uncertainties
   for ( unsigned int iscl = 0; iscl <= npoint; iscl++ ) {
      SetScaleFactorsMuRMuF(xmurs[iscl],xmufs[iscl]);
      CalcCrossSection();
      for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
         if ( iscl == 0 ) {
            XsUnc.xs.push_back(XSection[iobs]);
            XsUnc.dxsu.push_back(0);
            XsUnc.dxsl.push_back(0);
         } else {
            XsUnc.dxsu[iobs] = max(XsUnc.dxsu[iobs],XSection[iobs]-XsUnc.xs[iobs]);
            XsUnc.dxsl[iobs] = min(XsUnc.dxsl[iobs],XSection[iobs]-XsUnc.xs[iobs]);
         }
      }
   }

   //! Divide by cross section != 0 to give relative uncertainties
   for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
      if ( abs(XsUnc.xs[iobs]) > DBL_MIN ) {
         XsUnc.dxsu[iobs] = XsUnc.dxsu[iobs] / XsUnc.xs[iobs];
         XsUnc.dxsl[iobs] = XsUnc.dxsl[iobs] / XsUnc.xs[iobs];
      } else {
         XsUnc.dxsu[iobs] = 0.;
         XsUnc.dxsl[iobs] = 0.;
      }
      logger.debug["GetScaleUncertainty"]<<"iobs = " << iobs << "dxsl = " << XsUnc.dxsl[iobs] << ", dxsu = " << XsUnc.dxsu[iobs] <<endl;
   }

   logger.info["GetScaleUncertainty"]<<"Setting scale factors back to default of unity."<<endl;
   SetScaleFactorsMuRMuF(xmurs[0],xmufs[0]);

   return XsUnc;
}
