#ifndef __fnlogeneratorconstants__
#define __fnlogeneratorconstants__

#include <string>
#include <vector>

namespace fastNLO {

   struct GeneratorConstants {
      //! GeneratorConstants
      //!
      //! Collection of generator specific constants.
      //! These are:
      //!  - name of generator
      //!  - references for generator
      //!  - (additional information about generator may 
      //!    be included in References)
      std::string Name; //!< Name of generator
      std::vector<std::string> References; //!< References for generator. Include additional information here (e.g. 'run-mode' or process)
      std::vector<std::string > GetCodeDescription() {
	 //! Get 'CodeDescription' usable for fastNLO table
	 std::vector<std::string > CodeDescr(References.size()+1);
	 CodeDescr[0] = Name;
	 for ( unsigned int i = 0 ; i<References.size() ; i++ ) 
	    CodeDescr[i+1] = References[i];
	 return CodeDescr;
      }
   };


   struct ProcessConstants {
      //! ProcessConstants
      //!
      //! Collection of process specific constants.
      //! Please see fastNLO table format definition for 
      //! a detailed explanation.
      //! 
      //! AsymmetricProcesses may only be required if NPDFDim=1 is chosen.
      int LeadingOrder; //!< Order in alpha_s of leading order process
      int UnitsOfCoefficients; //!< Units of coeffients as passed to fastNLO (negative power of 10: pb->12, fb->15
      int NPDF; //!< Number of PDFs involved
      int NSubProcessesLO; //!< Number of subprocesses of the considered process in LO run.
      int NSubProcessesNLO; //!< Number of subprocesses of the considered process in NLO run.
      int NSubProcessesNNLO; //!< Number of subprocesses of the considered process in NNLO run.
      int IPDFdef1; //!< Define PDF linear combinations corresponding to partonic subprocesses (hadron-hadron: 3
      int IPDFdef2; //!< Flag to define PDF linear combinations (dependent on IPDFdef1. Use 1 for jet-production in pp/ppbar)
      int IPDFdef3LO; //!< Unique identifier (dependent on IPDFdef1, IPDFdef2) for specifying the PDF linear combinations.
      int IPDFdef3NLO; //!< Unique identifier (dependent on IPDFdef1, IPDFdef2) for specifying the PDF linear combinations.
      int IPDFdef3NNLO; //!< Unique identifier (dependent on IPDFdef1, IPDFdef2) for specifying the PDF linear combinations. 
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffLO; //! PDF Linear combinations for LO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffNLO; //! PDF Linear combinations for NLO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffNNLO; //! PDF Linear combinations for NNLO calculation (used only if IPDFdef2==0)
      int NPDFDim; //!< Internal way to store PDF linear combinations. Use 1 (half-matrix storage) or 2 (full-matrix storage) for hadron-hadron collisions.
      std::vector<std::pair<int,int> > AsymmetricProcesses; //!< (if NPDFDim=1) Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax
      std::string Name; //!< Name of process (e.g. pp->2jets) (add 'run-mode' or more details)
      std::vector<std::string> References; //!< References for process (also other plain text lines can be included here)
      std::vector<std::string > GetProcessDescription() {
	 //! Get 'CodeDescription' usable for fastNLO table
	 std::vector<std::string > ProcDescr(References.size()+1);
	 ProcDescr[0] = Name;
	 for ( unsigned int i = 0 ; i<References.size() ; i++ ) 
	    ProcDescr[i+1] = References[i];
	 return ProcDescr;
      }
   };
      
};

#endif
