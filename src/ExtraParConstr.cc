#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>

extern "C" { 
	
	// interface to MINUIT
	void mnpout_(int*,char*,double*,double*,double*,double*,int*,int); 

		#include "dimensions.h"
		const int nExtraParamMax=NEXTRAPARAMMAX_C;
		struct COMMON_ExtraPars_t {
			char Names[nExtraParamMax][128];
			double Value[nExtraParamMax], Step[nExtraParamMax], Min[nExtraParamMax], Max[nExtraParamMax], ConstrVal[nExtraParamMax], ConstrUnc[nExtraParamMax];
			int iExtraParamMinuit[nExtraParamMax], nExtraParam;
		};
		extern COMMON_ExtraPars_t extrapars_;

	// actual routine
	void getextraparsconstrchi2_(double& chi2) {
		int    len = 1024;
		char   parname[len];
		double par, unc, bound_l, bound_h;
		int    status = 0;

		chi2 = 0.0;
		for (int p = 0; p < extrapars_.nExtraParam; p++) {
		if (extrapars_.ConstrUnc[p] == 0.0) continue;

		mnpout_(extrapars_.iExtraParamMinuit + p,
				parname, &par, &unc, &bound_l, &bound_h, &status, len);
		if (status < 0) {
			std::printf("ERROR in GetExtraParsConstrChi2: parameter %d\n",
						extrapars_.iExtraParamMinuit[p]);
			std::exit(1);
		}
		const double deviation = (par - extrapars_.ConstrVal[p]) / extrapars_.ConstrUnc[p];
		chi2 += deviation * deviation;
		}
	}

	// actual routine
	void printminuitextrapars_(int& iflag) {

		int len=1024;
		char parname[len];
		double par;
		double unc;
		double bound_l;
		double bound_h;
		int status=0;

		int firstpar=1;
		for(int p=0; p<extrapars_.nExtraParam; p++) {
			mnpout_(extrapars_.iExtraParamMinuit+p, parname, &par, &unc, &bound_l, &bound_h, &status, len);		        


			if(status<=0) continue;

			if(firstpar) {
				firstpar=0;
				//printf("MINUIT EXTRA PARAMETERS:\n");
				printf("%5s%32s%9s%9s%9s%9s%9s%9s\n", "NO", "NAME", "VALUE", "STEP", "LIM L", "LIM H","SHIFT","REDUCT");
			}
			parname[std::string(parname).find(' ')]='\0';
			printf("%5d%32s%9.4f%9.4f", extrapars_.iExtraParamMinuit[p], parname, par, unc);
			if(bound_l!=0.0 || bound_h!=0.0) {
				printf("%9.4f%9.4f", bound_l, bound_h);
			}
			else {
				printf("%9s%9s", "", "");
			}
			if(extrapars_.ConstrUnc[p]==0.0) {
				printf("%9s%9s", "", "");
			} else {
				// SHIFT is the SIGNED pull (par-val)/unc, not its square: the sign
				// carries information and older outputs printed it this way.
				const double shift     = (par - extrapars_.ConstrVal[p]) / extrapars_.ConstrUnc[p];
				const double reduction = unc / extrapars_.ConstrUnc[p];
				std::printf("%9.4f%9.4f", shift, reduction);
			}
			printf("\n");
		}
    if(iflag == 3)
    {
      const int ndigcomma = 8;
      printf("----- Parameters in YAML format (can copy-paste into parameters.yaml):\n");
      for(int p=0; p<extrapars_.nExtraParam; p++)
      {
        mnpout_(extrapars_.iExtraParamMinuit+p, parname, &par, &unc, &bound_l, &bound_h, &status, len);
        //if(status<=0) continue;
        parname[std::string(parname).find(' ')]='\0';
        printf("  %s : [ %.*f, %.*f", parname, ndigcomma, par, ndigcomma, unc);
        if(bound_l!=0.0&&bound_h!=0.0) {
          printf(", %.*f, %.*f", ndigcomma, bound_l, ndigcomma, bound_h);
        }
        if (extrapars_.ConstrUnc[p] != 0.0) {
          // signed pull, matching the SHIFT column of the table above
          const double shift     = (par - extrapars_.ConstrVal[p]) / extrapars_.ConstrUnc[p];
          const double reduction = unc / extrapars_.ConstrUnc[p];
          std::printf(", %.*f, %.*f", ndigcomma, shift, ndigcomma, reduction);
        }
        printf(" ]\n");
      }
      printf("----- End of parameters in YAML format\n");
    }
    fflush(stdout);
	}
}

