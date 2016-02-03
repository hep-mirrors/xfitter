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
			char Names[nExtraParamMax][32];
			double Value[nExtraParamMax], Step[nExtraParamMax], Min[nExtraParamMax], Max[nExtraParamMax], ConstrVal[nExtraParamMax], ConstrUnc[nExtraParamMax];
			int iExtraParamMinuit[nExtraParamMax], nExtraParam;
		};
		extern COMMON_ExtraPars_t extrapars_;

	// actual routine
	void getextraparsconstrchi2_(double& chi2) {

		int len=1024;
		char parname[len];
		double par;
		double unc;
		double bound_l;
		double bound_h;
		int status=0;

		chi2=0.0;
		for(int p=0; p<extrapars_.nExtraParam; p++) {
			if(extrapars_.ConstrUnc[p]==0.0) continue;
			mnpout_(extrapars_.iExtraParamMinuit+p, parname, &par, &unc, &bound_l, &bound_h, &status, len);
			if(status<0) {
				printf("ERROR in GetExtraParsConstrChi2: something is wrong with parameter %d\n", extrapars_.iExtraParamMinuit[p]);
				exit(1);
			}
			chi2+=pow((par-extrapars_.ConstrVal[p])/extrapars_.ConstrUnc[p], 2.0);
		}
		//printf("chi2: %e\n", chi2);
	}

	// actual routine
	void printminuitextrapars_() {

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
				printf("%5s%20s%9s%9s%9s%9s%9s%9s\n", "NO", "NAME", "VALUE", "STEP", "LIM L", "LIM H","SHIFT","REDUCT");
			}
			parname[std::string(parname).find(' ')]='\0';
			printf("%5d%20s%9.4f%9.4f", extrapars_.iExtraParamMinuit[p], parname, par, unc);
			if(bound_l!=0.0&&bound_h!=0.0) {
				printf("%9.4f%9.4f", bound_l, bound_h);
			}
			else {
				printf("%9s%9s", "", "");
			}
			if(extrapars_.ConstrUnc[p]==0.0) {
				printf("%9s%9s", "", "");
			}
			else {
				double shift=(par-extrapars_.ConstrVal[p])/extrapars_.ConstrUnc[p];
				double reduction=unc/extrapars_.ConstrUnc[p];
				printf("%9.4f%9.4f", shift, reduction);
			}
			printf("\n");
		}
	}
}

