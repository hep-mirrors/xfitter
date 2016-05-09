#ifndef PLOT_H
#define PLOT_H

#include <string>
#include <vector>
#include <fstream>
#include "ConfigFile.h"
#include "crosssection.h"
#include <ctime>
#include <gsl/gsl_mode.h>


using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;


enum TOPLOT_param {e_qT, e_y, e_l, e_mu};
enum TOPLOT_err {e_val, e_muP, e_muM, e_Delta, e_DeltaP, e_DeltaM};
enum TOPLOT_val{
	e_erg = e_mu+1,
	e_LL = e_erg+e_DeltaM+1,
	e_LLO = e_LL+e_DeltaM+1,
	e_LO = e_LLO+e_DeltaM+1
};

enum TOPLOT_Last {e_Last = e_LO+e_DeltaM};

enum MATCHING{ MATCHING_NAIV=0, MATCHING_CONST, MATCHING_VAR, MATCHING_CONST_NULL, MATCHING_VAR_NULL};


template <class T> T max(std::vector<T>& Vec){
	T tmp = Vec[0];
	for(uint i=1; i<Vec.size(); ++i) tmp= GSL_MAX(tmp,Vec[i]);
	return tmp;
}
template <class T> T min(std::vector<T>& Vec){
	T tmp = Vec[0];
	for(uint i=1; i<Vec.size(); ++i) tmp= GSL_MIN(tmp,Vec[i]);
	return tmp;
}


template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}


class PlotPoint{
public:
	vector<double> values;

protected:
	CS_params* p;
	Hardfunc Hf;
	DD_CS* D_LL;
	DD_CS* D_LLO;
	DD_CS* D_LO;


	bool dy;
	MATCHING mscheme;

	vector<double> ergs;
	vector<int> toPlot;


	vector<vector<double> > muerrors;
	vector<double> muval;

	inline void init_muerr(){
		//cout <<"PlotPoint::init_muerr()" << endl;
		muval = vector<double>(2*nMuErr+1 , 0.0);
		muval[0]=1.;
		for(int i= 1;i<=nMuErr; ++i){
			muval[i]=1.+(i*(MuFac-1.))/((double)nMuErr);
			muval[2*nMuErr+1-i]=1.+(1.-MuFac)/(MuFac*nMuErr)*i;
		}
		//cout << "\tmuval\t= " << muval << endl;
		muerrors.resize(4);
		for(uint i=0; i<muerrors.size(); ++i) muerrors[i] = vector<double>(2*nMuErr+1 , 0.0);
		//cout << "\tmuerrors\t= " << muerrors << endl;
	}


	inline void standard_deviation(TOPLOT_val res, vector<double> &erg){
		double S0 =0.0;
		double norm = 1./(erg.size()-1.);
		for(uint i=1; i<erg.size(); ++i){
			S0+=erg[i];
		}
		S0*=norm;
		erg[0]=S0;
		//cout	<< " norm = " << norm << "\tsize = " << erg.size()
			//	<< "\t S0 = " << S0 << endl;

		//cout << values[res+e_Delta] << " -> ";
		for(uint i=1; i<erg.size(); ++i){
			//cout << "+= " << gsl_pow_2(S0-erg[i]) << endl;
			values[res+e_Delta]+=gsl_pow_2(S0-erg[i]);
		}
		//cout << values[res+e_Delta] << " -> ";
		values[res+e_Delta]=sqrt(norm*values[res+e_Delta]);
		//cout << values[res+e_Delta] << endl;
		//p.error[1]=p.error[2]=0.;
	}

	inline void hessian_errors(TOPLOT_val res, vector<double> &erg){
		double S0= erg[0];
		for(int i=1; i<=(0.5*(erg.size()-1.)); ++i){
			//cout << "i=" << i << endl;
			values[res+e_Delta]+=gsl_pow_2(erg[2*i-1]-erg[2*i]);
			values[res+e_DeltaP]+=gsl_pow_2(gsl_max(erg[2*i-1]-S0,gsl_max(erg[2*i]-S0,0.0)));
			values[res+e_DeltaM]+=gsl_pow_2(gsl_max(S0-erg[2*i-1],gsl_max(S0-erg[2*i],0.0)));
		}
		values[res+e_Delta]=0.5*sqrt(values[res+e_Delta]);
		values[res+e_DeltaP]=sqrt(values[res+e_DeltaP]);
		values[res+e_DeltaM]=sqrt(values[res+e_DeltaM]);
	}

	inline void _calcPDFErr(TOPLOT_val res, vector<double> &erg){
		values[res+e_Delta]=values[res+e_DeltaP]=values[res+e_DeltaM]=0.0;
		if(p->PDFerr==PDF_ERR_HESSIAN) hessian_errors(res,erg);
		else standard_deviation(res,erg);
	}


	inline void _calcRes(TOPLOT_val res, vector<double> &erg, double sign=+1.){
		//cout << "PlotPoint::_calcRes("<<res<<","<<erg<<", "<<sign<<")" << endl;
		if(p->PDFn) _calcPDFErr(res,erg);
		for(int i=0; i<=p->PDFn; i++) ergs[i]+=sign*erg[i];
		values[res]=erg[0];
		/*cout << "\t ergs\t= " << ergs << endl;
		cout << "\t values\t= " << values << endl;
		cout << "\t muerrors\t= " << muerrors << endl;*/
	}

	inline void _calcMuErr(TOPLOT_val res, DD_CS& CS){
		CS.p.PDFn=0;
		//muerrors[RES_MATCHED][0]+=muerrors[p.resummation][0];
		muerrors[CS.p.res][0]=CS.erg[0];
		for(int i=1; i<=2*nMuErr; ++i){
			//tep.update(p->qT,p->mu*muval[i],p->y);
			CS.update(p->qT,p->mu*muval[i],p->y);
			//cout << "CS("<<CS.p.mu << ") = " << CS[0] << endl;
			muerrors[CS.p.res][i]=CS[0];
		}
		CS.p.PDFn=p->PDFn;
		CS.p.update(p->qT,p->mu,p->y);
		//cout << muerrors << endl;
		values[res+e_muP]=max(muerrors[CS.p.res])-muerrors[CS.p.res][0];
		values[res+e_muM]=muerrors[CS.p.res][0]-min(muerrors[CS.p.res]);
	}

public:
	inline void init(CS_params &p_){
		//cout << "->PlotPoint::init(CS_params &p_)" << p_ << endl;
		p=(&p_);
		//cout << p << endl;
		init_muerr();
		values.resize(toPlot.size());
		ergs.resize(p->PDFn+1);
		delete D_LL;
		delete D_LLO;
		delete D_LO;



		CS_params tmpParam;
		if(p->res == RES_LL || p->res==RES_MATCHED){
		//if(toPlot[e_LL]){
			CUTOFF tmp = p->cutoff;
			if(p->Lambda==0.) tmp=CUTOFF_NO;
			tmpParam.init(p->_order,p->improved,p->exponent,p->piquadrat,p->col,p->prod,RES_LL,p->ckm,p->gen,p->PDFerr,p->PDFn, tmp ,p->Lambda);
			tmpParam.init(p->s,p->M,p->mu_t,p->mu_h);
			//cout << "p=" << p << endl << "p_LL = " << p_LL << endl;
			if(dy)	D_LL = new DD_CS(tmpParam);
			else D_LL = new D_CS(tmpParam);
		}
		if(p->res != RES_LL && p->res!=RES_LO){
		//if(toPlot[e_LLO]){
			tmpParam.init(p->_order,false,false,false,p->col,p->prod,RES_LLO,p->ckm,p->gen,p->PDFerr,p->PDFn);
			tmpParam.init(p->s,p->M,p->mu_t,p->mu);
			if(dy)	D_LLO = new DD_CS(tmpParam);
			else D_LLO = new D_CS(tmpParam);
		}
		if(p->res != RES_LL && p->res!=RES_LLO){
		//if(toPlot[e_LO]){
			tmpParam.init(p->_order,false,false,false,p->col,p->prod,RES_LO,p->ckm,p->gen,p->PDFerr,p->PDFn);
			tmpParam.init(p->s,p->M,p->mu_t,p->mu);
			if(dy)	D_LO = new DD_CS(tmpParam);
			else D_LO = new D_CS(tmpParam);
		}
		if(mscheme){
			uint n = p->_order;
			if(mscheme==MATCHING_CONST_NULL || mscheme==MATCHING_VAR_NULL) n=0;
			tmpParam.init(n,p->improved,p->exponent,p->piquadrat,p->col,p->prod,RES_LL,p->ckm,p->gen,p->PDFerr,p->PDFn, p->cutoff ,p->Lambda);
			tmpParam.init(p->s,p->M,p->mu_t,p->mu_h);
			//cout << "PlotPoint::init create Hf(" << (void*)&tmpParam << endl;
			Hf.init(tmpParam);
		}
		//cout << "<-PlotPoint::init(CS_params &p_)" << endl;
	}

	/*inline PlotPoint():p(NULL)
	,p_LL(NULL),p_LLO(NULL),p_LO(NULL),Hf(NULL),D_LL(NULL),D_LLO(NULL),D_LO(NULL){
	}*/

	inline PlotPoint(CS_params &p_, bool dy_, MATCHING mscheme_, vector<int>& toPlot_)
	:Hf(Hardfunc(p_)),D_LL(NULL),D_LLO(NULL),D_LO(NULL)
	,dy(dy_),mscheme(mscheme_),toPlot(toPlot_){
		init(p_);
	}

	inline ~PlotPoint(){
		delete D_LL;
		delete D_LLO;
		delete D_LO;
	}

	inline void update(){
		//cout << "->PlotPoint::update()" << endl;
		ergs.assign(p->PDFn+1,0.);
		if(p->res == RES_LL || p->res==RES_MATCHED){
		//if(toPlot[e_LL]){//cout << "PlotPoint::update()\t p_LL=" << endl;
			//p_LL->update(*p);
			D_LL->update(*p);
			_calcRes(e_LL,D_LL->erg);
			if(toPlot[e_erg+e_muP]) _calcMuErr(e_LL,*D_LL);
		}
		if(p->res != RES_LL && p->res!=RES_LO){
			//if(toPlot[e_LLO]){//cout << "PlotPoint::update()\t p_LLO=" << endl;
			//p_LLO->update(*p);
			D_LLO->update(*p);
			_calcRes(e_LLO,D_LLO->erg,-1);
			if(toPlot[e_erg+e_muP]) _calcMuErr(e_LLO,*D_LLO);
		}
		if(p->res != RES_LL && p->res!=RES_LLO){
			//if(toPlot[e_LO]){//cout << "PlotPoint::update()\t p_LO=" << endl;
			//p_LO->update(*p);
			D_LO->update(*p);
			_calcRes(e_LO,D_LO->erg);
			if(toPlot[e_erg+e_muP]) _calcMuErr(e_LO,*D_LO);
		}

		double H=1.;
		if(mscheme) H = Hf(p->qT,p->mu).von(1.);
		//cout << H;
		values[e_erg] = values[e_LL]+H*(values[e_LO]-values[e_LLO]);
		for(int i=0; i<(2*nMuErr+1); ++i){
			if(mscheme==MATCHING_VAR || mscheme==MATCHING_VAR_NULL) H = Hf(p->qT,p->mu*muval[i]).von(1.);
			//cout << "\t" << H;
			muerrors[RES_MATCHED][i]=muerrors[RES_LL][i]+H*(muerrors[RES_LO][i]-muerrors[RES_LLO][i]);
		}
		//cout << endl;
		values[e_erg+e_muP]=max(muerrors[RES_MATCHED])-muerrors[RES_MATCHED][0];
		values[e_erg+e_muM]=muerrors[RES_MATCHED][0]-min(muerrors[RES_MATCHED]);
		if(p->PDFn) _calcPDFErr(e_erg,ergs);
	}
};


class Plot{
public:
	ConfigFile config;

	CS_params p;

	vector<double> lambda;
	vector<double> q;
	vector<double> mu;
	vector<double> y;
	vector<int> toPlot;

	PlotPoint* PP;

	bool PDFErr;
	bool muErr;
	double qstar;
	bool LHAinit;
	MATCHING matchingScheme;

	clock_t start, tmp, end;

	string infileName;
	string outfileName;
	string pdfName;
	string FL;


	bool fail;
	string fail_string;

	Plot(){LHAinit=false;fail=false;};
	~Plot(){
		delete PP;
	}
	int init(string configFileName);
	int readin(string configFileName);
	Plot(string configFileName){LHAinit=false;init(configFileName);};
	int createPlot();


	void print(){
		cout << p;

		cout << "q_T\t: " << q << endl;
		cout << "mu\t: " << mu << endl;
		cout << "lambda\t: " << lambda << endl;
	}

	inline int _setOutfileName(string suffix){
		size_t found = infileName.find_last_of(".");
		outfileName = infileName.substr(0,found) + suffix + ".out";
		//cout << "outfileName: \""<< outfileName << "\"" << endl;
		return 1;
	}

	inline int _setOutfileName(){
		string tmp;
		return _setOutfileName(tmp);
	}

	inline void plot(){
		p.update(p.qT,p.mu,p.y);
		PP->update();

		ofstream of;
		of.open(outfileName.c_str(), ofstream::app);
		int n = static_cast<int>(log10(10./p.gen));
		of.width(n+2);
		of.precision(n);
		for(int i=0; i<=e_Last; ++i){
			if(toPlot[i]) of << PP->values[i] << "\t";
		}
		of << endl;
		of.close();
	}


	inline void _selecty(){
		if(!y.size()) plot();
		for(uint i=0; i<y.size(); ++i){
			PP->values[e_y]=p.y=y[i];
			plot();
		}
	}

	inline void _selectmu(){
		for(uint i=0; i<mu.size(); ++i){
			PP->values[e_mu]=p.mu=mu[i];
			_selecty();
		}
		if(mu.size()==0){
			PP->values[e_mu]=p.mu=p.qT+qstar;
			_selecty();
		}
	}

	inline void _selectqT(){
		start = tmp = clock();
		for(uint i=0; i<q.size(); ++i){
			PP->values[e_qT]=p.qT=q[i];
			_selectmu();
			end =clock();
			cout <<  endl << "qT = " << q[i] << endl << "\t-- Ticks =" << (end-tmp)*1e-6 << "m\t= " << (end-tmp)/CLOCKS_PER_SEC << "sec" << endl;
			tmp = end;
		}
		end = clock();
		cout << "Ticks =" << (end-start)*1e-6 << "m\t= " << (end-start)/CLOCKS_PER_SEC << "sec" << endl;
	}


	inline void print_standard_in(){
		ofstream of;
		of.open("example.in");
		of		<< "#\texample infile" << endl
				<< "#\t<- comments" << endl
				<< "#\tAll mass afflicted parameters in GeV" << endl
				<< "" << endl
				<< "#--------- Process -----------" << endl
				<< "collision\t=\tPP\t\t\t#[PP; PPbar; PbarPbar]" << endl
				<< "production\t=\tH\t\t\t#[DY; H; Z; W]" << endl
				<< "M\t\t=\t125\t\t\t#[double] GeV" << endl
				<< "sqrts\t\t=\t8000\t\t \t#[double] GeV" << endl
				<< "#ckm\t\t=\tVud Vus Vub ... Vtb\t#[double ... double]" << endl
				<< "#--------- Plot -----------" << endl
				<< "q_T\t\t=\t0.1 0.2 0.3 0.4 0.6 0.8 #[double ... double] GeV" << endl
				<< "\t\t\t1. 1.5 2. 3. 4. 5. 5.5" << endl
				<< "\t\t\t6. 6.5 7. 7.5 8. 8.5 9." << endl
				<< "\t\t\t9.5 10. 10.5 11. 11.5" << endl
				<< "\t\t\t12. 12.5 13. 13.5 14." << endl
				<< "\t\t\t14.5 15. 16. 17. 18. 19." << endl
				<< "\t\t\t20. 22.5 25. 27.5 30." << endl
				<< "\t\t\t32.5 35. 40. 45. 50." << endl
				<< "\t\t\t55. 60.\t" << endl
				<< "#y\t\t=\t-1 -0.5 0. 0.5 1\t#[double ... double]" << endl
				<< "accuracy\t=\t1e-3\t\t\t#[double]" << endl
				<< "#--------- Scales -----------" << endl
				<< "#mu\t\t=\t2 3 4 5 6 7 8 9 10 \t#[double ... double] GeV" << endl
				<< "#q_star\t\t=\t7.7\t\t\t#[double] GeV" << endl
				<< "#mu_h\t\t=\t125\t\t\t#[double] GeV" << endl
				<< "mu_t\t\t=\t172.6\t\t\t#[double] GeV" << endl
				<< "#ScaleError\t=\tfalse\t\t\t#[true; false]" << endl
				<< "#--------- Resummation ----------" << endl
				<< "#resummation\t=\tMATCHED\t\t\t#[RES, RES_FO, FO, M_CORR, MATCHED]" << endl
				<< "order\t\t=\t1\t\t\t#[int]" << endl
				<< "improved\t=\ttrue\t\t\t#[true; false]" << endl
				<< "PiQuadrat\t=\ttrue\t\t\t#[true; false]" << endl
				<< "#exponent\t=\tfalse\t\t\t#[true; false]" << endl
				<< "#MScheme\t=\tNAIVE\t\t\t#[NAIVE; CONST; VAR]" << endl
				<< "#--------- Hadronic effects ----------" << endl
				<< "#cutoff\t\t=\tGAUSS\t\t\t#[NO; HARD; DIPOLE; GAUSS]" << endl
				<< "#LambdaNP\t=\t0.\t0.3\t0.6\t#[double ... double] GeV" << endl
				<< "#--------- PDFs ----------" << endl
				<< "PDFset\t\t=\tMSTW2008nnlo68cl \t#[string]" << endl
				<< "#PDFerror\t=\tNO\t\t\t#[NO; HESSIAN; GAUSSIAN]" << endl
			<< "" << endl;

		of.close();
	}


	inline string _FirstLine(){


		time_t now = time(0);
		char* dt = ctime(&now);

		FL.clear();
		FL += "######################################################\n";
		FL += "##  ";											FL += "\n";
		FL += "##  ";					FL += PACKAGE_STRING;	FL += "\n";
		FL += "##  ";					FL += dt;
		FL += "##  ";											FL += "\n";
		FL += "######################################################\n";
		FL += "##  ";											FL += "\n";
		FL += "##  Infile: ";			FL += infileName;		FL += "\n";
		FL += "##  ";											FL += "\n";
		FL += "######################################################\n";


		std::ifstream in;
		in.open(infileName.c_str());
		string tmp;
		while(!getline(in,tmp).eof()){
			FL += "##\t";
			FL += tmp;
			FL += "\n";
		}
		in.close();
		FL += "######################################################\n";
		FL += "##  Outfile: ";									FL += "\n";
		FL += "######################################################\n";

		//string
		FL += "# ";


		for(int i=0; i<=e_Last; ++i){
			if(i==e_qT && q.size()){
				toPlot[e_qT]=1;
				FL +=  "q_T\t";
			}
			if(i==e_y && y.size()){
				toPlot[e_y]=1;
				FL +=  "y\t";
			}
			if(i==e_l && lambda.size()){
				toPlot[e_l]=1;
				FL +=  "l\t";
			}
			if(i==e_mu){
				toPlot[e_mu]=1;
				FL +=  "mu\t";
			}
			if(i==e_erg){
				toPlot[e_erg]=1;
				if(y.size()) FL +=  "d2s/dqdy\t";
				else         FL +=  "ds/dq\t";
			}

			if(i==(e_erg+e_Delta) && PDFErr){
				toPlot[(e_erg+e_Delta)]=1;
				 FL += "D PDF\t";
			}
			if(i==(e_erg+e_DeltaP) && toPlot[(e_erg+e_Delta)] && p.PDFerr==PDF_ERR_HESSIAN){
				toPlot[(e_erg+e_DeltaP)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_erg+e_DeltaM) && toPlot[(e_erg+e_DeltaP)]){
				toPlot[(e_erg+e_DeltaM)]=1;
				 FL += "D PDF-\t";
			}
			//if(i==(e_erg+e_muP) && nMuErr){
			if(i==(e_erg+e_muP) && muErr){
				toPlot[e_erg+e_muP]=1;
				 FL += "D mu+\t";
			}
			if(i==(e_erg+e_muM) && toPlot[e_erg+e_muP]){
				toPlot[e_erg+e_muM]=1;
				 FL += "D mu-\t";
			}

			if(i==e_LL && (//p.res==RES_LL ||
							p.res==RES_MATCHED)){
				toPlot[e_LL]=1;
				//for(uint n=0; n<p.order; n++) FL += "N";
				 FL += "RES\t";
			}

			if(i==(e_LL+e_Delta) && toPlot[(e_erg+e_Delta)] && toPlot[e_LL]){
				toPlot[(e_LL+e_Delta)]=1;
				 FL += "D PDF\t";
			}
			if(i==(e_LL+e_DeltaP) && toPlot[(e_erg+e_DeltaP)] && toPlot[e_LL]){
				toPlot[(e_LL+e_DeltaP)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_LL+e_DeltaM) && toPlot[(e_LL+e_DeltaP)]){
				toPlot[(e_LL+e_DeltaM)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_LL+e_muP) && toPlot[e_erg+e_muP] && toPlot[e_LL]){
				toPlot[(e_LL+e_muP)]=1;
				 FL += "D mu+\t";
			}
			if(i==(e_LL+e_muP) && toPlot[(e_LL+e_muP)]){
				toPlot[(e_LL+e_muM)]=1;
				 FL += "D mu-\t";
			}

			if(i==e_LLO && //!(p.res==RES_LL || p.res==RES_LO)
							(p.res==RES_MATCHED || p.res==RES_MATCHING_CORR)
																){
				toPlot[e_LLO]=1;
				//for(int n=0; n<p._order; n++) FL += "N";
				 FL += "RES_FO\t";
			}

			if(i==(e_LLO+e_Delta) && toPlot[(e_erg+e_Delta)] && toPlot[e_LLO]){
				toPlot[(e_LLO+e_Delta)]=1;
				 FL += "D PDF\t";
			}
			if(i==(e_LLO+e_DeltaP) && toPlot[(e_erg+e_DeltaP)] && toPlot[e_LLO]){
				toPlot[(e_LLO+e_DeltaP)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_LLO+e_DeltaM) && toPlot[(e_LLO+e_DeltaP)]){
				toPlot[(e_LLO+e_DeltaM)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_LLO+e_muP) && toPlot[e_erg+e_muP] && toPlot[e_LLO]){
				toPlot[(e_LLO+e_muP)]=1;
				 FL += "D mu+\t";
			}
			if(i==(e_LLO+e_muP) && toPlot[(e_LLO+e_muP)]){
				toPlot[(e_LLO+e_muM)]=1;
				 FL += "D mu-\t";
			}

			if(i==e_LO && //!(p.res==RES_LL || p.res==RES_LLO)
					(p.res==RES_MATCHED || p.res==RES_MATCHING_CORR)
														){
				toPlot[e_LO]=1;
				//for(int n=0; n<p._order; n++) FL += "N";
				 FL += "FO\t";
			}

			if(i==(e_LO+e_Delta) && toPlot[(e_erg+e_Delta)] && toPlot[e_LO]){
				toPlot[(e_LO+e_Delta)]=1;
				 FL += "D PDF\t";
			}
			if(i==(e_LO+e_DeltaP) && toPlot[(e_erg+e_DeltaP)] && toPlot[e_LO]){
				toPlot[(e_LO+e_DeltaP)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_LO+e_DeltaM) && toPlot[(e_LO+e_DeltaP)]){
				toPlot[(e_LO+e_DeltaM)]=1;
				 FL += "D PDF+\t";
			}
			if(i==(e_LO+e_muP) && toPlot[e_erg+e_muP] && toPlot[e_LO]){
				toPlot[(e_LO+e_muP)]=1;
				 FL += "D mu+\t";
			}
			if(i==(e_LO+e_muP) && toPlot[(e_LO+e_muP)]){
				toPlot[(e_LO+e_muM)]=1;
				 FL += "D mu-\t";
			}
		}

		return FL;
	}


};

extern std::ostream& operator<<( std::ostream& os, const Plot& P );

#endif
