#include "plot.h"

std::ostream& operator<<( std::ostream& os, const Plot& P )
{
	// Save a triplet to os
	//os << P.p;

	os << "q_T\t: " << P.q << endl;
	os << "mu\t: " << P.mu << endl;
	os << "lambda\t: " << P.lambda << endl;
	return os;
}

int Plot::createPlot(){
	ofstream of;
	if(fail){
		of.open(outfileName.c_str());
		of << FL << endl;
		of << fail_string << endl;
		of.close();
	}
	for(uint i=0; i<lambda.size(); ++i){
		p.Lambda=lambda[i];
		_setOutfileName("_l" + NumberToString(i));
		of.open(outfileName.c_str());
		of << FL << endl;
		of.close();
		PP = new PlotPoint(p,y.size(),matchingScheme,toPlot);
		PP->values[e_l]=p.Lambda;
		_selectqT();
	}
	if(lambda.size()==0){
		p.Lambda=0.;
		of.open(outfileName.c_str());
		of << FL << endl;
		of.close();
		PP = new PlotPoint(p,y.size(),matchingScheme,toPlot);
		PP->values[e_l]=p.Lambda;
		_selectqT();
	}


	return 1;
}

int Plot::readin(string configFileName){
	fail = false;
	fail_string.clear();
	infileName = configFileName;
	_setOutfileName();
	cout << "\t--\tOutfileName: " << outfileName << endl;


	string tmps;
	config = ConfigFile(infileName);


	if(!config.readInto(tmps,"collision")){
		cout << "\t--\tNo \"collision\" argument [PP; PPbar; PbarPbar] -> exit(1)" << endl;
		fail_string +="\t--\tNo \"collision\" argument [PP; PPbar; PbarPbar] -> exit(1)\n";
		fail = true;//exit(1);
	}
	if(tmps=="PP"){	p.col=COL_PP;
	}else{
		if(tmps=="PPbar"){	p.col=COL_PbarP;
		}else{
			if(tmps=="PbarPbar"){	p.col=COL_PbarPbar;
				}else{
					cout << "\t--\tWrong \"collision\" argument [PP; PPbar; PbarPbar] -> exit(1)" << endl;
					fail_string +="\t--\tWrong \"collision\" argument [PP; PPbar; PbarPbar] -> exit(1)\n";
					fail = true;//exit(1);
	}	}	}

	if(!config.readInto(tmps,"cutoff")){
		cout << "\t--\tNo \"cutoff\" argument [NO; HARD; DIPOLE; GAUSS]-> No cutoff" << endl;
		p.cutoff=CUTOFF_NO;
		//lambda.resize(1);
	}else{
		if(tmps=="NO"){
			p.cutoff=CUTOFF_NO;
			if(config.readInto(lambda,"LambdaNP")){
				cout << "\t--\t\"cutoff\" argument = NO" << endl
					 << "\t--\t-> \"LambdaNP\" argument has no effect" << endl;
				//lambda.clear();
			}
			//lambda.resize(1);
		}else{
			if(tmps=="HARD"){p.cutoff=CUTOFF_HARD;
			}else{
				if(tmps=="DIPOLE"){p.cutoff=CUTOFF_DIPOL;
				}else{
					if(tmps=="GAUSS"){p.cutoff=CUTOFF_GAUSS;
					}else{
						cout << "\t--\tWrong \"cutoff\" argument [NO; HARD; DIPOLE; GAUSS] -> exit(1)" << endl;
						fail_string +="\t--\tWrong \"cutoff\" argument [NO; HARD; DIPOLE; GAUSS] -> exit(1)\n";
						fail = true;//exit(1);
	}	}	}	}	}


	if(p.cutoff!=CUTOFF_NO && !config.readInto(lambda,"LambdaNP")){
		cout << "\t--\tNo \"LambdaNP\" argument [double] for cutoff -> exit(1)" << endl;
		fail_string +="\t--\tNo \"LambdaNP\" argument [double] for cutoff -> exit(1)\n";
		fail = true;//exit(1);
	}


	//if(!config.readInto(PDFErr,"PDFerror",false)){
	//		cout << "\t--\tNo \"PDFerror\" argument [true; false]-> No PDFerror" << endl;
	//}
	if(!config.readInto(tmps,"PDFerror")){
		cout << "\t--\tNo \"PDFerror\" argument [NO; HESSIAN; GAUSSIAN]-> NO" << endl;
		p.PDFerr=PDF_ERR_NO;
	}else{
		p.PDFerr=PDF_ERR_WRONG_STATEMENT;
		if(tmps=="NO") p.PDFerr=PDF_ERR_NO;
		if(tmps=="HESSIAN") p.PDFerr=PDF_ERR_HESSIAN;
		if(tmps=="GAUSSIAN") p.PDFerr=PDF_ERR_GAUSSIAN;
		if(p.PDFerr==PDF_ERR_WRONG_STATEMENT){
			cout << "\t--\tWrong \"PDFerror\" argument [NO; HESSIAN; GAUSSIAN] -> exit(1)" << endl;
			fail_string +="\t--\tWrong \"PDFerror\" argument [NO; HESSIAN; GAUSSIAN] -> exit(1)\n";
			fail = true;//exit(1);
		}
	}

	if(!config.readInto(muErr,"ScaleError",false)){
		cout << "\t--\tNo \"ScaleErr\" argument [true; false]-> ScaleErr=false" << endl;
	}

	if(!config.readInto(p.exponent, "exponent",false)){
		cout << "\t--\tNo \"exponent\" argument [true; false]-> exponent=false" << endl;
	}

	if(!config.readInto(p.gen, "accuracy",1e-4)){
		cout << "\t--\tNo \"accuracy\" argument [double]-> accuracy=1e-3" << endl;
	}

	if(!config.readInto(p.improved, "improved",true)){
		cout << "\t--\tNo \"improved\" argument [true; false]-> improved=true" << endl;
	}

	if(!config.readInto(p._order,"order",1)){
		cout << "\t--\tNo \"order\" argument [int]-> order=1" << endl;
	}

	if(!config.readInto(p.piquadrat,"PiQuadrat",true)){
		cout << "\t--\tNo \"PiQuadrat\" argument [true; false]-> PiQuadrat=true" << endl;
	}

	if(!config.readInto(tmps, "production")){
		cout << "\t--\tNo \"production\" argument [DY; H; Z, W]" << endl;
		fail_string +="\t--\tNo \"production\" argument [DY; H; Z, W]\n";
		fail = true;//exit(1);
	}else{
		p.C_B = C_F;
		if(tmps=="DY"){
			p.prod=PROD_PHOTON;
		}else{
			if(tmps=="H"){
				p.C_B = C_A;
				p.prod=PROD_H;
			}else{
				if(tmps=="Z"){
					p.prod=PROD_Z;
				}else{
					if(tmps=="W"){
						p.prod=PROD_W;//cout << "The Boson you have called is temporary not available -> exit(1)" << endl;
						//exit(1);
					}else{
						cout << "\t--\tWrong \"production\" argument [DY; H; Z, W] -> exit(1)" << endl;
						fail_string +="\t--\tWrong \"production\" argument [DY; H; Z, W] -> exit(1)\n";
						fail = true;//exit(1);
	}	}	}	}	}

	if(!config.readInto(q,"q_T")){
		cout << "\t--\tNo \"q_T\" argument [double ... double] -> exit(1)" << endl;
		fail_string +="\t--\tNo \"q_T\" argument [double ... double] -> exit(1)\n";
		fail = true;//exit(1);
	}

	if(!config.readInto(p.s,"sqrts")){
		cout << "\t--\tNo \"sqrts\" argument [double] -> exit(1)" << endl;
		fail_string +="\t--\tNo \"sqrts\" argument [double] -> exit(1)\n";
		fail = true;//exit(1);
	}else{ p.s*=p.s;}

	if(!config.readInto(p.M,"M")){
		cout << "\t--\tNo \"M\" argument [double] -> exit(1)" << endl;
		fail_string +="\t--\tNo \"M\" argument [double] -> exit(1)\n";
		fail = true;//exit(1);
	}

	if(!config.readInto(pdfName, "PDFset")){
			cout << "\t--\tNo \"PDFset\" argument [string] -> exit(1)" << endl;
			fail_string +="\t--\tNo \"PDFset\" argument [string] -> exit(1)\n";
			fail = true;//exit(1);
		}

	//	if(!config.readInto(p.PDFcenter, "PDFcenter")){
			//cout << "\t--\tNo \"PDFcenter\" argument [true; false] -> exit(1)" << endl;
		//	exit(1);
	//	}

	//Needs M,Prod,gen, PDFset
	if(!config.readInto(mu,"mu")){
		cout << "\t--\tNo \"mu\" argument [double ... double] -> mu = qT + q_star" << endl;
		if(!config.readInto(qstar, "q_star")){
			qstar = NAN;
			cout	<< "\t--\tNo \"q_star\" argument [double] ->  default q_star(M) = " << endl;
		}
	}else if(config.readInto(qstar, "q_star"))
		cout << "\t--\t\"mu\" argument given -> \"q_star\" argument has no effect" << endl;





	if(!config.readInto(y,"y")){
		cout << "\t--\tNo \"y\" argument -> ds/dq" << endl;
	}else{
		cout << "\t--\t\"y\" argument found -> d2s/(dqdy)" << endl;
	}

	if(!config.readInto(tmps,"resummation")){
		cout << "\t--\tNo \"resummation\" argument [RES, RES_FO, FO, M_CORR, MATCHED] -> MATCHED" << endl;
		p.res= RES_MATCHED;
	}else{
		p.res = RES_WRONG_STATEMENT;
		if(tmps=="RES")	p.res= RES_LL;
		if(tmps=="RES_FO")	p.res=RES_LLO;
		if(tmps=="FO")	p.res=RES_LO;
		if(tmps=="M_CORR")	p.res=RES_MATCHING_CORR;
		if(tmps=="MATCHED")	p.res=RES_MATCHED;
		if(p.res==RES_WRONG_STATEMENT){
			cout << "\t--\tWrong \"resummation\" argument [RES, RES_FO, FO, M_CORR, MATCHED] -> exit(1)" << endl;
			fail_string +="\t--\tWrong \"resummation\" argument [RES, RES_FO, FO, M_CORR, MATCHED] -> exit(1)\n";
			fail = true;//exit(1);
	}	}



	if(!config.readInto(p.mu_h,"mu_h")){
		cout << "\t--\tNo \"mu_h\" argument -> mu_h=M" << endl;
		p.mu_h=p.M;
	}else{
		if(p.res!=RES_LL){
			cout << "\t--\tmu_h does not effect fixed-order calculations" << endl;
		}
	}

	if(!config.readInto(p.mu_t,"mu_t") && p.prod==PROD_H){
		cout << "\t--\tNo \"mu_t\" argument -> mu_t=mu_h" << endl;
		p.mu_t=p.mu_h;
	}

	vector<double> tmpckm;
	if(config.readInto(tmpckm,"CKM")){
		if(tmpckm.size()!=9){
			cout << "\t--\tWRONG \"CKM\" argument size double[9] -> exit(1)" << endl;
			fail_string +="\t--\tWRONG \"CKM\" argument size double[9] -> exit(1)\n";
			fail = true;//exit(1);
		}
		p.ckm[lha_u][lha_d]=tmpckm[0];
		p.ckm[lha_u][lha_s]=tmpckm[1];
		p.ckm[lha_u][lha_b]=tmpckm[2];

		p.ckm[lha_c][lha_d]=tmpckm[3];
		p.ckm[lha_c][lha_s]=tmpckm[4];
		p.ckm[lha_c][lha_b]=tmpckm[5];

		p.ckm[lha_t][lha_d]=tmpckm[6];
		p.ckm[lha_t][lha_s]=tmpckm[7];
		p.ckm[lha_t][lha_b]=tmpckm[8];
	}else{
	  /*
		p.ckm[lha_u][lha_d]=0.97427;
		p.ckm[lha_u][lha_s]=0.22534;
		p.ckm[lha_u][lha_b]=0.00351;

		p.ckm[lha_c][lha_d]=0.22520;
		p.ckm[lha_c][lha_s]=0.97344;
		p.ckm[lha_c][lha_b]=0.0412;

		p.ckm[lha_t][lha_d]=0.00867;
		p.ckm[lha_t][lha_s]=0.0404;
		p.ckm[lha_t][lha_b]=0.999146;
	  */
	        p.ckm[lha_u][lha_d]=ckm_matrix_.Vud;
	        p.ckm[lha_u][lha_s]=ckm_matrix_.Vus;
	        p.ckm[lha_u][lha_b]=ckm_matrix_.Vub;

		p.ckm[lha_c][lha_d]=ckm_matrix_.Vcd;
		p.ckm[lha_c][lha_s]=ckm_matrix_.Vcs;
		p.ckm[lha_c][lha_b]=ckm_matrix_.Vcb;

		p.ckm[lha_t][lha_d]=ckm_matrix_.Vtd;
		p.ckm[lha_t][lha_s]=ckm_matrix_.Vts;
		p.ckm[lha_t][lha_b]=ckm_matrix_.Vtb;
		if(p.prod==PROD_W){
			cout	<< "\t--\tNo \"CKM\" argument -> CKM =" << endl
					<< "\t--\t\t(" << p.ckm[lha_u][lha_d] << ",\t"
					         << p.ckm[lha_u][lha_s] << ",\t"
					         << p.ckm[lha_u][lha_b] << ")" << endl
					<< "\t--\t\t(" << p.ckm[lha_c][lha_d] << ",\t"
					         << p.ckm[lha_c][lha_s] << ",\t"
					         << p.ckm[lha_c][lha_b] << ")" << endl
					<< "\t--\t\t(" << p.ckm[lha_t][lha_d] << ",\t"
					         << p.ckm[lha_t][lha_s] << ",\t"
					         << p.ckm[lha_t][lha_b] << ")" << endl ;
		}
	}



	if(!config.readInto(tmps, "MScheme")){
		if(p.res==RES_MATCHED || p.res==RES_MATCHING_CORR)
		cout << "\t--\tNo \"MScheme\" argument [NAIVE; CONST; VAR] -> NAIVE " << endl;
		matchingScheme = MATCHING_NAIV;
	}else{
		if(tmps=="NAIVE") matchingScheme = MATCHING_NAIV;
		else 	if(tmps=="CONST") matchingScheme = MATCHING_CONST;
		else 	if(tmps=="VAR") matchingScheme = MATCHING_VAR;
		else 	if(tmps=="CONST0") matchingScheme = MATCHING_CONST_NULL;
		else 	if(tmps=="VAR0") matchingScheme = MATCHING_VAR_NULL;
		else {
			cout << "\t--\tWRONG \"MScheme\" argument [NAIVE; CONST; VAR] -> exit(1)" << endl;
			fail_string +="\t--\tWRONG \"MScheme\" argument [NAIVE; CONST; VAR] -> exit(1)\n";
			fail = true;//exit(1);
		}
	}

	return 1;
}
int Plot::init(string configFileName){

	lambda.clear();
	q.clear();
	mu.clear();
	y.clear();
	//values = vector<double>(e_Last+1,0);
	toPlot = vector<int>(e_Last+1,0);
	PDFErr=false;
	muErr=false;





	cout << "Plot::init -> Plot::readin(" <<configFileName<<")"<< endl;
	Plot::readin(configFileName);

	//	LHAPDF::initPDFSet(pdfName, LHAPDF::LHGRID);
	LHAinit=true;

	if(PDFErr){
	  p.PDFn = 0;//LHAPDF::numberPDF();
	}
	else p.PDFn=0;

	p.init(p._order,p.improved,p.exponent,p.piquadrat,p.col,p.prod,p.res,p.ckm,p.gen,p.PDFerr,p.PDFn,p.cutoff,p.Lambda);
	p.init(p.s,p.M,p.mu_t,p.mu_h);
	_FirstLine();

	if (qstar!=qstar){
		qstar = p.q_star();
		cout << "qstar = " << qstar;
	}


	return 1;
}
