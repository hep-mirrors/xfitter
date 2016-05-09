#ifndef FIXED_h
#define FIXED_h
#include "param.h"
#include <gsl/gsl_integration.h>

class _FO:public Pol<double>{
protected:
	CS_params* p;

	double qT;
	double y;
	double c[2];//[E^y,E^-y]

	double M2;
	double sqrts;
	double mp;
	int n_f;
	double gen;
	double a;
	double mu_F;

	double x[2];
	VectorPhi Phi[2];

	double shat;
	double that;
	double uhat;

	inline void shat_(double xb){
		shat=gsl_pow_2(sqrts)*x[0]*xb;
	}

	inline void that_(){
		that = (M2 - sqrts*mp*c[0]*x[0]);
	}

	inline void uhat_(double xb){
		uhat = (M2 - sqrts*mp*c[1]*xb);
	}

	inline void xb_(double t){
		x[1] = (t/(sqrts*mp*c[1]-gsl_pow_2(sqrts)*x[0]));
	}

	inline double x_min_(){
		return GSL_MAX(0,(sqrts*mp*c[1]-M2)
				/(gsl_pow_2(sqrts)-sqrts*mp*c[0]));
	}

	inline void updateparameter(double dx){
		x[0] = dx;
		that_();
		xb_(that);
		shat_(x[1]);
		uhat_(x[1]);
		Phi[0].update(x[0],mu_F);
		Phi[1].update(x[1],mu_F);
		//cout << "x=" << dx << "\tx2=" << x[1] << "\ts=" << shat << "\tt=" << that
		//		<< "\tu=" << uhat << "\txmin=" << x_min_() << endl;
	}

	inline virtual void v_update(){
		cout << "*************************************************************************" << endl
			 << "**                     class _FO::v_update not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}

	inline virtual void v_init(){
		cout << "*************************************************************************" << endl
			 << "**                     class _FO::v_init not defined -> exit(1)        **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}

public:

	inline void init(){
		mp = sqrt(M2+gsl_pow_2(p->qT));
		n_f = p->n_f;
		qT = p->qT;
		mu_F = p->mu;
		a = p->a;
		v_init();
	}

	inline void init(CS_params &p_){
		p=(&p_);
		order(p->order);
		M2 = gsl_pow_2(p->M);
		gen = p->gen;
		sqrts = sqrt(p->s);
		Phi[0] = (VectorPhi(p->sign[0]));
		Phi[1] = (VectorPhi(p->sign[1]));
	}

	inline void init(double qT_, double mu_){
		p->update(qT_,mu_);
		init();
	}

	inline void update(){
		y = p->y;
		c[0]=exp(p->y);
		c[1]=1./c[0];
		swap(c[0],c[1]);
		v_update();
		setX(a);
	}

	inline void update(double y_){
		p->update(y_);
		update();
	}

	inline _FO(CS_params &p_):p(&p_){}

	inline virtual ~_FO(){}
};

class FO_H: public _FO{
protected:

	/*inline void that_(){
		that = (M2 - sqrts*mp*c[1]*x[0]);
	}

	inline void uhat_(){
		uhat = (M2 - sqrts*mp*c[0]*x[1]);
	}

	inline void xb_(){
		x[1] = (-that/(gsl_pow_2(sqrts)*x[0]-sqrts*mp*c[0]));
	}

	inline double x_min_(){
		return GSL_MAX(0,(sqrts*mp*c[0]-M2)
				/(gsl_pow_2(sqrts)-sqrts*mp*c[1]));
	}*/

	inline double g_gg(){
		return N_c*(gsl_pow_4(M2) + gsl_pow_4(shat) + gsl_pow_4(that) + gsl_pow_4(uhat))/(shat*that*uhat);
	}

	inline double g_qg(){
		return C_F*(gsl_pow_2(that) + gsl_pow_2(shat))/(-uhat);
	}

	inline double g_gq(){
		return C_F*(gsl_pow_2(uhat) + gsl_pow_2(shat))/(-that);
	}

	inline double g_qbq(){
		return 2.*gsl_pow_2(C_F)*(gsl_pow_2(that) + gsl_pow_2(uhat))/(shat);
	}

	inline double f(){

		double result=0.;

		//------ GQ + QG--------
		result+=Phi[1][e_g]*g_qg()*Phi[0].sumofquarks(n_f);

		result+=Phi[0][e_g]*g_gq()*Phi[1].sumofquarks(n_f);

	#ifndef ONLY_MIXED
		//------ GG--------
		result+=Phi[0][e_g]*g_gg()*Phi[1][e_g];

		//------ QQ--------
		double tmp=0.;
		for(e_parton q=e_(n_f); q>e_g; --q){
			tmp+=Phi[0][q]*Phi[1][bar(q)];
			tmp+=Phi[0][bar(q)]*Phi[1][q];
		}
		result+=g_qbq()*tmp;
	#endif

		return result;
	}

	/*inline void updateparameter(double dx){
		x[0] = dx;
		that_();
		xb_();
		shat_(x[1]);
		uhat_();
		Phi[0].update(x[0],mu_F);
		Phi[1].update(x[1],mu_F);
	}*/


	inline static double integrand(double dx, void* parameter){
		FO_H* FO = (FO_H*)parameter;
		FO->updateparameter(dx);
		return (-FO->f()/(FO->x[0]*FO->that));
	}

	inline double FO_H_a(){
		double result, error;
		int status;


		double xmin=x_min_();
		if(xmin>1.) return 0.;

		int intMem=500000;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (intMem);
		gsl_function k;
		k.function = &integrand;
		k.params = this;


		//FIXMEstatus = gsl_integration_qags (&k, xmin, 1, gen/100., gen/100., intMem, w, &result, &error);
		status = gsl_integration_qag (&k, xmin, 1.0,  gen, gen, intMem, 3, w, &result, &error);


		if(status!=GSL_SUCCESS){
			cout << "----------- x Integral\t="<< result << " +- " << error << "\t D= "<< error/result*100 << "%" << endl
					<< "                      \tqT = " << p->qT << "\tmu = " << p->mu << endl;
			GSL_ERROR_CODE(status);
		}

		gsl_integration_workspace_free (w);

		return 2.*result/M2;
	}

	inline void v_update(){
		//swap(c[0],c[1]);
		if(_order){
			_P[1]=FO_H_a();
		}
	}

	inline void v_init(){}
public:

	inline FO_H(CS_params& p_):_FO(p_){
		init(p_);
	}
	inline ~FO_H(){}
};

class FO_DY:public _FO{
protected:
	Q_Matrix Q;
	Charge* C;

	inline double A_qq(double s, double t, double u){
		return (u*u+t*t+2.*s*M2)/(u*t);
	}

	inline double g_qq(){
		return 2.*C_F*A_qq(shat,that,uhat);
	}

	inline double g_qg(){
		return (-A_qq(uhat,that,shat));
	}

	inline double g_gq(){
		return (-A_qq(that,uhat,shat));
	}

	inline double jacobian(){
		return 1./(x[0]*sqrts*sqrts-sqrts*mp*c[1]);
	}

	inline double f(){
		double f=0.;
		double Lumi_GQ = Phi[0][e_g]*Phi[1].sumofQphi(Q,n_f);
		double Lumi_QG = Phi[0].sumofQphi(Q,n_f) *Phi[1][e_g];
		double Lumi_QQ = Phi[0]*Q*Phi[1];

#ifdef ONLY_MIXED
		Lumi_QQ=0.;
#endif
		f=Lumi_GQ*g_gq() + Lumi_QG*g_qg() + Lumi_QQ*g_qq();

		return f/(x[0]*x[1]);
	}

	inline static double integrand(double dx, void* parameter){
		FO_DY* FO = (FO_DY*)parameter;
		FO->updateparameter(dx);
		//cout << "f = " << FO->f() << "\tjac = " << FO->jacobian() << "\tdI = " << (FO->f()*FO->jacobian()) << endl;
		return (FO->f()*FO->jacobian());
	}

	inline double FO_DY_a(){
		double result, error;
		int status;


		double xmin=x_min_();
		if(xmin>1.) return 0.;

		int intMem=500000;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (intMem);
		gsl_function k;
		k.function = &integrand;
		k.params = this;


		//FIXMEstatus = gsl_integration_qags (&k, xmin, 1, gen/100., gen/100., intMem, w, &result, &error);
		status = gsl_integration_qag (&k, xmin, 1.0,  gen, gen, intMem, 3, w, &result, &error);


		if(status!=GSL_SUCCESS){
			cout << "----------- x Integral\t="<< result << " +- " << error << "\t D= "<< error/result*100 << "%" << endl
					<< "                      \tqT = " << p->qT << "\tmu = " << p->mu << endl;
			GSL_ERROR_CODE(status);
		}
		//cout << "Integral = " << result << endl;

		gsl_integration_workspace_free (w);

		return result;
	}

	inline void v_update(){
		if(_order){
			_P[1]=FO_DY_a();
		}
	}


	inline void v_init(){
		C->init(p->prod,p->ckm);
		Q = Q_Matrix(p->prod,p->ckm);
	}

public:
	inline FO_DY(CS_params& p_):_FO(p_),C(NULL){
		C = new Charge(p_.prod,p->ckm);
		init(p_);
	}
	inline ~FO_DY(){
		delete C;
	}
};








#endif
