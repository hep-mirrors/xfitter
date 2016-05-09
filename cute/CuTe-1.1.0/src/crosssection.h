#ifndef CROSSSECTION_H
#define CROSSSECTION_H

#include "Hardfunc.h"
#include "resumed.h"
#include "fixed.h"


class DD_Kernel:public Pol<double>{
protected:
	CS_params* p;

	inline virtual void v_init(){}

	inline virtual void v_init(CS_params &p_){}

	inline virtual void v_update(){}

	inline DD_Kernel(CS_params &p_):p(&p_){}

public:
	inline void init(){
		//cout << "DD_Kernel::init()" << endl;
		v_init();
	}

	inline void init(double qT_, double mu_){
		p->update(qT_,mu_);
		init();
	}

	inline void init(CS_params &p_){
		p=(&p_);
		order(p->order);
		v_init(p_);
	}

	inline void update(){
		//cout << "DD_Kernel::update()" << endl;
		v_update();
	}

	inline void update(double y){
		p->update(y);
		update();
	}

	inline virtual ~DD_Kernel(){
	}
};



class DD_Kernel_LL:public DD_Kernel{
protected:
	_Hardfunc* H;
	K_Integral* K;
	BQB* B;

	inline void v_init(){
		H->update();
		K->update();
		//cout << "H_new=" << H << endl;
		//cout << "K_new=" << K << endl;
		B->init();
	}

	inline void v_init(CS_params &p_){
		delete H;
		if(p->prod==PROD_H) H = new Hardfunc_H(*p);
		else H = new Hardfunc_V(*p);

		delete K;
		if(p->res==RES_LL) K = new K_LL(*p);
		else K = new K_LLO(*p);

		delete B;
		if(p->prod!=PROD_W && p->_order<2) B = new BQB_fast(*p);
		else B = new BQB(*p);
	}

	inline void v_update(){
		//cout << "\tH_new\t=" << *H << endl;
		//cout << "\tK_new\t=" << K << endl;
		B->update();
		//cout << "\tE_new\t=" << K->E() << endl;
		//cout << "\tB_new\t=" << *B << endl;
		Pol<double> tmp(K->replaceL_E(*B));
		//cout << "\tKEB_new\t=" << tmp << endl;
		tmp*=(*H);
		for(uint i=0; i<=_order; ++i) _P[i]=tmp[i];
	}

public:

	inline DD_Kernel_LL(CS_params &p_):DD_Kernel(p_),H(NULL),K(NULL),B(NULL){
		init(p_);
	}

	inline ~DD_Kernel_LL(){
		delete H;
		delete K;
		delete B;
	}
};

class DD_Kernel_LO:public DD_Kernel {
protected:
	_FO* F;

	inline void v_init(){
		F->init();
	}

	inline void v_init(CS_params &p_){
		delete F;
		if(p->prod==PROD_H){
			F = new FO_H(*p);
		}else{
			F = new FO_DY(*p);
			//cout << "FO not implemented -> exit(1)" << endl;
			//exit(1);
		}
	}

	inline void v_update(){
		F->update();
		for(uint i=0; i<=_order; ++i){	_P[i]=(*F)[i];	}
	}

public:

	inline DD_Kernel_LO(CS_params &p_):DD_Kernel(p_),F(NULL){
		init(p_);
	}

	inline ~DD_Kernel_LO(){
		delete F;
	}

};

class Sigma_0{
protected:
	CS_params* p;
	double sigma_0;
	static constexpr double BARN = 0.389379323e9;

	inline virtual void v_update(){
		sigma_0=0.;
	}
public:
	inline void init(CS_params &p_){
		p=(&p_);
	}

	inline void update(){
		//cout << "Sigma_0::update()" << endl;
		//p->print(cout);
		v_update();
		sigma_0*=BARN*2.*p->qT;
		//cout << "\tsigma_0 = " << sigma_0 << endl;//
	}

	inline void update(double qT, double mu){
		p->update(qT,mu);
		update();
	}

	inline double operator ()(){return sigma_0;}

	inline Sigma_0(CS_params &p_):p(&p_){}

	inline virtual ~Sigma_0(){}
};

class Sigma_0_Z:public Sigma_0{
protected:
	inline void v_update(){
		sigma_0=4*PI2*ALPHA_EM/(N_c*p->s);
	}
public:
	inline Sigma_0_Z(CS_params &p_):Sigma_0(p_){}
};

class Sigma_0_gamma:public Sigma_0{
protected:
	inline void v_update(){
		sigma_0=4*M_PI/(3*N_c*p->s)*gsl_pow_2(ALPHA_EM/p->M);
	}
public:
	inline Sigma_0_gamma(CS_params &p_):Sigma_0(p_){}
};

class Sigma_0_H:public Sigma_0{
protected:
	inline complex<double> f_q(double x){
		if(x>=1.0) return(pow(asin(1.0/sqrt(x)),2));
		complex<double> f_q=sqrt(1-x);
		f_q=(f_q+1.)/(f_q-1.);
		f_q=log(f_q);
		return -0.25*(f_q-complex<double>(0.0,M_PI));
	}

	inline double A_q(){
		complex<double> A;
		double x;
		//double m_q[]  ={0.0,0.0,1.4,0.0,172.6,4.75};//FIXME 0-6
		/*for(int q=1; q<=6; ++q){
			//x=gsl_pow_2(2.0*getThreshold(p.thread,q)/p.M);
			x=gsl_pow_2(2.0*m_q[q]/p.M);
			std::cout << " x=" << x << std::endl;//FIXME
			//std::cout << "m_" << q << "=" << getThreshold(p.thread,q);
			//std::cout << "\tx_" << q << "=" << x;
			if(x>0){
				A+=(f_q(x)*(1.0-x)+1.0)*1.5*x;
				//std::cout << "\tA=" << A;
			}
			//std::cout << std::endl;
		}*/
		x=gsl_pow_2(2.0*172.6/p->M);
		A+=(f_q(x)*(1.0-x)+1.0)*1.5*x;
		A=A*conj(A);
		return (real(A));
	}

	inline void v_update(){
		//cout << "Sigma_0_H::v_update()" << endl;
		sigma_0=G_F*M_SQRT1_2*gsl_pow_2(4.*M_PI*p->a*p->M)/(288.*M_PI*p->s);
		//cout << "\tsigma_0 = " << sigma_0 << " * " << A_q() << " = ";
		sigma_0*=A_q();
	}
public:
	inline Sigma_0_H(CS_params &p_):Sigma_0(p_){}
};

class DD_CS{
public:
	CS_params p;
	vector<double> erg;
protected:
	Sigma_0* A;
	DD_Kernel* DD;


	inline virtual double v_erg(){
		DD->update();
		//cout << "\tDD_new = " << *DD << endl;
		return (*DD).von(1.);
	}

public:

	inline void init(CS_params &p_){
		//cout << "DD_CS::init(CS_params &p_)" << endl;
		p = CS_params(p_);
		//cout << p << p_ << endl;

		delete A;
		switch (p.prod) {
		case PROD_H: A = new Sigma_0_H(p); 	break;
		case PROD_PHOTON: A = new Sigma_0_gamma(p); 	break;
		case PROD_W: A = new Sigma_0_Z(p); 	break;
		case PROD_Z: A = new Sigma_0_Z(p); 	break;
		default:cout << "Wrong Production -> exit(1)" << endl; exit(1); break;
		}

		delete DD;
		if(p.res==RES_LO) DD = new DD_Kernel_LO(p);
		else DD= new DD_Kernel_LL(p);

		erg.resize(p.PDFn+1,0.);
	}

	inline DD_CS(CS_params &p_):A(NULL),DD(NULL){
		init(p_);
	}

	inline virtual ~DD_CS(){
		delete A;
		delete DD;
	}


	inline void update(double qT, double mu, double y=0.){
		//cout << "DiffCS::update("<<p.y<<")" << endl;
		p.update(qT,mu,y);
		//		LHAPDF::initPDF(0);
		A->update();

		DD->init();
		if(p.PDFn==0 || p.PDFerr!=PDF_ERR_GAUSSIAN){
		//if(p.PDFcenter || p.PDFn==0){
			erg[0]=v_erg()*(*A)();
		}

		for(int i=1; i<=p.PDFn; ++i){
			std::cout << "-" << flush;
			//LHAPDF::initPDF(i);
			erg[i]=v_erg()*(*A)();
		}
		//cout << "erg = " << erg << endl;
		//cout << "sigma_0 = " << (*A)() << endl;

	}
	inline void update(CS_params const &p_){
		update(p_.qT,p_.mu,p_.y);
	}

	inline double const operator [] (uint const i)const{ return erg[i];}

};



class D_CS:public DD_CS{
protected:

	inline static double integrand(double y, void* params){
		D_CS* DCS = (D_CS*) params;
		DCS->DD->update(y);
		return DCS->DD->von(1.);
	}

	inline double v_erg(){
		int mem = 20000;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (mem);
		double result, error;
		gsl_function k;

		k.function = &integrand;
		k.params = this;
		//cout << "-\tgsl_integration_qag (&k, 0., "<<p.ymax<<", "<<p.gen<<", "<<p.gen<<", "<<mem<<",2, w, &result, &error);" << endl;
		int status = gsl_integration_qag (&k, 0., p.y_max, p.gen, p.gen, mem,2, w, &result, &error);
		if(status!=GSL_SUCCESS) GSL_ERROR_CODE(status);
		gsl_integration_workspace_free (w);
		return 2.*result;
	}


public:

	inline D_CS(CS_params &p_):DD_CS(p_){
	}

	inline ~D_CS(){
	}
};



#endif
