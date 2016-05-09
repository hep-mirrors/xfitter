#ifndef HARDFUNC_H
#define HARDFUNC_H


#include "param.h"
//#define FIXEDRGE


class U_exp: public Pol< complex<double> >{
protected:
	A_strong a_strong;
	static constexpr double epsilon=1e-6;

	inline Pol<double> Gamma(int n_f){
		Pol<double> Gammapol(3,0.);
		Gammapol[0]=Gamma_0;
		Gammapol[1]=Gamma_1(n_f);
		Gammapol[2]=Gamma_2(n_f);
		Gammapol[3]=Gamma_3(n_f);
		return Gammapol;
	}

	inline Pol<double> gammaS(int n_f){
		Pol<double> gammaH(2,0.);
		gammaH[0]=gamma_0S;
		gammaH[1]=gamma_1S(n_f);
		gammaH[2]=gamma_2S(n_f);
		return gammaH;
	}

	inline Pol<double> gammat(int n_f){
		Pol<double> gammaH(2,0.);
		gammaH[0]=0.0;
		gammaH[1]=-2.*beta_1(n_f);
		gammaH[2]=-2.*2.*beta_2(n_f);
		return gammaH;
	}

	inline Pol<double> gammaq(int n_f){
		Pol<double> gammaq(2,0.);
		gammaq[0]=gamma_0q;
		gammaq[1]=gamma_1q(n_f);
		gammaq[2]=gamma_2q(n_f);
		return gammaq;
	}

	inline Pol< complex<double> > a_gamma(complex<double> a1, complex<double> a2, Pol< double > gamma, uint n_f){
		Pol< complex<double> > erg(3,0.);
		erg[0]= log(a2/a1)*gamma(0);
		if(_order>0) erg[1]=(a2-a1)*(gamma(1)-gamma(0)*beta_1_0(n_f));
		if(_order>1) erg[2]=(a2*a2-a1*a1)*0.5*(gamma(2)
							-gamma(0)*beta_2_0(n_f)
							-beta_1_0(n_f)*(gamma(1)-gamma(0)*beta_1_0(n_f)));
		erg.order(_order);
		erg*=0.5*beta__0(n_f);
		return erg;
	}

	inline Pol< complex<double> > S(complex<double> a_1, complex<double> a_2, int n_f){
		complex<double> r=a_2/a_1;
		complex<double> l=log(r);
		Pol< complex<double> > S(2,0.);
		S[0]=(r-1.0-r*l)/(r*a_1);
		S[0]+=(Gamma_1_0(n_f)-beta_1_0(n_f))*(-r+l+1.)+0.5*beta_1_0(n_f)*l*l;
		S[1]=(a_1*((beta_1_0(n_f)*Gamma_1_0(n_f)-beta_2_0(n_f))*(1.-r+r*l)+(gsl_pow_2(beta_1_0(n_f))-beta_2_0(n_f))*(1.-r)*l
					-0.5*(gsl_pow_2(beta_1_0(n_f))-beta_2_0(n_f)-beta_1_0(n_f)*Gamma_1_0(n_f)+Gamma_2_0(n_f))*(1.-r*2.+r*r)));
		S*=(0.25*Gamma_0*gsl_pow_2(beta__0(n_f)));
		S[2]=(-1.833142541410085*pow(a_2,2) + (2.506993010267918*pow(a_2,3))/a_1 -
				   0.7322252658394445*pow(a_1,2) + 0.5084142962680396*pow(a_2,2)*l -
				   4.421483041233234*pow(a_1,2)*l + 0.05837479698161122*pow(a_2,2)*r);
				//(-0.058374796981611354 - 3.9130687449651957*l + 0.058374796981611354*r)*pow(a_2,2);
		S.order(_order);
		/*cout << "\t\tS_new("<<a_1<<", "<<a_2<<", "<<n_f<<")"
			 << "\t=" << S << endl;
		Spol(a_1,a_2.real(),n_f,_order);*/
		return S;
	}

	inline virtual Pol< complex<double> > u_exp(complex<double> a_h, complex<double> a_s, complex<double> L_h, uint n_f){
		cout << "*************************************************************************" << endl
			 << "**                      class U_exp::u_exp not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
		return Pol< complex<double> >();
	}

public:

	/*inline void update(double mu_h, complex<double> a_h, int n_f_h, double mu_s, complex<double> a_s, int n_f_s, complex<double> L_h){

	}*/
	inline void update(complex<double> M, complex<double> mu_h, complex<double> mu){
		for(uint i=0;i<=_order; ++i) _P[i]=0.;

		int n_f_h = get_n_f( mu_h);
		complex<double> a_h = a_strong(mu_h);
		complex<double> L_h = 2.*log(M/mu_h);

		int n_f_s = get_n_f( mu);
		complex<double> a_s = a_strong(mu);

		double a;
		for(lha_parton q = n_f_h; q>n_f_s; --q){
			double mu=get_threshold(q);
			a=a_strong(mu+epsilon);
			*this+=u_exp(a_h,a,L_h,q);
			//a_h=a_von(mu-epsilon,alpha_s);
			//cout << "U("<<a_h<<","<<a<<","<<L_h<<","<<q<<") = " << u_exp(a_h,a,L_h,q)*2. << endl;
			a_h=a_strong(mu-epsilon);
			L_h=2.*log(M/mu);
		}
		//a_s = a_von(mu_s,alpha_s);
		//a_s = a_strong(mu_s);
		//cout << "U("<<a_h<<","<<a_s<<","<<L_h<<","<<n_f_s<<") = " << u_exp(a_h,a_s,L_h,n_f_s)*2. << endl;
		*this+=u_exp(a_h,a_s,L_h,n_f_s);
	}

	inline U_exp(uint order_,A_strong &a_strong_)
	:a_strong(a_strong_){
		order(order_);
	}

	/*inline U_exp(uint order_, Alpha_s const &alpha_s_){init(order_,alpha_s_);}
	inline U_exp(uint order_, double mu0, double alpha0, uint order0, double gen, uint n_f){
		init(order_, mu0, alpha0, order0, gen, n_f);}
	inline U_exp(_AllParams const &p){init(p);}*/

	virtual ~U_exp(){}
};

class U_exp_V:public U_exp{
protected:
	inline Pol< complex<double> > u_exp(complex<double> a_h, complex<double> a_s, complex<double> L_h, uint n_f){
		return ( S(a_h,a_s,n_f)*2.*C_F
				 -a_gamma(a_h,a_s,gammaq(n_f),n_f)*2.
				 -a_gamma(a_h,a_s,Gamma(n_f),n_f)*L_h*C_F );
	}
public:
	inline U_exp_V(uint order_,A_strong &a_strong_):U_exp(order_,a_strong_){}
};

class U_exp_S:public U_exp{
protected:
	inline Pol< complex<double> > u_exp(complex<double> a_h, complex<double> a_s, complex<double> L_h, uint n_f){
		return (  S(a_h,a_s,n_f)*2.*C_A
				-a_gamma(a_h,a_s,gammaS(n_f),n_f)
				 -a_gamma(a_h,a_s,Gamma(n_f),n_f)*L_h*C_A );
	}
public:
	inline U_exp_S(uint order_,A_strong &a_strong_):U_exp(order_,a_strong_){}
};

class U_exp_t:public U_exp{
protected:
	inline Pol< complex<double> > u_exp(complex<double> a_h, complex<double> a_s, complex<double> L_h, uint n_f){
		return (-a_gamma(a_h,a_s,gammat(n_f),n_f));
	}
public:
	inline U_exp_t(uint order_,A_strong &a_strong_):U_exp(order_,a_strong_){}
};

/*class U_exp_H:public U_exp{
//U_exp_S + U_exp_t
protected:
	inline Pol< complex<double> > u_exp(complex<double> a_h, complex<double> a_s, complex<double> L_h, uint n_f){
		return ( S(a_h,a_s,n_f)*2.*C_A						//C_S
				 -a_gamma(a_h,a_s,gammaS(n_f),n_f)			//C_S
				 -a_gamma(a_h,a_s,gammat(n_f),n_f)			//C_t
				 -a_gamma(a_h,a_s,Gamma(n_f),n_f)*L_h*C_A );//C_S
	}
public:
	inline U_exp_H(uint order_, Alpha_s const &alpha_s_){ init(order_,alpha_s_);}
};*/


class C_exp:protected Pol<Pol<complex<double> > >{
friend class C_H_exp;
protected:
	uint n_f_h;
	A_strong a_strong;

	inline virtual void init(){
		cout << "*************************************************************************" << endl
			 << "**                      class C_exp::init() not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}

	inline C_exp(uint order_,A_strong const &a_strong_):n_f_h(0),a_strong(a_strong_){
		order(order_);
	};

public:
	/*inline Pol<Pol<complex<double> > >const operator() (complex<double> a_h, complex<double> L_h, uint n_f_){
		//cout << "\tC_exp::operator()" << endl;
		if(n_f_h!=n_f_) {
			n_f_h=(n_f_);
			init();
		}
		Pol<Pol<complex<double> > > tmp(*this);
		//cout << "\ttmp\t\t=" << tmp << endl;
		tmp.setX(a_h);
		for(uint i=0; i<=_order; ++i) tmp[i].setX(L_h);
		//cout << "\ttmp("<< a <<"," <<L_h <<")\t=" << tmp << endl;
		return tmp;
	}*/
	inline Pol<Pol<complex<double> > >const operator() (complex<double> M, complex<double> mu_h){
		//cout << "\tC_exp::operator()" << endl;
		uint n = get_n_f(mu_h);
		if(n_f_h!=n) {
			n_f_h=(n);
			init();
		}
		Pol<Pol<complex<double> > > tmp(*this);

		complex<double> a_h = a_strong(mu_h);
		complex<double> L_h = 2.*log(M/mu_h);

		tmp.setX(a_h);
		for(uint i=0; i<=_order; ++i) tmp[i].setX(L_h);
		//cout << "\ttmp("<< a <<"," <<L_h <<")\t=" << tmp << endl;
		return tmp;
	}

	inline virtual ~C_exp(){}
};

class C_V_exp: public C_exp{
protected:

	inline Pol<complex<double> > _HVF(){
		Pol<complex<double> > _HVF(4,0.);
		_HVF[0]=0.125*255+0.5*7*PI2-83.0*PI4/360.0-30*zeta_3;
		_HVF[1]=-0.5*45-0.5*3*PI2+24*zeta_3;
		_HVF[2]=0.5*25-PI2/6.0;
		_HVF[3]=-3;
		_HVF[4]=0.5;
		return _HVF;
	}

	inline Pol<complex<double> > _HVA(){
		Pol<complex<double> > _HVA(4,0.);
		_HVA[0]=- 51157.0/648.0 - 337.0*PI2/108.0+11.0*PI4/45.0+313.0*zeta_3/9.0;
		_HVA[1]=2545.0/54.0+11.0*PI2/9.0-26*zeta_3;
		_HVA[2]=-233.0/18.0+PI2/3.0;
		_HVA[3]=11.0/9.0;
		_HVA[4]=0.0;
		return _HVA;
	}


	inline Pol<complex<double> > _HVf(){
		Pol<complex<double> > _HVf(4,0.);
		_HVf[0]=4085.0/162.0+23.0/27.0*PI2+4.0/9.0*zeta_3;
		_HVf[1]=-418.0/27.0-4.0*PI2/9.0;
		_HVf[2]=38.0/9.0;
		_HVf[3]=-4.0/9.0;
		_HVf[4]=0.0;
		return _HVf;
	}

	inline Pol<complex<double> > _CV1(){
		Pol<complex<double> > _CV1(2,0.);
		_CV1[0]=C_F*(-8+PI2/6.0);
		_CV1[1]=C_F*3;
		_CV1[2]=-C_F;
		return _CV1;
	}


	inline Pol<complex<double> > _CV2(int n_f){
		return ((_HVF()*C_F + _HVA()*C_A + _HVf()*T_F*n_f)* C_F );
	}


	inline void init(){
		uint n=_order;
		order(2);
		_P[0]=Pol<complex<double> >(0,0.);
		_P[1]=_CV1();
		_P[2]=_CV2(n_f_h);
		//cout << "\tC_V_exp::update() -> " << *this << endl;
		order(n);
	}

public:
	inline C_V_exp(uint order_,A_strong const &a_strong_):C_exp(order_,a_strong_){}
};

class C_S_exp: public C_exp{
friend class C_H_exp;
protected:

	inline Pol<complex<double> > _HSA(){
		Pol<complex<double> > _HSA(4,0.);
		_HSA[0] = 5105./162. + 67.*PI2/36. + PI4/72. - 143.*zeta_3/9.;
		_HSA[1] = 80./27. - 11.*PI2/9. - 2.*zeta_3;
		_HSA[2] = -67./9. + PI2/6.;
		_HSA[3] = 11./9.;
		_HSA[4] = 0.5;
		return _HSA;
	}

	inline Pol<complex<double> > _HSF(){
		Pol<complex<double> > _HSF(1,0.);
		_HSF[0] = -67./3.+16.*zeta_3;
		_HSF[1] = 4.;
		return _HSF;
	}

	inline Pol<complex<double> > _HSAF(){
		Pol<complex<double> > _HSAF(3,0.);
		_HSAF[0] = (-1832./9.-5.*PI2-92*zeta_3)/9.;
		_HSAF[1] = (104./3.+4.*PI2)/9.;
		_HSAF[2] = 20./9.;
		_HSAF[3] = -4./9.;
		return _HSAF;
	}

	inline Pol<complex<double> > _CS0(){
		Pol<complex<double> > _CS0(0,1.);
		return _CS0;
	}

	inline Pol<complex<double> > _CS1(){
		Pol<complex<double> > _CS1(2,0.);
		_CS1[0]=C_A*PI2/6.0;
		_CS1[1]=0.0;
		_CS1[2]=-C_A;
		return _CS1;
	}

	inline Pol<complex<double> > _CS2(int n_f){
		return (_HSA()*gsl_pow_2(C_A) + (_HSAF()*C_A + _HSF()*C_F)*T_F*n_f);
	}

	inline void init(){
		uint n=_order;
		order(2);
		_P[0]=Pol<complex<double> >(0,0.);
		_P[1]=_CS1();
		_P[2]=_CS2(n_f_h)-_P[1]*_P[1]*0.5;
		order(n);
	}

public:
	inline C_S_exp(uint order_,A_strong const &a_strong_):C_exp(order_,a_strong_){}
};

class C_t_exp: public C_exp{
friend class C_H_exp;
protected:

	inline Pol< complex<double> > _HtA(int n_f){
		Pol< complex<double> > _HtA(1,0.);
		_HtA[0] = 1063./36*C_A - 5./6.*T_F - 47./9.*T_F*n_f;
		_HtA[1] = -7.*C_A;
		return _HtA;
	}

	inline Pol< complex<double> > _HtF(int n_f){
		Pol< complex<double> > _HtF(1,0.);
		_HtF[0] = 0.5*27.*C_F - 4./3*T_F - 5.*T_F*n_f;
		_HtF[1] =-8.*T_F*n_f;
		return _HtF;
	}

	inline Pol< complex<double> > _HtAF(){
		Pol< complex<double> > _HtAF(1,0.);
		_HtAF[0] = -100./3.;
		_HtAF[1] = 11.;
		return _HtAF;
	}

	inline Pol< complex<double> > _Ct0(){
		Pol< complex<double> > _Ct0(0,1.);
		return _Ct0;
	}

	inline Pol< complex<double> > _Ct1(){
		Pol< complex<double> > _Ct1(0,0.);
		_Ct1[0]=5.*C_A-3.*C_F;
		return _Ct1;
	}

	inline Pol< complex<double> > _Ct2(int n_f){
		return (_HtA(n_f)*C_A + _HtAF()*C_A*C_F + _HtF(n_f)*C_F );
	}


	inline void init(){
		uint n=_order;
		order(2);
		_P[0]=Pol<complex<double> >(0,0.);
		_P[1]=_Ct1();
		_P[2]=_Ct2(n_f_h)-_P[1]*_P[1]*0.5;
		order(n);
	}

public:
	inline C_t_exp(uint order_,A_strong const &a_strong_):C_exp(order_,a_strong_){}
};

/*class C_H_exp: public C_exp{
protected:

	C_S_exp CS;
	C_t_exp Ct;

	inline void update(){
		CS.C_exp::update(n_f);
		Ct.C_exp::update(n_f);
		for(uint i=0; i<_order; ++i) _P[i]=CS[i]+Ct[i];
	}

public:
	inline C_H_exp(uint order_, uint n_f_)
	:CS(_order,n_f)
	,Ct(_order,n_f){
		C_exp::init(order_,n_f_);
	}
};*/


class _Hardfunc: public Pol< double >{
protected:
	CS_params* p;

	U_exp *U;
	C_exp *C;

	complex<double> M;
	complex<double> mu_h;
	complex<double> mu;


	//bool exponent;
	//A_strong a_strong;
	//uint n_f;

	inline _Hardfunc():U(NULL),C(NULL){}

	inline virtual void v_init(CS_params &p_){
		cout << "*************************************************************************" << endl
			 << "**                  class Hardfunc::v_init not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}
	inline virtual void v_update(){
		cout << "*************************************************************************" << endl
			 << "**                class Hardfunc::v_update not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}

	inline _Hardfunc(CS_params &p_):U(NULL),C(NULL){
		//init(p_);
	}

public:
	inline void init(CS_params &p_){
		p=(&p_);
		//cout << "->_Hardfunc::init("<<(void*)&p_<<") -> " << (void*)p << endl;
		M = complex<double>(0.,p->M);

		order(p->_order);
		v_init(p_);
		//cout << "_order="<<_order << endl;
	}

	inline void update(){
		//cout << "->_Hardfunc::update()" << endl << *p << endl;
		order(p->_order);
		if(p->piquadrat) mu_h = complex<double>(0.,p->mu_h);
		else mu_h = p->mu_h;
		mu = p->mu;
		//cout << "\torder = " << _order << "\tmu_h = " << mu_h << "\tmu = " << mu << endl;
		v_update();
		if(p->exponent) {
			_P[0]=von(1.);
			order(0);
		}
		_exp();
		if(p->improved) Pol_a_Improve(*this);
		order(p->order);
		//cout << "H(mu="<<p->mu<<",\tmu_h="<<p->mu_h<<"\tmu_t="<<p->mu_t << ")  a_h=" << p->a_h << "\tn_h=" << p->n_f_h << "\ta=" << p->a << "\tn=" << p->n_f
	//		<</* "\tH=" << *this <<*/ " = " << von(1) << endl;
		//cout << "H_new=" <<von(1.)<< endl;
		//cout << "<-Hardfunc::update()" << endl;
		//cout << "<-_Hardfunc::update()" << endl;
	}

	inline void update(double qT_, double mu_){
		//cout << "->_Hardfunc::update("<<qT_<<", "<<mu_<<")" << (void*)p << endl;
		p->update(qT_,mu_);
		//cout << p << endl;
		update();
	}

	inline virtual ~_Hardfunc(){
		delete U;
		delete C;
	}

	/*inline Pol<double> & operator ()(){
		update();
		return *this;
	}*/



};

std::ostream& operator<<( std::ostream& os, const _Hardfunc* H );

class Hardfunc_V: public _Hardfunc{
protected:
	inline void v_init(CS_params &p_){
		delete U;
		delete C;
		U=new U_exp_V(_order,p->a_s);
		C=new C_V_exp(_order,p->a_s);
	}

	inline void v_update(){
		//Pol< Pol<complex<double> > > tmp((*C)(p->a_h,p->L_h,p->n_f));
		//U->update(p->mu_h,p->a_h,p->n_f_h,p->mu,p->a,p->n_f,p->L_h);
		Pol< Pol<complex<double> > > tmp((*C)(M,mu_h));
		U->update(M,mu_h,mu);
		for(uint i=0; i<=_order; ++i) _P[i]=(tmp[i].von(1.)+(*U)[i]).real()*2.;
	}

public:
	inline Hardfunc_V(CS_params &p_):_Hardfunc(p_){
		init(p_);
	}

	inline ~Hardfunc_V(){
	}
};

class Hardfunc_H: public _Hardfunc{
protected:
	U_exp *Ut;
	C_exp *Ct;


	inline void v_init(CS_params &p_){
		delete U;
		delete C;
		delete Ut;
		delete Ct;
		U=new U_exp_S(_order,p->a_s);
		C=new C_S_exp(_order,p->a_s);
		Ut=new U_exp_t(_order,p->a_s);
		Ct=new C_t_exp(_order,p->a_s);
	}

	inline void v_update(){
		//cout << "->Hardfunc_H::v_update()" << endl;
		Pol< Pol<complex<double> > > tmp((*C)(M,mu_h)+(*Ct)(M,p->mu_t));
		U->update(M,mu_h,mu);
		Ut->update(M,p->mu_t,mu);

		for(uint i=0; i<=_order; ++i) _P[i]=(tmp[i].von(1.)+(*U)[i]+(*Ut)[i]).real()*2.;
		//cout << "<-Hardfunc_H::v_update()" << endl;
	}

public:
	inline Hardfunc_H(CS_params &p_):_Hardfunc(p_),Ut(NULL),Ct(NULL){
		init(p_);
	}

	inline ~Hardfunc_H(){
		delete Ut;
		delete Ct;
	}
};

class Hardfunc{
protected:
	_Hardfunc* H;
public:
	CS_params p;
	inline void init(CS_params const &p_){
		p.init(p_);
		//cout << "->Hardfunc::init("<<(void*)&p_<<") ->" << (void*)&p << endl;
		delete H;
		if(p.prod==PROD_H) H = new Hardfunc_H(p);
		else H = new Hardfunc_V(p);
		//cout << "<-Hardfunc::init("<<(void*)&p_<<") ->" << (void*)&p << endl;
	}

	inline Pol< double >& operator ()(double qT, double mu){
		//cout << "->Hardfunc("<<qT<<", "<<mu<<")" << (void*)&p << endl;
		H->update(qT,mu);
		//cout << H << endl;
		return *H;
	}

	inline Hardfunc(CS_params const &p_):H(NULL){	init(p_);}
};

#endif
