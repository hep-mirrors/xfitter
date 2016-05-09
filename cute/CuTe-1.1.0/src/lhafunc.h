/**
* @file lhafunc.h
* @author Daniel Wilhelm
* @date	26.12.2012
* @brief Functions to handle the LHAPDF library
*/

#ifndef LHAFUNC_h
#define LHAFUNC_h

#include "konst.h"
#include "vecmath.h"
#include "polynom.h"
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_roots.h>
//#define ONLY_QUARKS
//#define ONLY_GLUONS
//#define ONLY_MIXED

#include "xfitter_cpp.h"

using std::cout;
using std::endl;
using std::complex;




inline double get_threshold(lha_parton q){
	//herafitter
	return flav_threshold_(q);
	/*
	double getThreshold = LHAPDF::getThreshold(q);
	if(getThreshold==0) return LHAPDF::getQMass(q);
	return getThreshold;
	*/
}

inline uint get_n_f(double mu){
	mu=fabs(mu);
	//herafitter
	//	return flav_number_(mu);
	double threshold;
	for(lha_parton q=lha_t; q>lha_g; --q){
		threshold = get_threshold(q);
		//cout << "threshold_"<<q << "("<<mu<<") = " << threshold << endl;
		if(threshold==-1.) continue;
		if(mu>threshold) return (q);
	}
	return 0;
}

inline uint get_n_f(complex<double> mu){
	return get_n_f(mu.real()+mu.imag());
}

class A_strong{
public:
	double Lambda_QCD;
protected:
	static const uint oMax=3;
	double mu_0;
	double a_0;
	double gen;
	uint n_f;
	uint order;

	inline Pol< complex<double> > apol(complex<double> L)const{
		complex<double> LL = log(L);
		Pol< complex<double> > apol(oMax,0.);
		double b10 = beta_1_0(n_f);
		double b20 = beta_2_0(n_f);
		double b30 = beta_3_0(n_f);
		apol[0]=1.;
		apol[1]=-b10*LL;
		apol[2]=(b20+(-1.+LL*(-1.+LL))*gsl_pow_2(b10));
		apol[3]=0.5*(b30-6.*b10*b20*LL+gsl_pow_3(b10)*(-1.+LL*(4. +LL*(5.-2.*LL))) );

		L=1./(L*beta_0(n_f));
		LL=L;
		for(uint i=0; i<=oMax; ++i){
			apol[i]*=LL;
			LL*=L;
		}
		apol.order(order);
		return apol;
	}

	inline Pol< double > Dapol(double mu, double Lambda)const{
		double L=log(gsl_pow_2(mu/Lambda));
		double LL=log(L);
		Pol< double > Dapol(oMax,0.);
		double b10 = beta_1_0(n_f);
		double b20 = beta_2_0(n_f);
		double b30 = beta_3_0(n_f);
		Dapol[0]=1.;
		Dapol[1]=b10*(1.-2.*LL);
		Dapol[2]=3.*b20+gsl_pow_2(b10)*(-2.+LL*(-5.+3.*LL));
		Dapol[3]=3.*b10*b20+2.*b30-12.*b10*b20*LL+gsl_pow_3(b10)*(-4.+LL*(3.+LL*(13.-4.*LL)));

		LL=2./(L*Lambda);
		L=1./(L*beta_0(n_f));
		for(uint i=0; i<=oMax; ++i){
			LL*=L;
			Dapol[i]*=LL;
		}
		Dapol.order(order);
		return Dapol;
	}

	inline Pol< double > bpol(double mu, double Lambda)const{
		double L=log(gsl_pow_2(mu/Lambda));
		double LL=log(L);
		Pol< double > bpol(oMax,0.);
		double b10 = beta_1_0(n_f);
		double b20 = beta_2_0(n_f);
		double b30 = beta_3_0(n_f);
		bpol[0]=1.;
		bpol[1]=b10*(1.-2.*LL);
		bpol[2]=3.*b20+gsl_pow_2(b10)*(-2.+LL*(-5.+3.*LL));
		bpol[3]=3.*b10*b20+2.*b30-12.*b10*b20*LL+gsl_pow_3(b10)*(-4.+LL*(3.+LL*(13.-4.*LL)));

		LL=-2./(L);
		L=1./(L*beta_0(n_f));
		for(uint i=0; i<=oMax; ++i){
			LL*=L;
			bpol[i]*=LL;
		}
		bpol.order(order);
		return bpol;
	}

	inline static double f(double Lambda, void *params){
		A_strong *a = (A_strong*) params;
		complex<double>  L(log(gsl_pow_2(a->mu_0/Lambda)),0.);
		return (a->apol(L).von(1.).real() - a->a_0);
	}

	inline static double df(double Lambda, void *params){
		A_strong *a = (A_strong*) params;
		return a->Dapol(a->mu_0,Lambda).von(1.);
	}

	inline static void fdf (double Lambda, void *params, double *y, double *dy){
	     *y = f(Lambda,params);
	     *dy = df(Lambda,params);
	}

	inline int init_lambda(double estimate){
		int status;
		int iter = 0, max_iter = 1000;
		const gsl_root_fdfsolver_type *T;
		gsl_root_fdfsolver *s;
		double x0;
		double x = estimate;
		gsl_function_fdf FDF;

		FDF.f = &f;
		FDF.df = &df;
		FDF.fdf = &fdf;
		FDF.params = this;

		T = gsl_root_fdfsolver_steffenson;
		//T = gsl_root_fdfsolver_newton;
		s = gsl_root_fdfsolver_alloc (T);
		gsl_root_fdfsolver_set (s, &FDF, x);

		do{
			iter++;
			status = gsl_root_fdfsolver_iterate (s);
			GSL_ERROR_CODE(status);
			x0 = x;
			x = gsl_root_fdfsolver_root (s);
			status = gsl_root_test_delta (x, x0, 0, gen);
			GSL_ERROR_CODE(status);
		}while (status == GSL_CONTINUE && iter < max_iter);
		Lambda_QCD=x;
		gsl_root_fdfsolver_free (s);
		return status;
	}

	inline void init(double mu_0_){
		mu_0=(mu_0_);
		//		a_0=0.25*LHAPDF::alphasPDF(mu_0)*M_1_PI;
		a_0=0.25*appl_fnalphas_(mu_0)*M_1_PI;
		//order=LHAPDF::getOrderAlphaS();
		order=steering_.i_fit_order_ - 1;
		n_f=get_n_f(mu_0);
		init_lambda(0.2);
	}
public:


	inline A_strong(double gen_=1e-8):gen(gen_){}

	inline double operator() (double mu){
	  //		return 0.25*LHAPDF::alphasPDF(mu)*M_1_PI;
		return 0.25*appl_fnalphas_(mu)*M_1_PI;
	}

	inline complex<double> operator() (complex<double> mu){
		if(mu.imag()==0.) return operator ()(mu.real());
		else{
			init(mu.real()+mu.imag());
			return apol(2.*log(mu/Lambda_QCD)).von(1.);
		}
	}

	inline double b(double mu){
		init(mu);
		/*cout << "A_strong::init(" << mu << ")_"<< order<<" = " << Lambda_QCD
			 <<	"\ta = " << apol(2.*log(mu/Lambda_QCD)).von(1.) << " = " << 0.25*LHAPDF::alphasPDF(mu)*M_1_PI
			 << "\tb = " << bpol(mu,Lambda_QCD).von(1.) << endl;*/
		return bpol(mu,Lambda_QCD).von(1.);
	}
};


extern Matrix< 2, complex<double> > Pauli(int i);


typedef Matrix<N_FMAX+1,double> CKM;

template<class T>
inline T sgn(T t)
{
  if(t < T(0)) return T(-1);
  else if (t > T(0)) return T(1);
  return T(0);
}


class Charge{
protected:
	//static const int N=2*N_FMAX+1;
	//double* qf;
	//complex<double>** r;
	//complex<double>** l;
	vector< double > qf;
	vector< vector< complex<double> > > r;
	vector< vector< complex<double> > > l;
	Matrix<2,complex<double> > pauli_three;
	Matrix<2,complex<double> > pauli_plus;
	//Matrix<2,complex<double> > pauli_minus;
	CKM ckm;
	//Matrix<N_FMAX+1,double> ckm;

	inline void init_qf(){
		double tmp = 1./3.;
		qf[e_tbar] = -tmp*2.;
		qf[e_bbar] = tmp;
		qf[e_cbar] = -tmp*2.;
		qf[e_sbar] = tmp;
		qf[e_ubar] = -tmp*2.;
		qf[e_dbar] = tmp;
		qf[e_g] = 0.;
		qf[e_d] = -tmp;
		qf[e_u] = tmp*2.;
		qf[e_s] = -tmp;
		qf[e_c] = tmp*2.;
		qf[e_b] = -tmp;
		qf[e_t] = tmp*2.;
	}

	inline complex<double> l_Photon(lha_parton i, lha_parton j){
		if(i==-j) return qf[e_(j)];
		else return 0.;
	}
	inline complex<double> r_Photon(lha_parton i, lha_parton j){
		return l_Photon(i,j);
	}
	inline complex<double> r_W(lha_parton i, lha_parton j){
		return 0.;
	}
	inline complex<double> l_W(lha_parton i, lha_parton j){
		return pauli_plus[lha_uptype(i)][lha_uptype(j)]*M_SQRT1_2/SThetaW*ckm[abs(j)][abs(i)];
	}
	inline complex<double> l_Z(lha_parton i, lha_parton j){
		if(i==-j) return pauli_three[lha_downtype(i)][lha_downtype(j)]/sin(2.*ThetaW)*sgn((double)j)-qf[e_(j)]*tan(ThetaW);
		return 0.;
	}
	inline complex<double> r_Z(lha_parton i, lha_parton j){
		if(i==-j) return -qf[e_(j)]*tan(ThetaW);
		return 0.;
	}

	inline complex<double> r_von(PROD prod, lha_parton i, lha_parton j){
		if( i*j >= 0) return 0.;
		if( j<0) return -conj(l_von(prod,j,i));
		if(prod==PROD_PHOTON)	return r_Photon(i, j);
		if(prod==PROD_W)		return r_W(i, j);
		if(prod==PROD_Z) 		return r_Z(i, j);
		return 0.;
	}
	inline complex<double> l_von(PROD prod, lha_parton i, lha_parton j){
		if( i*j >= 0) return 0.;
		if( j<0) return -conj(r_von(prod,j,i));
		if(prod==PROD_PHOTON)	return l_Photon(i, j);
		if(prod==PROD_W)		return l_W(i, j);
		if(prod==PROD_Z) 		return l_Z(i, j);
		return 0.;
	}



public:
	inline void init(PROD prod){
		//cout << "->Charge::init(" << prod << ")" << endl;
		init_qf();
		for(int i=0; i<N_PARTONS; ++i) for(int j=0; j<N_PARTONS; ++j){
			r[i][j]=r_von(prod,lha_(i),lha_(j));
			l[i][j]=l_von(prod,lha_(i),lha_(j));

		}
		//cout << endl << "<-Charge::init(" << prod << ")" << endl;
	}
	inline void init(PROD prod, CKM ckm_){
		ckm = ckm_;
		init(prod);
	}

	inline Charge(PROD prod=PROD_WRONG_STATEMENT, CKM ckm_=CKM(0.) )
	//inline Charge(PROD prod, Matrix< N_FMAX+1,double > ckm_)
	:pauli_three(Pauli(3))
	,pauli_plus((Pauli(1)+Pauli(2)*complex<double>(0.,1.))*0.5)
	,ckm(ckm_){
		qf.resize(N_PARTONS);
		r.resize(N_PARTONS);
		l.resize(N_PARTONS);
		for(int i=0; i<N_PARTONS; i++){
			r[i].resize(N_PARTONS);
			l[i].resize(N_PARTONS);
		}
		init(prod);
	}
	//inline Charge(PROD prod):Charge(prod,Matrix<N_FMAX+1,double>(0.)){}
	//inline Charge():Charge(PROD_WRONG_STATEMENT){}

	inline ~Charge(){}

	inline double q(e_parton i){ return qf[i];}
	inline complex<double> R(e_parton i, e_parton j){ return r[i][j];}
	inline complex<double> L(e_parton i, e_parton j){ return l[i][j];}

	inline double RR(e_parton i, e_parton j){ return norm(r[i][j]);}
	inline double LL(e_parton i, e_parton j){ return norm(l[i][j]);}
	inline double Q(e_parton i, e_parton j){ return RR(i,j)+LL(i,j);}

};


class Q_Matrix:public Matrix<N_PARTONS,double>{
protected:
	//PROD const _prod;
	PROD _prod;
	Charge C;
	CKM ckm;


	inline void init(){
		for(e_parton i=0; i<N_PARTONS; ++i)for(e_parton j=0; j<N_PARTONS; ++j) {
			_M[i][j]=0.5*C.Q(i,j);
		}
		if(_prod==PROD_H) _M[e_g][e_g]=1.;
	}

public:

	inline Q_Matrix(PROD prod=PROD_WRONG_STATEMENT, CKM ckm_=CKM(0.))// Matrix<N_FMAX+1,double> ckm_=Matrix<N_FMAX+1,double>(0.))
	:Matrix<N_PARTONS,double>(0.)
	,_prod(prod)
	,C(prod,ckm_)
	,ckm(ckm_){
		init();
	}
	inline Q_Matrix(Q_Matrix const &Q):Matrix<N_PARTONS,double>(Q),_prod(Q.prod()){}

	inline PROD const prod()const{return _prod;}

	inline Vector<N_PARTONS,double> &operator [] (e_parton const i){ return _M[i];}
	inline Vector<N_PARTONS,double> const &operator [] (e_parton const i)const {return _M[i];}

	inline double const operator ()(lha_parton const i, lha_parton const j)const{
		return _M[e_(i)][e_(j)];
	}

};


class VectorPhi: public Vector< N_PARTONS,double >{
protected:
	int _sign;
public:
	inline void update(double x, double q){
	  //	        LHAPDF::xfx(x,q,_V);
	        appl_fnpdf_(x,q,_V);

		double c=1./x;
		for(e_parton i=0; i<N_PARTONS; ++i){
			_V[i]*=c;
		}
		if(_sign==-1)for(e_parton i=0; i<e_g; ++i){
			double tmp=_V[i];
			_V[i]=_V[bar(i)];
			_V[bar(i)]=tmp;
		}

#ifdef ONLY_QUARKS
		_V[e_g]=0.0;
#endif
#ifdef ONLY_GLUONS
		for(e_parton i=0; i<e_g; ++i){
			_V[i]=_V[bar(i)]=0.0;
		}
#endif
	}

	inline VectorPhi(int const sign=+1)
		:Vector< N_PARTONS,double >(0.)
		,_sign(sign){
	}

	inline VectorPhi(double x, double q, int const sign=1)
		:Vector< N_PARTONS,double >(0.)
		,_sign(sign){
		update(x,q);
	}

	inline double const &operator [] (e_parton const i)const {return _V[i];}
	//inline double const &operator [] (lha_parton const i)const {return _V[e_(i)];}

	inline double const operator ()(lha_parton const i)const{return _V[e_(i)];}
	//inline double const operator ()(e_parton const i)const{return _V[i];}
	inline double const operator ()(lha_parton const i, double x, double q)const
	{
	  //	  return (LHAPDF::xfx(x,q,_sign*i)/x);
	  double f[N_PARTONS];
	  appl_fnpdf_(x,q,f);
	  return f[_sign*i+6]/x;
	}
	//inline double const operator ()(e_parton const i, double x, double q)const{return (LHAPDF::xfx(x,q,_sign*lha_(i))/x);}

	inline double const sumofquarks(int n_f)const{
		double sum=0.;
		for(e_parton i=e_(n_f); i>e_g; --i) sum+=(_V[i]+_V[bar(i)]);
		return sum;
	}

	inline int const sign()const{ return _sign;};

	inline double const sumofQphi(Q_Matrix const &Q,int n_f)const{
		Vector<N_PARTONS,double> tmp(*this*Q);
		double sum=0.;
		for(e_parton i=e_(n_f); i>e_g; --i) sum+=(tmp[i]+tmp[bar(i)]);
		return sum;
	}
};

extern std::ostream& operator<<( std::ostream& os, VectorPhi const &V );



#endif
