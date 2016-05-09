#ifndef PARAM_H
#define PARAM_H

#include "lhafunc.h"

using std::complex;
using std::cout;
using std::endl;


/*static const std::string PDF_NAME = "MSTW2008nnlo68cl";
static const bool PDF_CENTER = true;*/
//|A_q|^2=35.1201

/*static const std::string PDF_NAME = "NNPDF21_as_0117_100";
static const bool PDF_CENTER = false;*/
//|A_q|^2=1.08458

/*static const std::string PDF_NAME = "JR09VFnnloE";
static const bool PDF_CENTER = true;*/
//|A_q|^2=1.07991

/*static const std::string PDF_NAME = "HERAPDF15NNLO_EIG";
static const bool PDF_CENTER = true;*/
//|A_q|^2=1.07889


//static vector<vector<double> > debug;


class Q_star{
protected:
	double q_star;
	A_strong a_strong;
	double G_0_B;

	double gen;
	double M;

	inline double eta(double mu){
		return G_0_B*a_strong(mu)*2.*log(M/mu);
	}

	inline double deta(double mu){
		return 2.*G_0_B*(log(M/mu)*a_strong.b(mu)-a_strong(mu))/mu;
	}

	inline static double f(double x, void* params){
		Q_star* p = (Q_star*) params;
		//debug[0].push_back(x);
		//debug[1].push_back(p->eta(x));
		//cout << "\tx = " << x << endl;
		//cout << "\ta = " << p->a_strong(x) << "\tb = " << p->a_strong.b(x) << endl;
		//cout << "\tf = " << p->eta(x)-1.;
		return (p->eta(x)-1.);
	}

	inline static double df(double x, void* params){
		Q_star* p = (Q_star*) params;
		//debug[2].push_back(p->deta(x));
		//cout << "\tdf = " << p->deta(x) << endl;
		return p->deta(x)*10;// Times 10, so walk 10 times slower, thus the danger to walk into alpha_s=inf is smaller.
	}

	inline static void fdf (double x, void *params, double *y, double *dy){
		*y = f(x,params);
	    *dy = df(x,params);
	}

	inline int NSolve(double estimate){
		//cout << "Q_star::NSolve("<<estimate<<","<<gen<<")" << endl;
	        //debug.resize(3);
		int status;
		int iter = 0, max_iter = 10000;
		const gsl_root_fdfsolver_type *T;
		gsl_root_fdfsolver *s;
		double x0;
		double x = estimate;
		gsl_function_fdf FDF;

		FDF.f = &f;
		FDF.df = &df;
		FDF.fdf = &fdf;
		FDF.params = (void*)this;

		T = gsl_root_fdfsolver_steffenson;
		//T = gsl_root_fdfsolver_newton;
		s = gsl_root_fdfsolver_alloc (T);
		gsl_root_fdfsolver_set (s, &FDF, x);
		do{
			iter++;
			//cout << iter << endl;
			status = gsl_root_fdfsolver_iterate (s);
			GSL_ERROR_CODE(status);
			x0 = x;
			x = gsl_root_fdfsolver_root (s);
			//x=fabs(x);
			status = gsl_root_test_delta (x, x0, 0, gen);
			GSL_ERROR_CODE(status);
		}while (status == GSL_CONTINUE && iter < max_iter);
		q_star=x;
		gsl_root_fdfsolver_free (s);
		return status;
	}

public:

	inline Q_star(){}

	inline double operator() (double M_, double G_0_B_, double gen_=1e-8){
		G_0_B = (G_0_B_);
		M = (M_);
		gen = (gen_);
		NSolve(2.*G_0_B/3.);
		return q_star;
	}
};





class CS_params{
protected:
	Q_star q_star_c;
public:
	A_strong a_s;
	double s;
	double M;
	double qT;
	double mu;
	double mu_h;
	double mu_t;
	double y;
	int _order;
	uint order;

	bool improved;
	bool exponent;
	bool piquadrat;
	COL col;
	PROD prod;
	RES res;
	CUTOFF cutoff;
	double Lambda;

	CKM ckm;
	//Matrix< N_FMAX+1, double > ckm;

//	init
	//int sign;
	int sign[2];
	double C_B;	// C_A || C_F

// update_y
	double xi[2];	//sTau*exp(+-y)

// update qT, mu
	uint n_f;
	double a;//alpha_s(mu)/(4Pi)
	double sTau; //sqrt((M^2+qT^2)/s)
	double eta; //Gamma_0*a*log(M^2/mu^2)
	double h_B;	// beta_0/C_A  ||  -gamma_0q/C_F
	double L_mu; //log(mu^2/qT^2)
	double y_max;

//update mu_h
	uint n_f_h;
	complex<double> a_h;	//alpha_s(mu_h)/(4Pi)

//update mu_t
	uint n_f_t;
	double a_t;	//alpha_s(mu_t)/(4Pi)

	//bool PDFcenter;
	PDF_ERR PDFerr;
	int PDFn;	//PDFcount
	double gen;


	inline void init( uint order_, bool improved_, bool exponent_, bool piquadrat_, COL col_,
						PROD prod_, RES res_, CKM const &ckm_, double gen_, PDF_ERR PDFerr_/*, bool PDFcenter_*/, int PDFn_, CUTOFF cutoff_=CUTOFF_NO, double Lambda_=0.){
						_order=(order_); improved=(improved_); exponent=(exponent_); piquadrat=(piquadrat_); col=(col_);
						prod=(prod_); res=(res_); ckm=(ckm_); gen=(gen_); PDFerr=(PDFerr_);/*PDFcenter=(PDFcenter_);*/ PDFn=(PDFn_); cutoff=(cutoff_); Lambda=(Lambda_);
		sign[0] = sign[1] = +1;
		n_f=5;
		if(col == COL_PbarP) sign[0] = -1;
		if(col == COL_PbarPbar) sign[0] = sign[1] = -1;
		if(prod==PROD_H){
			C_B=C_A;	//h_A in update(qT,mu)
			h_B=beta_0(n_f)/C_A;
		}
		else{
			C_B=C_F;
			h_B=-gamma_0q/C_F;
		}
		if(res_==RES_LLO || res_==RES_LO){
			improved=false;
			piquadrat=false;
			exponent=false;
			cutoff_=CUTOFF_NO;
		}
		order=(_order*(1+improved_));
	}

	inline void init(double s_, double M_, double mu_t_, double mu_h_){
		//cout << "CS_paramsinit("<<s_<<","<<M_<<","<<mu_t_<<","<<mu_h_<<")\t";
		s=(s_); M=(M_); mu_t=(mu_t_); mu_h=(mu_h_);

		n_f_t=get_n_f(mu_t);
		a_t=a_s(mu_t);
		n_f_h=get_n_f(mu_h);
		//cout << "n_f_h=" << n_f_h;
		if(piquadrat){
			a_h=a_s(complex<double>(0.,mu_h));
		}
		else{
			a_h=a_s(mu_h);
		}
		//cout << endl;
	}

	inline void update(double y_){
		y=(y_);
		xi[0] = sTau*exp(y);
		xi[1] = sTau*exp(-y);
	}

	inline void update(double qT_, double mu_){
		qT=(qT_);
		mu=(mu_);
		if(res==RES_LLO || res==RES_LO)	if(mu_h!=mu /*|| mu_t!=mu*/) init(s,M,mu_t,mu);

		n_f=get_n_f(mu);
		//cout << "n_f("<<mu<<") = "<<n_f<<endl;
		a=a_s(mu);

		sTau=sqrt((M*M+qT*qT)/s);
		eta=Gamma_0*a*2.*log(M/mu);
		if(prod==PROD_H) h_B=beta_0(n_f)/C_A;	//h_F in init
		L_mu=2.*log(mu/qT);
		y_max = -log( (M*M+s-sqrt(gsl_pow_2(M*M-s)-4*qT*qT*s )) / ( 2*sqrt((M*M+qT*qT)*s) ));
	}

	inline void update(double qT_, double mu_,double y_){
		update(qT_,mu_);
		update(y_);
	}

	inline void update(CS_params const &p_){
		update(p_.qT,p_.mu,p_.y);
	}

	inline CS_params(){
		init(0, false, false, false, COL_PP, PROD_WRONG_STATEMENT,RES_LL,CKM(0.),1e-4,PDF_ERR_WRONG_STATEMENT,/*true,*/0);
	}

	inline void init(CS_params const &p){
		init(p._order,p.improved,p.exponent,p.piquadrat,p.col,p.prod,p.res,p.ckm,p.gen,p.PDFerr/*,p.PDFcenter*/,p.PDFn, p.cutoff,p.Lambda);
		init(p.s,p.M,p.mu_t,p.mu_h);
	}

	inline CS_params(CS_params const &p){
		init(p);
	}

	inline void print(ostream& os)const{
		os << "--\t O=" << _order << "\timp=" << improved;
		if(improved) os << "\tO_eps=" << order;
		os << "\texp=" << exponent << "\tpi2=" << piquadrat << endl;

		os << "--\t";
		switch (col) {
		case COL_PP: os << "PP -> ";break;
		case COL_PbarP: os << "P\\bar(P) -> ";break;
		default:break;
		}
		switch (prod) {
		case PROD_H: os << "H\t";break;
		case PROD_PHOTON: os << "\\gamma\t";break;
		case PROD_W: os << "W\t";break;
		case PROD_Z: os << "Z\t";break;
		default:break;
		}
		for(uint i=1; i<=order; ++i) os << "N";
		switch (res) {
		case RES_LL: os << "LL\t";break;
		case RES_LLO: os << "LLO\t";break;
		case RES_LO: os << "LO\t";break;
		default:break;
		}
		os << endl;

		os 	<< "--\t s=" << s << "\tM=" << M << endl
			<< "--\tmu_t=" << mu_t << "\ta_t=" << a_t << "\tn_f_t=" << n_f_t << endl
			<< "--\tmu_h=" << mu_h << "\ta_h=" << a_h << "\tn_f_h=" << n_f_h << endl
			<< "--\t  mu=" << mu << "\t  a=" << a << "\t  n_f=" << n_f << "\tL_mu=" << L_mu << endl
			<< "--\tqT=" << qT << endl
			<< "--\ty=" << y << "\txi= (" << xi[0] << ", " << xi[1] << ")" << endl;
	}

	inline double q_star(){
		//cout << "q_star_c("<<M<<","<<Gamma_0*C_B<<","<<gen<<") = " << q_star_c(M,Gamma_0*C_B,gen)<< endl;
		return q_star_c(M,Gamma_0*C_B,gen);
	}
};


extern std::ostream& operator<<( std::ostream& os, CS_params const &p);

extern std::ostream& operator<<( std::ostream& os,CS_params const *p );


template < typename T >
inline Pol< T >& Pol_a_Improve(Pol< T > & P){
	uint n=P.order();
	P.order(2*n);
	for(uint i=n; i>0; --i){
		P[2*i]=P[i];
		P[i]-=P[i];
	}
	return P;
}

template < typename T >
inline Pol< Pol< T > >& Pol_a_L_Improve(Pol< Pol< T > > & P){
	uint n=P.order();
	//for(uint i=0; i<=P.order(); ++i) n=GSL_MAX_INT(n,2*i-P[i].order());
	T null = P[0][0]-P[0][0];
	Pol< Pol< T > > tmp;
	tmp.order(2*n);
	for(uint a=0; a<=n; ++a){
		uint L= P[a].order();
		for(uint l=0; l<=L; ++l){
			if(P[a][l]==null) continue;
			uint epsilon=GSL_MAX_INT(0,2*a-l);
			//cout << "tmp(a="<<a<<",l="<<l<<",e="<<epsilon<<")" << endl;
			if(tmp[epsilon].order()<l) tmp[epsilon].order(l);
			tmp[epsilon][l]=P[a][l];
			//cout << "\t=" <<tmp << endl;
		}
	}
	P=tmp;
	return P;
}
typedef double (*Integrand)(double, void*);



#endif
