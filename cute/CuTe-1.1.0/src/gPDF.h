#ifndef GPDF_H
#define GPDF_H

#include "param.h"
#include "apkernels.h"
#include "iostream"
#include <gsl/gsl_integration.h>
using std::ofstream;

using std::cout;
using std::endl;

//#define debBB





namespace gPDF_Kernels {
using namespace AP_Kernels;

//							------------------  R  ------------------

//------------ R_gg
inline double _R_ggf(double z, double fz, double f){
	return ( 0.0 );
}

inline double _R_ggf(double x, double f){
	return ( -C_A*PI2/6.*f);
}

inline double _R_ggf(double z, double x, double fz, double f){
	return ( _R_ggf(z, fz, f) +  _R_ggf(x, f)/(1.-x));
}

inline double R_ggf(double z, double x, double fz, double f, int n){
	return ( _R_ggf(z, fz, f) +  _R_ggf(x, f)/(1.-x));
}

//------------ R_gq

inline double R_gq(double z){
	return (2.*C_F*z);
}

inline double _R_gqf(double z, double x, double fz, double f){
	return (R_gq(z)*fz);
}

inline double R_gqf(double z, double x, double fz, double f, int n){
	return (R_gq(z)*fz);
}

//------------ R_qq
inline double _R_qqf(double z, double fz, double f){
	return ( C_F*(2.*(1.-z))*fz );
}

inline double _R_qqf(double x, double f){
	return ( -C_F*PI2/6.*f);
}

inline double _R_qqf(double z, double x, double fz, double f){
	return ( _R_qqf(z, fz, f) +  _R_qqf(x, f)/(1.-x));
}

inline double R_qqf(double z, double x, double fz, double f, int n){
	return ( _R_qqf(z, fz, f) +  _R_qqf(x, f)/(1.-x));
}

//------------ R_qg

inline double R_qg(double z){
	return (4.*T_F*z*(1.-z));
}

inline double _R_qgf(double z, double x, double fz, double f){
	return (R_qg(z)*fz);
}

inline double R_qgf(double z, double x, double fz, double f, int n){
	return (R_qg(z)*fz);
}
//							------------------  Q2  ------------------

//------------ Q_gg
inline double _Q_ggf(double z, double fz, double f,int n_f){
	return ( _D_ggf(z,fz,f,n_f) - 2.*beta_0(n_f)*_P_ggf(z,fz,f,n_f));
}

inline double _Q_ggf(double x, double f,int n_f){
	return ( _D_ggf(x,f,n_f) - 2.*beta_0(n_f)*_P_ggf(x,f,n_f));
}

inline double Q_ggf(double z, double x, double fz, double f,int n_f){
	return (D_ggf(z,x,fz,f,n_f) - 2.*beta_0(n_f)*P_ggf(z,x,fz,f,n_f) );
}

//------------ Q_qg

inline double _Q_gqf(double z, double fz, double f,int n_f){
	return ( _D_gqf(z,fz,f,n_f) - 2.*beta_0(n_f)*P_gq(z)*fz);
}

inline double _Q_gqf(double x, double f){
	return ( _D_gqf(x,f) );
}

inline double Q_gqf(double z, double x, double fz, double f,int n_f){
	return ( D_gqf(z,x,fz,f,n_f) - 2.*beta_0(n_f)*P_gqf(z,x,fz,f) );
}

//------------ Q_qq
inline double _Q_qqf(double z, double fz, double f,int n_f){
	return ( _D_qqf(z,fz,f) - 2.*beta_0(n_f)*_P_qqf(z,fz,f));
}

inline double _Q_qqf(double x, double f,int n_f){
	return ( _D_qqf(x,f) - 2.*beta_0(n_f)*_P_qqf(x,f ));
}

inline double Q_qqf(double z, double x, double fz, double f,int n_f){
	return ( D_qqf(z,x,fz,f) - 2.*beta_0(n_f)*P_qqf(z,x,fz,f) );
}

//------------ Q_qg

inline double _Q_qgf(double z, double fz, double f,int n_f){
	return ( _D_qgf(z,fz,f,n_f) - 2.*beta_0(n_f)*P_qg(z)*fz);
}

inline double _Q_qgf(double x, double f){
	return ( _D_qgf(x,f) );
}

inline double Q_qgf(double z, double x, double fz, double f,int n_f){
	return ( D_qgf(z,x,fz,f,n_f) - 2.*beta_0(n_f)*P_qgf(z,x,fz,f) );
}



}  // namespace gPDF_Kernels



using namespace gPDF_Kernels;


class VectorB: public Vector< N_PARTONS, Pol<Pol<double> > >{
friend class BQB_fast;

protected:
	//CS_params* p;
	PROD prod;
	uint order;
	int sign;
	VectorPhi phixi;
	VectorPhi phiz;
	double gen;
	//bool improved;


	double xi;
	double mu;
	double a;
	int n_f;

	DistFunc P_iq;
	DistFunc P_ig;
	e_parton q;


	inline static double _integrand_g(double z, void* params){
		double c=1./z;
		VectorB *B = (VectorB*) params;
		//if(z>1.-B->gen) return _integrand_g(1.-B->gen,params);
		double xi(B->xi);
		double x(c*xi);
		int n_f(B->n_f);
		B->phiz.update(x,B->mu);
		double soq_z(B->phiz.sumofquarks(n_f));
		double soq_xi(B->phixi.sumofquarks(n_f));
		double g_z(B->phiz[e_g]);
		double g_xi(B->phixi[e_g]);
	#ifdef ONLY_MIXED
		if(B->prod==PROD_H)
			return (B->P_iq(z, xi, c*soq_z,soq_xi,n_f));
		else
			return (B->P_ig(z, xi, c*g_z,g_xi,n_f));
	#endif
		return (B->P_iq(z, xi, c*soq_z,soq_xi,n_f)	+ B->P_ig(z, xi, c*g_z,g_xi,n_f));
	}

	inline static double _integrand_q(double z, void* params){
		double c=1./z;
		VectorB *B = (VectorB*) params;
		//if(z>1.-B->gen) return _integrand_q(1.-B->gen,params);
		double xi(B->xi);
		double mu(B->mu);
		double x(c*xi);
		int n_f(B->n_f);

		double q_z(B->phiz(lha_(B->q),x,mu));
		double q_xi(B->phixi[B->q]);
		double g_z(B->phiz(lha_(e_g),x,mu));
		double g_xi(B->phixi[e_g]);
		//cout << "VectorB::_integrand_q (" << z << ")" << endl
		//	 << "xi="<<xi<<" mu="<<mu<<" x="<<x<<" n_f="<<n_f<< endl
		//	 << "q_z="<<q_z<<" q_xi="<<q_xi<<" g_z="<<g_z<<" g_xi="<<g_xi<< endl;



	#ifdef ONLY_MIXED
		if(B->prod==PROD_H)
			return (B->P_iq(z, xi, c*q_z,q_xi,n_f));
		else
			return (B->P_ig(z, xi, c*g_z,g_xi,n_f));
	#endif
		//cout << "_integrand_q(" <<z<<") = "<<(B->P_iq(z, xi, c*q_z,q_xi,n_f)	+ B->P_ig(z, xi, c*g_z,g_xi,n_f)) << endl;
		return (B->P_iq(z, xi, c*q_z,q_xi,n_f)	+ B->P_ig(z, xi, c*g_z,g_xi,n_f));
	}



	inline double DistOtimesPhi(Integrand f){
		double DistOtimesPhi = 0.0;
		int intMem=500000;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (intMem);
		double result, error;
		int status;

		gsl_function k;
		k.function = f;
		k.params = (void*)(this);
		//status = gsl_integration_qags (&k, xi, 1.0, 0., gen, intMem, w, &result, &error);
		status = gsl_integration_qag (&k, xi, 1.0, 0., gen, intMem,6, w, &result, &error);
		DistOtimesPhi=result;

		if(status!=GSL_SUCCESS && status!=GSL_EROUND){
			cout << "-----------Vector B DistOtimesPhi_qag went wrong "<<endl
					<< "result = " << result << "\tErr = " << error/result << endl
					<< " f_"<<lha_(q)<<"("<<xi<<")= " << k.function(xi,k.params)
					<< "\tf(1.0)= " << k.function(1.,k.params) << endl;
			GSL_ERROR_CODE(status);
			//cout << *this << endl;
			cout << "----------------------------------------" << endl;
		}

		gsl_integration_workspace_free (w);
		return (DistOtimesPhi);
	}

	inline double Bi_a0_L0_e0(e_parton i){
		return phixi[i];
	}

	inline double Bq_a1_L1_e1(){
		//cout << "Bq_a1_L1_e1" << endl;
		P_ig=&P_qgf;
		P_iq=&P_qqf;
		double Bq_a1_L1_e1=DistOtimesPhi(_integrand_q);
		Bq_a1_L1_e1*=(-0.5);
		return Bq_a1_L1_e1;
	}

	inline double Bg_a1_L1_e1(){
		P_ig=&P_ggf;
		P_iq=&P_gqf;
		double Bg_a1_L1_e1=DistOtimesPhi(_integrand_g);
		Bg_a1_L1_e1*=(-0.5);
		return Bg_a1_L1_e1;
	}

	inline double Bq_a1_L0_e2(){
		//cout << "Bq_a1_L0_e2" << endl;
		P_ig=&R_qgf;
		P_iq=&R_qqf;
		double Bq_a1_L0_e2=DistOtimesPhi(_integrand_q);
		return Bq_a1_L0_e2;
	}

	inline double Bg_a1_L0_e2(){
		P_ig=&R_ggf;
		P_iq=&R_gqf;
		double Bg_a1_L0_e2=DistOtimesPhi(_integrand_g);
		return Bg_a1_L0_e2;
	}

	inline double Bq_a2_L2_e2(){
		//cout << "Bq_a2_L2_e2" << endl;
		P_ig=&D_qgf;
		P_iq=&D_qqf;
		double Bq_a2_L2_e2=DistOtimesPhi(_integrand_q);
		Bq_a2_L2_e2*= 0.125;
		Bq_a2_L2_e2+= 0.5*beta_0(n_f)*_V[q][1][1]/a;
		return Bq_a2_L2_e2;
	}

	inline double Bg_a2_L2_e2(){
		P_ig=&D_ggf;
		P_iq=&D_gqf;
		double Bg_a2_L2_e2=DistOtimesPhi(_integrand_g);
		Bg_a2_L2_e2*= 0.125;
		Bg_a2_L2_e2+= 0.5*beta_0(n_f)*_V[e_g][1][1]/a;
		return Bg_a2_L2_e2;
	}

	inline void Bq(){
		//cout << "VectorB::Bq()\tn_f=" << n_f << " Phixi = " << *phixi << endl;
		for(q=e_(-n_f); q<=e_(n_f); q++){
			if(q==e_g)continue;
			_V[q] = Pol< Pol<double> >(order,Pol<double>(0));

			_V[q][0][0]=Bi_a0_L0_e0(q);
			if(order>0){
				_V[q][1].order(1);
				//if(!improved){
					_V[q][1][0]=Bq_a1_L0_e2()*a;
					_V[q][1][1]=Bq_a1_L1_e1()*a;
				//}else{
				//	_V[q][1][1]=Bq_a1_L1_e1()*a;
				//}
				if(order>1){
					_V[q][2].order(2);
					//if(!improved){
						_V[q][2][2]=Bq_a2_L2_e2()*gsl_pow_2(a);
					//}else{
					//	_V[q][2][0]=Bq_a1_L0_e2()*a;
					//	_V[q][2][2]=Bq_a2_L2_e2()*gsl_pow_2(a);
	}	}	}	}	//}

	inline void Bg(){
		_V[e_g] = Pol< Pol<double> >(order,Pol<double>(0));

		_V[e_g][0].order(0);
		_V[e_g][0][0]=Bi_a0_L0_e0(e_g);

		if(order>0){
			_V[e_g][1].order(1);
			//if(!improved){
				_V[e_g][1][0]=Bg_a1_L0_e2()*a;
				_V[e_g][1][1]=Bg_a1_L1_e1()*a;
			//}else{
			//	_V[e_g][1][1]=Bg_a1_L1_e1()*a;
			//}
			if(order>1){
				_V[e_g][2].order(2);
				//if(!improved){
					_V[e_g][2][2]=Bg_a2_L2_e2()*gsl_pow_2(a);
				//}else{
				//	_V[e_g][2][0]=Bg_a1_L0_e2()*a;
				//	_V[e_g][2][2]=Bg_a2_L2_e2()*gsl_pow_2(a);
				//}
			}
		}
	}

	inline void update(){
		if(prod==PROD_H) Bg();
		else Bq();
	}

public:
	inline void update(double xi_,double mu_, double a_, uint n_f_){
		//cout << "VectorB::update("<<xi_<<","<<mu_<<","<<a_<<","<<n_f_<<")" << endl;
		xi=xi_; a=a_; mu=mu_; n_f=n_f_;
		phixi.update(xi_,mu_);
		update();
	}

	inline VectorB(){}
	inline VectorB(CS_params const &p, int const sign_=1)
	:Vector< N_PARTONS,Pol< Pol<double> > >(Pol< Pol<double> >(0,Pol<double>(0)))
	,prod(p.prod)
	,order(p.order)
	,sign(sign_)
	,phixi(sign)
	,phiz(sign)
	,gen(p.gen){
	}
};



class BQB: public Pol< Pol<double> >{
protected:
	CS_params* p;
	Q_Matrix _Q;
	double mu;
	double a;
	VectorB _B[2];
	uint n_f;

	inline virtual void v_update(){
		_B[0].update(p->xi[0],mu,a,n_f);
		_B[1].update(p->xi[1],mu,a,n_f);


		Pol< Pol<double> >tmp(_B[0]*_Q*_B[1]);
		tmp.order(_order);
		for(uint i=0; i<=_order;++i)_P[i]=tmp[i];
#ifdef ONLY_MIXED
		_P[0][0]=0.;
#endif
	}

	inline virtual void v_init(){
	}

	inline virtual void v_init(CS_params& p_){
	}

public:

	inline void update(){
		//cout << "BQB::update("<<p->y<<"["<<p->xi[0]<<","<<p->xi[1]<<"])" << endl;
		v_update();
		if(p->improved) Pol_a_L_Improve(*this);
		order(p->order);
	}

	inline void update(double y){
		p->update(y);
		update();
	}

	inline void init(){
		mu=p->mu;
		a=p->a;
		n_f=p->n_f;
		v_init();
	}

	inline void init(double qT_, double mu_){
		p->update(qT_,mu_);
		init();
	}


	inline void init(CS_params &p_){
		p=(&p_);
		_Q = Q_Matrix(p->prod,p->ckm);
		_B[0]= VectorB(*p,p->sign[0]);
		_B[1]= VectorB(*p,p->sign[1]);
		order(p->order);
		v_init(p_);
	}

	inline BQB(CS_params &p_){
		init(p_);
	}


	inline Pol< Pol<double> > & operator() (double y){
		update(y);
		return *this;
	}
	inline Pol< Pol<double> > & operator() (){
		update();
		return *this;
	}

	inline virtual ~BQB(){}
};

class BQB_fast:public BQB{
protected:
	VectorPhi _phixi[2];
	VectorPhi _phiz[2];


	double _xi[2];
	double gen;

	DistFunc P_iq;
	DistFunc P_ig;

	//double div;


	inline double PhiOtimesDIST_Q_DistOtimesPhi(){
		//cout << "->BQB_fast::PhiOtimesDIST_Q_DistOtimesPhi()" << endl;
		double DistOtimesPhi = 0.0;

		for(int i=0; i<=1; ++i){
			_B[i].phixi.update(_xi[i],mu);
			_B[i].P_ig=P_ig;
			_B[i].P_iq=P_iq;
			_B[i].mu=mu;
			_B[i].xi=_xi[i];
			_B[i].n_f=n_f;
		}

		if(p->prod==PROD_H){
			//cout << "\tPROD_H" << endl;
			Integrand f = &VectorB::_integrand_g;
			DistOtimesPhi= _B[0].DistOtimesPhi(f)*_Q[e_g][e_g]*_B[1].DistOtimesPhi(f);
			//_B[0].update(_xi[0],mu,a,n_f);
		}else{
			Integrand f = &VectorB::_integrand_q;
			for(e_parton q=e_(n_f); q>e_g; --q){
				_B[0].q=q;	_B[1].q=bar(q);
				DistOtimesPhi+= _B[0].DistOtimesPhi(f)*_Q[q][bar(q)]*_B[1].DistOtimesPhi(f);
				_B[0].q=bar(q);	_B[1].q=q;
				DistOtimesPhi+= _B[0].DistOtimesPhi(f)*_Q[bar(q)][q]*_B[1].DistOtimesPhi(f);
			}
		}

		//cout << "<-BQB_fast::PhiOtimesDIST_Q_DistOtimesPhi()" << endl;
		return (DistOtimesPhi);
	}

	inline static double _integrand_q(double z, void* params){
		double erg=0.0;
		BQB_fast *P = (BQB_fast*) params;
		//if(z>1.-P->gen) return _integrand_q(1.-P->gen,params);
		double c=1./z;
		double mu=P->mu;
		uint n_f=P->n_f;
		double BB00 = (*P)[0][0];
		VectorPhi *phixi = P->_phixi;
		VectorPhi *phiz = P->_phiz;
		double *xi = P->_xi;

		for(int i=0; i<=1; ++i){
			int j=1-i%2;
			if(z>xi[i]){
				phiz[i].update(c*xi[i],mu);
				double sQPj=phixi[j].sumofQphi(P->_Q,n_f);
#ifndef ONLY_MIXED
				erg+=P->P_iq(z,xi[i],c*(phiz[i]*P->_Q*phixi[j]),BB00,n_f);
#endif
				erg+=P->P_ig(z,xi[i],c*(phiz[i](lha_g)*sQPj),phixi[i](lha_g)*sQPj,n_f);
			}
		}
		return erg;
	}

	inline static double _integrand_g(double z, void* params){
		double erg=0.0;
		BQB_fast *P = (BQB_fast*) params;

		double c=1./z;
		double mu=P->mu;
		uint n_f=P->n_f;
		double BB00 = (*P)[0][0];
		VectorPhi *phixi = P->_phixi;
		VectorPhi *phiz = P->_phiz;
		double *xi = P->_xi;

		for(int i=0; i<=1; ++i){
			int j=1-i%2;
			if(z>xi[i]){
				phiz[i].update(c*xi[i],mu);
				double g_i =phiz[i][e_g];
				double g_j = phixi[j][e_g];
				double soq_z = phiz[i].sumofquarks(n_f);
				double soq_xi = phixi[i].sumofquarks(n_f);
				erg+=P->P_iq(z,xi[i],c*soq_z*g_j,soq_xi*g_j,n_f);
#ifndef ONLY_MIXED
				erg+=P->P_ig(z,xi[i],c*(g_i*g_j),BB00,n_f);
#endif
			}
		}

		return erg;
	}

	inline double PhiOtimesDEL_Q_DistOtimesPhi(Integrand f){
		double erg = 0.0;
		int intMem=900000;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (intMem);
		double result, error;
		int status;

		gsl_function k;
		k.function = f;
		k.params = (void*)(this);
		//status = gsl_integration_qags (&k, GSL_MIN(_xi[0],_xi[1]), GSL_MAX(_xi[0],_xi[1]), 0., gen, intMem, w, &result, &error);
		status = gsl_integration_qag (&k, GSL_MIN(_xi[0],_xi[1]), GSL_MAX(_xi[0],_xi[1]), 0., gen, intMem,6, w, &result, &error);
		erg=result;
		if(status!=GSL_SUCCESS && status!=GSL_EROUND){
					cout << "-----------BQB_Z:: DistOtimesDel_DistOtimesPhi_qag went wrong y="<< p->y<<endl
							<< "result = " << result << "\tErr = " << error/result << "\tgen = " << gen << endl
							<< " f("<<GSL_MIN(_xi[0],_xi[1])<<")= " << k.function(GSL_MIN(_xi[0],_xi[1]),k.params)
							<< "\tf("<<GSL_MAX(_xi[0],_xi[1])<<")= " << k.function(GSL_MAX(_xi[0],_xi[1]),k.params) << endl;
					GSL_ERROR_CODE(status);
					cout << "----------------------------------------" << endl;
				}
		/*div=k.function(1.,k.params);
		 cout << div << " = ";
		if(div!=div){
			cout << "NAN" << endl;
		}else{
			cout << endl;
		}*/
		//status = gsl_integration_qags (&k, GSL_MAX(_xi[0],_xi[1]), 1., 0., gen, intMem, w, &result, &error);
		status = gsl_integration_qag (&k, GSL_MAX(_xi[0],_xi[1]),1. , 0., gen, intMem,6, w, &result, &error);
		erg+=result;
		if(status!=GSL_SUCCESS && status!=GSL_EROUND){
			cout << "-----------BQB_Z:: DistOtimesDel_DistOtimesPhi_qag went wrong y="<< p->y<<endl
					<< "result = " << result << "\tErr = " << error/result << "\tgen = " << gen << endl
					<< " f("<<GSL_MAX(_xi[0],_xi[1])<<")= " << k.function(GSL_MAX(_xi[0],_xi[1]),k.params)
					<< "\tf(1.)= " << k.function(1.,k.params) << endl;
			GSL_ERROR_CODE(status);
			cout << "----------------------------------------" << endl;
		}

		gsl_integration_workspace_free (w);
		//div=0.;
		return (erg);
	}

	inline void update_Z(){
		_phixi[0].update(_xi[0],mu);
		_phixi[1].update(_xi[1],mu);

		Integrand f= &_integrand_q;

		_P[0].order(0);
#ifndef ONLY_MIXED
		_P[0][0]=_phixi[0]*_Q*_phixi[1];
#endif
		if(!_order) return;
		_P[1].order(1);
		P_iq=&R_qqf;	P_ig=&R_qgf;	_P[1][0]= PhiOtimesDEL_Q_DistOtimesPhi(f);
		P_iq=&P_qqf;	P_ig=&P_qgf;	_P[1][1]= -0.5*PhiOtimesDEL_Q_DistOtimesPhi(f);
		if(_order>1 ){
			_P[2].order(2);
			P_iq=&D_qqf;	P_ig=&D_qgf;
			_P[2][2]=0.125*(PhiOtimesDEL_Q_DistOtimesPhi(f)-2.*beta_0(n_f)*(-2.*_P[1][1]));
			P_iq=&P_qqf;	P_ig=&P_qgf;	_P[2][2]+=0.25*PhiOtimesDIST_Q_DistOtimesPhi();
		}
		setX(a);
	}

	inline void update_H(){
		//cout << "->BQB_fast::update_H()" << endl;
		_phixi[0].update(_xi[0],mu);
		_phixi[1].update(_xi[1],mu);

		Integrand f= &_integrand_g;

		_P[0].order(0);
#ifndef ONLY_MIXED
		_P[0][0]=_phixi[0]*_Q*_phixi[1];
#endif
		if(!_order) return;
		_P[1].order(1);
		P_iq=&R_gqf;	P_ig=&R_ggf;	_P[1][0]= PhiOtimesDEL_Q_DistOtimesPhi(f);
		P_iq=&P_gqf;	P_ig=&P_ggf;	_P[1][1]= -0.5*PhiOtimesDEL_Q_DistOtimesPhi(f);
		if(_order>1){
			_P[2].order(2);
			P_iq=&D_gqf;	P_ig=&D_ggf;
			//debug.clear();
			_P[2][2]=0.125*(PhiOtimesDEL_Q_DistOtimesPhi(f)-2.*beta_0(n_f)*(-2.*_P[1][1]));
			//debug[0].push_back(p->y);
			//debug[1].push_back(_P[2][2]);
			P_iq=&P_gqf;	P_ig=&P_ggf;	_P[2][2]+=0.25*PhiOtimesDIST_Q_DistOtimesPhi();
		}
		setX(a);
		//cout << "<-BQB_fast::update_H()" << endl;
	}

	inline void v_update(){
		_xi[0]=p->xi[0];
		_xi[1]=p->xi[1];
		if(p->prod==PROD_H) update_H();
		else update_Z();
	}

	inline void v_init(){
	}

	inline void v_init(CS_params &p_){
		gen=p->gen;
		_phixi[0]=VectorPhi(p->sign[0]);
		_phixi[1]=VectorPhi(p->sign[1]);
		_phiz[0]=VectorPhi(p->sign[0]);
		_phiz[1]=VectorPhi(p->sign[1]);
	}



public:
	inline BQB_fast(CS_params &p_):BQB(p_){v_init(p_);}
};

#endif
