#ifndef RESUMED_H
#define RESUMED_H

#include "gPDF.h"
#include <gsl/gsl_sf_bessel.h>



class Epol{
protected:
	CS_params* p;
	Pol< Pol< Pol<double> > > Eexp;
	Pol< complex<double> > _k;
	Pol< Pol<double> > _g;
	uint n_f;
	double eta;
	double a;

	uint order;

	inline void init(){
			//cout<<"Epol::init()\tCB="<<p->C_B<<"\thB="<<p->h_B<<"\tn_f="<<p->n_f <<endl;
			n_f=(p->n_f);
			order=(p->order);
			a=0.;
			eta=0.;


			//Eexp[eta][a][L]
			Eexp.order(1);
			Eexp[0].order(3);
			Eexp[1].order(3);

			Eexp[0][0].order(0);
			Eexp[1][0].order(1);
			Eexp[0][0][0]=0.;
			Eexp[1][0][0]=0.;
			Eexp[1][0][1]=-1.;		//\in \hat{k}

			Eexp[0][1].order(2);
			Eexp[1][1].order(2);
			Eexp[0][1][0]=0.;
			Eexp[0][1][1]=2.*p->h_B;		//\in \hat{g}
			Eexp[0][1][2]=-0.5*Gamma_0;		//\in \hat{k}
			Eexp[1][1][0]=-d_2(n_f);		//\in \hat{g}
			Eexp[1][1][1]=-Gamma_1_0(n_f);		//\in \hat{g}
			Eexp[1][1][2]=-0.5*beta_0(n_f);		//\in \hat{k}

			Eexp[0][2].order(3);
			Eexp[1][2].order(3);
			Eexp[0][2][0]=0.;
			Eexp[0][2][1]=0.;
			Eexp[0][2][2]=-(0.5*Gamma_1(n_f)-beta_0(n_f)*p->h_B);		//\in \hat{g}
			Eexp[0][2][3]=-beta_0(n_f)*Gamma_0/3.;		//\in \hat{g}
			Eexp[1][2][0]=0.;
			Eexp[1][2][1]=0.;
			Eexp[1][2][2]=-(beta_0(n_f)*Gamma_1_0(n_f)+beta_1(n_f)/2.);		//\in \hat{g}
			Eexp[1][2][3]=-gsl_pow_2(beta_0(n_f))/3.;		//\in \hat{g}

			Eexp[0][3].order(4);
			Eexp[1][3].order(4);
			Eexp[0][3][0]=0.;
			Eexp[0][3][1]=0.;
			Eexp[0][3][2]=0.;
			Eexp[0][3][3]=0.;
			Eexp[0][3][4]=-Gamma_0*gsl_pow_2(beta_0(n_f))*0.25;		//\in \hat{g}
			//cout << "Eexp[0][3][4] = "<<Eexp[0][3][4]<< endl;
			Eexp[1][3][0]=0.;
			Eexp[1][3][1]=0.;
			Eexp[1][3][2]=0.;
			Eexp[1][3][3]=0.;
			Eexp[1][3][4]=-gsl_pow_3(beta_0(n_f))*0.25;		//\in \hat{g}
			//cout << "Eexp[1][3][4] = "<<Eexp[1][3][4]<< endl;

			Eexp[0]*=p->C_B;
			Eexp[1]*=p->C_B;
		}

public:

	inline void init(CS_params &p_){
		p=(&p_);
		init();
	}

	inline void update(){
		if(n_f!=p->n_f) init();
		if(eta!=p->eta || a!=p->a){
			eta=(p->eta);
			a=(p->a);

			if(p->res==RES_LL){
				_g = Eexp.von(eta);
				_g.setX(a);
				_k=Pol< complex<double> >(2,0.);
				_k[1]=_g[0][1]; _g[0][1]=0.;
				_k[2]=_g[1][2]; _g[1][2]=0.;
				if(p->improved) Pol_a_L_Improve(_g);
				_g.order(order);
				if(p->exponent){
					_k += _g.von(1.);
					_g.order(0);
					_g[0].order(0);
					_g[0][0]=0.;
				}
			}else{
				_g=(Eexp[1].timesX()*eta);
				_g[1][1]=0.;
				_g+=Eexp[0];
				_g.order(order);
				_g.setX(a);
				_k.order(0);
				_k[0]=0.;
			}
			//cout	<< "Epol::update()\t eta="<<eta<<"\ta="<<a<<"\tCB="<<p->C_B<<"\thB="<<p->h_B<<"\tn_f="<<n_f <<endl
			//		<< "\t g\t=" <<_g << endl;
		}
	}

	inline void update(double qT_, double mu_){
		p->update(qT_,mu_);
		update();
	}

	inline Epol(CS_params &p_){ init(p_);}

	inline Pol< complex< double> > const k(){
		update();
		return _k;
	}
	inline Pol< Pol< double> > const g(){
		update();
		return _g;
	}
	inline Pol< Pol< double> > const E(){
		return exp(g());
	}

	inline Pol< Pol<double> > const operator [](uint i)const{ return Eexp[i];}
};

extern std::ostream& operator<<( std::ostream& os, const Epol& E );

class K_Integral{
protected:
	CS_params* p;

	Epol* _E;
	uint n_f;
	double M;

	uint N;
	double qT;
	double mu;
	double a;
	double eta;


	/*inline virtual void v_init(){
		cout << "*************************************************************************" << endl
			 << "**              class K_Integral::v_init not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}*/
	inline virtual void v_init(CS_params &p_){
		cout << "*************************************************************************" << endl
			 << "**              class K_Integral::v_init not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}
	inline virtual void v_update(){
		cout << "*************************************************************************" << endl
			 << "**              class K_Integral::init not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}


public:
	inline Pol<complex< double > > const k(){ return _E->k();}
	inline Pol<Pol< double > > const g(){ return _E->g();}
	inline Pol<Pol< double > > const E(){ return _E->E();}

	/*inline void init(){
		_E->init();
		v_init();
	}*/

	inline void init(CS_params &p_){
		p=(&p_);
		delete _E;
		_E=new Epol(p_);
		M=p->M;
		v_init(p_);
	}

	inline void update(){
		//cout << "K_Integral::update()" << endl;
		n_f=p->n_f;
		qT=(p->qT);
		mu=p->mu;
		a=p->a;
		eta=p->eta;
		v_update();
	}

	inline void update(double qT_, double mu_){
		p->update(qT_,mu_);
		update();
	}

	inline virtual Pol<double> const replaceL(Pol< Pol<double> > const &P)const{
		cout << "*************************************************************************" << endl
			 << "**              class K_Integral::replaceL not defined -> exit(1)      **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
		return Pol< double >(0,0.);
	}

	inline Pol<double> const replaceL_E(Pol< Pol<double> > const &P){
		Pol<Pol<double> > tmp(P*E());
		tmp.order(p->order);
		return replaceL(tmp);
	}

	inline virtual void print( std::ostream &os)const{
		cout << "*************************************************************************" << endl
			 << "**              class K_Integral::print not defined -> exit(1)         **" << endl
			 << "*************************************************************************" << endl;
		exit(1);
	}

	inline K_Integral(CS_params &p_):_E(NULL){
	}

	inline virtual ~K_Integral(){
		delete _E;
	}
};

extern std::ostream& operator<<( std::ostream& os, const K_Integral* K );


class K_LL:public K_Integral{
protected:
	std::vector< double > _P;

	double Lambda;
	CUTOFF cutoff;
	double gen;

	uint n;
	Pol< complex<double> > _k_exp;

	inline static double k_n(double x, void* params){
		K_LL *K = (K_LL*)params;
		int n = K->n;
		Pol< complex<double> > k_exp = K->_k_exp;

		double L_T=2*log(K->mu*x/b_0);

		double _kT =x;
		_kT*=0.5*gsl_sf_bessel_J0(K->qT*x);
		if(K->cutoff == CUTOFF_GAUSS){
			_kT*=exp(-2*gsl_pow_2(K->Lambda*x));
		}else if(K->cutoff == CUTOFF_DIPOL){
			_kT*=gsl_pow_2(1./(1.+gsl_pow_2(K->Lambda*x)));
		}
		k_exp.von(L_T);
		_kT*=exp(k_exp.von(L_T).real());

		if(n)_kT*=gsl_pow_int(L_T,n);
		if(gsl_finite(_kT))	return _kT;
		else return 0.0;
	}

	inline static double kwick_n(double z, void* params){
		K_LL *K = (K_LL*)params;
		Pol< complex<double> > k_exp = K->_k_exp;
		complex< double > n(K->n,0.);

		complex< double > L_T(2*log(K->mu*z/b_0),PI);

		double _kT =z;
		_kT*=(exp(k_exp.von(L_T)+log(L_T)*n)).imag();
		_kT*=(-1.0*M_1_PI*gsl_sf_bessel_K0(K->qT*z));
		return _kT;
	}

	inline double K_n(Integrand f){
		double K_n = 0.0;
		int intMem=500000;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (intMem);
		double result, error;
		int status;

		gsl_function k;
		k.function = f;
		k.params = (void*)(this);
		//status = gsl_integration_qags (&k, params.xi, 1.0, params.gen, params.gen, intMem, w, &result, &error);
		//status = gsl_integration_qag (&k, xi, 1.0, 0., gen, intMem,6, w, &result, &error);
		if(cutoff == CUTOFF_HARD) {
			status = gsl_integration_qag (&k, 0.0, 1.0/Lambda,  0.0, gen, intMem, 6, w, &result, &error);
		}else{
			status = gsl_integration_qagiu (&k, 0.0, 0.0, gen, intMem, w, &result, &error);
		}
		K_n=result;

		if(status!=GSL_SUCCESS){
			cout << "-----------K_LL K_" << n << " went wrong "<<endl
					<< "result = " << result << "\tErr = " << error/result << endl
					<< " f(xmin)= " << k.function(0.,k.params)
					<< "\tf(xmax)= " << k.function(100.,k.params) << endl;
			GSL_ERROR_CODE(status);
			cout << "----------------------------------------" << endl;
		}

		gsl_integration_workspace_free (w);
		return K_n;
	}

	inline void v_update(){
		Integrand k;
		_k_exp = _E->k();
		if(cutoff==CUTOFF_NO) k=&kwick_n;
		else k=&k_n;
		for(n=0; n<=N; ++n) _P[n]=K_n(k);
	}

	/*inline void v_init(){
	}*/

	inline void v_init(CS_params &p_){
		N=p->order*(2+p->improved);
		_P.resize(N+1);
		Lambda=p->Lambda;
		cutoff=p->cutoff;
		gen=p->gen;
	}

public:

	inline virtual Pol<double> const replaceL(Pol< Pol<double> > const &P)const{
		uint n=P.order();
		Pol< double > tmp(n,0.);
		for(uint i=0; i<=P.order(); ++i){
			//uint m = GSL_MIN_INT(_N,P[i].order());
			uint m = P[i].order();
			for(uint j=0; j<=m; ++j){
				tmp[i]+=P[i][j]*_P[j];
			}
		}
		return tmp;
	}

	inline void print( std::ostream &os)const{
		os << "( " << _P[0];
		for(uint i=1; i<=N; ++i) os << " , " << _P[i];
		os << " )";
	}
	inline K_LL(CS_params &p_):K_Integral(p_){
		init(p_);}

};

class K_LLO:public K_Integral{
protected:

	vector< Pol<double> > _Kla;
	double Lmu;

	inline void v_update(){
		Lmu=p->L_mu;
		//cout << "K_LLO::update() \t_qT=" << _qT << "\t_eta=" <<  _eta
		//	 <<	"\t_Lmu=" << _Lmu << "\t_LM=" << _LM << endl;
		if(p->prod==PROD_H) eta*=C_A;
		else eta*=C_F;
		//double _Lmu=(2.*log(_mu/_qT));
		//double _LM=(2*log(_M/_qT));
		_Kla[0][0]=0.;
		_Kla[0][1]=1.;
		_Kla[0][2]=-Lmu;
		_Kla[0][3]=0.5*gsl_pow_2(Lmu);
		_Kla[0][4]=(-gsl_pow_3(Lmu)+4.*zeta_3)/6.;
		_Kla[0][5]=(gsl_pow_4(Lmu)-16.*Lmu*zeta_3)/24.;
		_Kla[0][6]=(-gsl_pow_5(Lmu)- 2.*(-24.*zeta_5) + 40.*gsl_pow_2(Lmu)*zeta_3)/120.;
		for(uint n=1; n<=N; ++n) for(uint m=0; m<N; ++m) _Kla[n][m]=-(m+1.)*_Kla[n-1][m+1];

		double tmp = 1./gsl_pow_2(qT);
		for(uint m=0; m<N; ++m){
			for(uint n=0; n<=N; ++n) _Kla[n][m]*=tmp;
			tmp*=eta;
		}
	}

	/*inline void v_init(){
	}*/

	inline void v_init(CS_params &p_){
		N=(7);
		_Kla.resize(N+1,Pol<double>(N,0.));
	}


public:

	inline virtual Pol<double> const replaceL(Pol< Pol<double> > const &P)const{
		uint order=P.order();
		Pol<double> tmp(order,0.);
		for(uint n=0; n<=order; ++n){
			if(P[n].order()>=N){
				cout << "**************************************************************" << endl
					 << "**                 L_\\perp^n : n>K_LLO::_N                **" << endl
					 << "**class K_LLO:: please extend update() function -> exit(1)**" << endl
					 << "**************************************************************" << endl;
			}
			for(uint a=0; a<=n;++a)for(uint l=0; l<=P[a].order(); ++l){
			tmp[n]+=P[a][l]*_Kla[l][n-a];
			//cout << "tmp(n="<<n<<",a="<<a<<",l="<<l<<",n-a="<<n-a<<")\t=" << tmp << endl;
			}
		}
		return tmp;
	}


	inline void print( std::ostream &os)const{
		os << "( " << _Kla[0];
		for(uint i=1; i<N; ++i) os << " , " << _Kla[i];
		os << " )";
	}

	inline K_LLO(CS_params &p_):K_Integral(p_){
		init(p_);	}
};


#endif
