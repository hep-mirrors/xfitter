#ifndef POLYNOM_h
#define POLYNOM_h

#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

/*#ifndef uint
#define uint unsigned int
#endif*/
typedef unsigned int uint;

	inline uint fac(uint n){
		if(n<=1) return 1;
		return n*fac(n-1);
	}

	inline double bin(int n, int k){
		return ((double)fac(n)/((double)(fac(k)*fac(n-k))));
	}

template < typename T>
class Pol;


template < typename T >
Pol< T > const exp(Pol< T > P);

template < typename T>
std::ostream& operator<<( std::ostream& os, Pol<T> const &P );


template < typename T>
class Pol {
protected:
	uint _order;
	std::vector< T > _P;

public:
//----------------------          Constructor
	inline Pol< T > ():_order(0),_P(1){
		_P.resize(1);
	}

	inline Pol< T > (T const &t):_order(0),_P(1,t){
		//std::cout << "Pol< T > (T const &t) = " << _P[0] << std::endl;
	}

	/*inline Pol< T > (uint const order):_order(order),_P(order+1){
		_P.resize(order+1);
	}*/

	inline Pol< T > (uint const order, T const &t):_order(order),_P(order+1,t){
		_P.resize(order+1);
	}

	inline Pol< T > (Pol<T> const &A):_order(A.order()),_P(A._P){}

//----------------------          get
	inline uint const order() const {return _order;}
	inline T const operator ()(uint const i)const{ return ((_order>=i) ? _P[i] : T()) ;}
//----------------------          set
	inline void order(uint const order_){
		_order=order_;
		_P.resize(order_+1);
	}
	inline T& operator [] (uint const i){ return _P[i];}
	inline T const&  operator [] (uint const i)const {return _P[i];	}

/*	// ---------------         skalar Mult

	//template < typename T2 >
	inline Pol< T >& operator *= (T const &a){
		std::cout << "***************************\tP*S" << std::endl;
		for(uint i=0; i<=_order; ++i) _P[i]*=a;
		return *this;
	}

	//template < typename T2 >
	inline Pol< T >& operator /= (T const &a){
		for(uint i=0; i<=_order; ++i) _P[i]/=a;
		return *this;
	}

	//template < typename T2 >
	inline Pol< T >const operator * (T const &a)const{
		Pol< T > tmp(*this);
		tmp*=a;
		return tmp;
	}

	//template < typename T2 >
	inline Pol< T >const operator / (T const &a)const{
		Pol< T > tmp(*this);
		tmp/=a;
		return tmp;
	}*/

	inline Pol< T > const operator -()const{
		Pol< T > tmp(*this);
		for(uint i=0; i<=_order; ++i) tmp[i]= -tmp[i];
		return tmp;
	}

	// ---------------         Pol addition

	template < typename T2 >
	inline Pol< T >& operator += (Pol< T2 > const &A){
		for(uint i=0; i<=_order; ++i) _P[i]+=A(i);
		for(uint i=(_order+1); i<=A.order(); ++i) _P.push_back(A(i));
		_order=(_P.size()-1);
		return *this;
	}
	template < typename T2 >
	inline Pol< T >& operator -= (Pol< T2 >  const &A){
		for(uint i=0; i<=_order; ++i) _P[i]-=A(i);
		for(uint i=(_order+1); i<=A.order(); ++i) _P.push_back(-A(i));
		_order=(_P.size()-1);
		return *this;
	}

	template < typename T2 >
	inline Pol< T >const  operator + (Pol< T2 > A)const{
		A+=(*this);
		return A;
	}
	template < typename T2 >
	inline Pol< T >const operator - (Pol< T2 > A)const{
		A-=(*this);
		return (-A);
	}

/*	// ---------------         skalar Addition
	inline Pol< T >& operator += (T const &a){
		_P[0]+=a;
		return *this;
	}
	inline Pol< T >& operator -= (T const &a){
		_P[0]-=a;
		return *this;
	}
	inline Pol< T >const operator + (T const &a)const{
		Pol< T > tmp(*this);
		tmp+=a;
		return tmp;
	}
	inline Pol< T >const operator - (T const &a)const{
		Pol< T > tmp(*this);
		tmp-=a;
		return tmp;
	}*/

	// ---------------         Pol Mult


	template < typename T2 >
	inline Pol< T >& mul(Pol< T2 > const P,uint order_){
		Pol<T > tmp(*this);
		order(order_);
		for(uint i=0; i<=_order; ++i){
			_P[i]= tmp(0)*P(i-0);
			for(uint k=1; k<=i; ++k) {
				_P[i]+= tmp(k)*P(i-k);
			}
		}
		return *this;
	}

	template < typename T2 >
	inline Pol< T >& operator *=(T2 const &t){
		//std::cout << "***************************\tP*P("<< t << ")" << std::endl;
		Pol P(t);
		mul(P,_order+P.order());
		return *this;
	}

	template < typename T2 >
	inline Pol< T >const operator *(T2 const &t)const{
		Pol< T > tmp(*this);
		tmp*=t;
		return tmp;
	}

	/*template < typename T2 >
	inline Pol< T >& operator *=(Pol< T2 > const &P){
		std::cout << "***************************\tP*P("<< P << ")" << std::endl;
				mul(P,_order+P.order());
		return *this;
	}

	template < typename T2 >
	inline Pol< T >const operator *(Pol< T2 > const &P)const{
		Pol< T > tmp(*this);
		tmp.mul(P,_order+P.order());
		return tmp;
	}*/


	// --------------------  von

	template < typename T2 >
	inline T const von (T2 const &x, uint const order)const{
		int m = ((_order>=order) ? order : _order) ;
		T von = _P[m];
		for(int i=(m-1); i>=0; --i){
			von*=x;
			von+=_P[i];
			//if(i==3)cout <<"von + " << x << " * " << _P[i] << endl;
		}
		return von;
	}

	template < typename T2 >
	inline T const von (T2 const &x)const{
		return von(x,_order);
	}

	template < typename T2 >
	inline Pol< T >& setX(T2 const &x){
		T X(x);
		for(uint i=1; i<=_order; ++i){
			//cout <<"setx " << x << " * " << _P[i] << endl;
			_P[i]*=X;
			X*=x;
		}
		return *this;
	}

	inline Pol< T > & timesX(){
		order(_order+1);
		for(uint i=_order; i>0; --i) _P[i]=_P[i-1];
		_P[0]-=_P[0];
		return *this;
	}


// --------------------  exp

	inline Pol< T > & _exp(){
		*this = exp(*this);
		return *this;
	}

	//---------------   PolImRe
/*	inline void ReIm(T const &ix, Pol< T > &R, Pol< T > &I, uint const order_)const{
		R.order(order_);
		I.order(order_);
		T X =-ix*ix;
		double f;
		T XR;
		for(uint n=0; n<=order_;++n){
			XR=X;
			R[n]=(*this).operator ()(n);
			I[n]=(*this).operator ()(n+1)*(n+1)*ix;
			for(uint l=1; (n+2*l)<=order_; l++){
				f=bin(n+2*l,n);
				R[n]+=(*this).operator ()(n+2*l)*XR*f;
				I[n]+=(*this).operator ()(n+2*l+1)*XR*ix*f*((double)(n+2*l+1))/((double)(2*l+1));
				XR*=X;
			}
		}
	}

	inline void ReIm(T const &ix, Pol< T > &R, Pol< T > &I)const{
		return ReIm(ix,R,I,_order);
	}*/
};

template < typename T>
std::ostream& operator<<( std::ostream& os, Pol<T> const &P )
{
	os << "(" << " " << P[0];
	for(uint i=1; i<= P.order(); ++i) os << " + " << P[i] << "*x^" << i;
	os << ")";
	return os;
}

/*//---------------   skalare addition


template < typename T >
inline Pol< T > const operator +(T const & t, Pol< T > P){
	P+=t;
	return P;
}*/


/*//---------------   skalare multiplikation

template < typename T, typename T2  >
inline Pol< T > const operator *(T2 const & t, Pol< T > P){
	P*=t;
	return P;
}*/





//--------------------------- complex

/*template < typename T >
inline void operator =(Pol< T > &lhs,Pol< T > const &rhs){
	uint n= rhs.order();
	lhs.order(n);
	for(uint i=0; i<=n; ++i) lhs[i]= rhs[i];
}*/

template < typename T >
inline Pol< T > const real(Pol< T > const &P){
	uint n = P.order();
	Pol< T > Real;
	Real.order(n);
	for(uint i=0; i<=n; ++i)Real[i]=real(P[i]);
	return Real;
}

template < typename T >
inline Pol< T > const imag(Pol< T > const &P){
	uint n = P.order();
	Pol< T > Imag(n);
	for(uint i=0; i<=n; ++i)Imag[i]=imag(P[i]);
	return Imag;
}

template < typename T >
inline Pol< T > const exp(Pol< T > P){
	uint n=P.order();
	T P0 = P[0];//FIXME implicit ??
	P[0]-=P0;
	Pol< T > Exp(P);
	//Exp[0]-=P0;
	Pol< T > tmp(P);
	for(uint i=1; i<n; ++i){
		tmp.mul(P,n);
		Exp+=tmp*Pol< T >(1./(fac(i+1)));//FIXME implicit ??
	}
	P0=exp(P0);
	for(uint i=1; i<=n; ++i)Exp[i]*=P0;
	//Exp*=P0;
	Exp[0]+=P0;
	return Exp;
}

/*template < typename T >
inline Pol< complex< T > > const PolvonRe(complex< T > X,Pol< complex< T > > const &P){
	Pol< complex< T > > tmp;
	complex< T > null=(X-X);
	uint N=P.order();
	tmp.order(N);

	tmp[N]=complex< T >(P[N].real(),null);

	for(uint n=(N-1); n>=0; --n){
		tmp[n]=complex< T >(P[n].real(),null);


		T reX = X.real();
		T imX = X.imag();
		for(uint l=1; (n+2*l)<=N; ++l){
			T b = bin(n+2*l,n);
			tmp[n]
		}
	}


	return tmp;
}*/

#endif
