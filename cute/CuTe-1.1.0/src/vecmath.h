#ifndef VECMATH_H
#define VECMATH_H

/*#ifndef uint
#define uint unsigned int
#endif*/

#include <vector>
#include <ostream>
typedef unsigned int  uint;

template <uint N, typename T>
class Matrix;

template <uint N, typename T>
class Vector{
protected:
	T *_V;


	inline T const min(T const &a, T const &b)const{return ( (a<b) ? a : b );}

	inline T const max(T const &a, T const &b)const{return ( (a>b) ? a : b );}

	inline T const power(T const &x, uint const n)const{
		if(n==1) return x;
		if(n%2) return x*power(x,n-1);
		T tmp(power(x,n/2));
		return tmp*tmp;
	}

public:
//---------------Konstruktoren
	inline Vector<N, T >()
		:_V(new T[N]){
	}

	inline Vector<N, T >(T const &t)
		:_V(new T[N]){
		for(uint i=0; i<N; ++i)_V[i]=t;
	}

	inline Vector< N,T >(Vector<N,T> const& V)
		:_V(new T[N]){
		for(uint i=0; i<N; ++i)_V[i]=V[i];
	}

	inline Vector< N,T >(T (&array)[N])
		:_V(new T[N]){
		for(uint i=0; i<N; ++i)_V[i]=array[i];
	}

	inline Vector< N,T >(std::vector<T> const& v)
		:_V(new T[N]){
		for(uint i=0; i<N; ++i)_V[i]=v[i];
	}

	~Vector< N,T >(){
		delete[] _V;
	}
//----------------	get
	inline int const n(){return N;}
	inline T &operator [] (int const i){ return _V[i];}
	inline T const &operator [] (int const i)const {return _V[i];	}
//----------------	set
	inline Vector< N,T >& operator =(Vector<N, T > const& V){
		for(uint i=0; i<N; ++i) _V[i]=V[i];
		return *this;
	}
//----------------	funcs
	inline T const p_norm(uint const p)const{
		T p_norm = power(_V[0],p);
		for(uint i=1; i<N; ++i) p_norm+=power(_V[0],p);
		return p_norm;
	}
//-----------------  standard calc
	template< typename T2 >
	inline Vector< N,T >& operator +=(Vector< N,T2 > const& V){
		for(uint i=0; i<N; ++i) _V[i]+=V[i];
		return *this;
	}

	template< typename T2 >
	inline Vector< N,T >& operator -=(Vector< N,T2 > const& V){
		for(uint i=0; i<N; ++i) _V[i]-=V[i];
		return *this;
	}

	template< typename T2 >
	inline Vector< N,T > const operator +(Vector< N,T2 > const& V)const{
		Vector< N,T > tmp;
		for(uint i=0; i<N; ++i) tmp[i]=_V[i]+V[i];
		return tmp;
	}

	template< typename T2 >
	inline Vector< N,T > const operator -(Vector< N,T2 > const& V)const{
		Vector< N,T > tmp(*this);
		tmp-=V;
		return tmp;
	}

//-----------------  scalar multiplikation
	inline Vector< N,T >& operator*=(T const &t){
		//std::cout << "***************************\tV*S" << std::endl;
		for(uint i=0; i<N; ++i) _V[i]*=t;
		return *this;
	}

	inline Vector< N,T > const operator*(T const &t)const{
		Vector< N,T > tmp(*this);
		tmp*=t;
		return tmp;
	}

	/*inline Vector< N,T >& operator*=(double t){
		std::cout << "***************************\tV*S" << std::endl;
		for(uint i=0; i<N; ++i) _V[i]*=t;
		return *this;
	}

	inline Vector< N,T > const operator*(double t)const{
		std::cout << "***************************\tV*S" << std::endl;
		Vector< N,T > tmp(*this);
		tmp*=t;
		return tmp;
	}*/

	inline Vector< N,T >& operator/=(T const &t){
		//std::cout << "***************************\tV/S" << std::endl;
		for(uint i=0; i<N; ++i) _V[i]/=t;
		return *this;
	}

	inline Vector< N,T > const operator/(T const &t)const{
		Vector< N,T > tmp(*this);
		tmp/=t;
		return tmp;
	}

	/*inline Vector< N,T >& operator/=(double t){
		std::cout << "///////////////////////////\tV/S" << std::endl;
		for(uint i=0; i<N; ++i) _V[i]/=t;
		return *this;
	}

	inline Vector< N,T > const operator/(double t)const{
		std::cout << "///////////////////////////\tV/S" << std::endl;
		Vector< N,T > tmp(*this);
		tmp/=t;
		return tmp;
	}*/

//------------------------	skalarproduct

	template < typename T2 >
	inline T const operator *(Vector< N,T2 > const& V)const{
		//std::cout << "***************************\tV*V" << std::endl;
	//inline Vector< N,T >& operator *=(T const& t){
		T tmp = _V[0]*V[0];
		for(uint i=1; i<N; ++i) tmp+=_V[i]*V[i];
		return tmp;
	}

//------------------------	matrix calc

	/*template < typename T2 >
	inline Vector< N,T >& operator *=(Matrix< N,T2 > M);
	template < typename T2 >
	inline Vector< N,T > const operator *(Matrix< N,T2 > const &M)const;*/
};

template < uint N, typename T >
inline std::ostream& operator<<(std::ostream& os,  Vector< N,T > const & V)
{
	os << "( " << V[0] << ",";
    for(uint i=1; i<(N-1); ++i) os << "\t" << V[i] << ",";
    os << "\t" << V[N-1] << " )";
    return os;
}


/*template < uint N, typename T, typename T2>
inline Vector< N,T > const operator *(T2 const& t, Vector< N,T > V){
	V*=t;
	return V;
}*/
//////////////////////////////////////////////////////////
/////////////		Matrix	/////////////
////////////////////////////////////////////////////////////

template < uint N, typename T >
class Matrix{
protected:
	//Vector< N,Vector< N,T > > _M;
	Vector< N,T > *_M;
	//T **_M;

public:
//---------------Konstruktoren
	inline Matrix< N,T >()
		:_M(new Vector< N,T>[N]){
		//:_M(new T[N][N]){
	}

	inline Matrix< N,T >(T const &t)
		:_M(new Vector< N,T>[N]){
		//:_M(new T[N][N]){
		for(uint i=0; i<N; ++i)for(uint j=0; j<N; ++j) _M[i][j]=t;
	}

	inline Matrix< N,T >(Matrix< N,T > const& M)
		:_M(new Vector< N,T>[N]){
		//:_M(new T[N][N]){
		for(uint i=0; i<N; ++i)for(uint j=0; j<N; ++j) _M[i][j]=M[i][j];
	}

	inline Matrix< N,T >(T (&array)[N][N])
		:_M(new Vector< N,T>[N]){
		//:_M(new T[N][N]){
		for(uint i=0; i<N; ++i)for(uint j=0; j<N; ++j) _M[i][j]=array[i][j];
	}

	~Matrix< N,T >(){
		delete[] _M;
	}


//---------------	get
	inline int const n()const{return N;}
	inline Vector< N,T >& operator [] (uint const i){ return _M[i];}
	inline Vector< N,T > const&  operator [] (uint const i)const {return _M[i];	}

//---------------	set
	inline Matrix< N,T >& operator =(Matrix< N,T > const& M){
		for(uint i=0; i<N; ++i)for(uint j=0; j<N; ++j) _M[i][j]=M[i][j];
		return *this;
	}

	inline void identity(T const &eins = 1., T const &null = 0.){
		for(uint i=0; i<N;++i)for(uint j=0; j<N;++j)
			if(i==j)_M[i][j]=eins;
			else _M[i][j]=null;
	}
//---------------	funcs

	inline Matrix< N,T >& transpose(){
		//T tmp;
		for(uint i=0; i<N; ++i) for(uint j=i+1; j<N; j++){
			std::swap(_M[i][j],_M[j][i]);
			/*tmp =_M[i][j];
			_M[i][j]=_M[j][i];
			_M[j][i]=tmp;*/
		}
		return *this;
	}

//---------------	standard calc
	template< typename T2 >
	inline Matrix< N,T >& operator +=(Matrix< N,T2 > const& M){
		for(uint i=0; i<N; ++i) _M[i]+=M[i];
		return *this;
	}
	template< typename T2 >
	inline Matrix< N,T >const operator +(Matrix< N,T2 > M)const{
		M+=(*this);
		return M;
	}

	template< typename T2 >
	inline Matrix< N,T >& operator -=(Matrix< N,T2 > const& M){
		for(uint i=0; i<N; ++i) _M[i]-=M[i];
		return *this;
	}
	template< typename T2 >
	inline Matrix< N,T >const operator -(Matrix< N,T2 > M)const{
		M-=(*this);
		return M;
	}

	template< typename T2 >
	inline Matrix< N,T >& operator *=(Matrix< N,T2 > M){
		M.transpose();
		Matrix< N,T > tmp(*this);
		for(uint i=0; i<N; ++i) for(uint j=0; j<N; j++) {
			_M[i][j]=tmp[i]*M[j];
		}
		return *this;
	}
	inline Matrix< N,T >const operator *(Matrix< N,T > M)const{
		M*=(*this);
		return M;
	}

	//template < typename T2 >
	inline Matrix< N,T >& operator *=(T const& t){
		for(uint i=0; i<N; ++i) _M[i]*=t;
		return *this;
	}
	//template < typename T2 >
	inline Matrix< N,T >const operator *(T const& t)const{
		Matrix< N,T > tmp(*this);
		tmp*=t;
		return tmp;
	}

	//template < typename T2 >
	inline Matrix< N,T >& operator /=(T const& t){
		for(uint i=0; i<N; ++i) _M[i]/=t;
		return *this;
	}
	//template < typename T2 >
	inline Matrix< N,T >const operator /(T const& t)const{
		Matrix< N,T > tmp(*this);
		tmp/=t;
		return tmp;
	}

};

template < uint N, typename T >
inline std::ostream& operator<<(std::ostream& os,  Matrix< N,T > const & M)
{
	for(uint i=0; i<(N-1); ++i){
		os << M[i] << std::endl;
	}
	os << M[N-1];

    return os;
}


template < uint N, typename TV, typename TM >
inline Vector< N,TV > const operator *(Matrix< N,TM > const &M, Vector< N,TV > const &V){
	//std::cout << "***************************\tM*V" << std::endl;
	Vector< N,TV > tmp;
	for(uint i=0; i<N; ++i){
		//std::cout << "tmp["<<i<<"] = " <<V[0]<<"*"<<M[i][0];
		tmp[i]=V[0]*M[i][0];
		for(uint j=1; j<N; ++j){
			//std::cout << " + " <<V[j]<<"*"<<M[i][j];
			tmp[i]+=V[j]*M[i][j];
		}
		//std::cout <<  " = " << tmp[i] << endl;
	}
	return tmp;
}

template < uint N, typename TV, typename TM >
inline Vector< N,TV > const operator *(Vector< N,TV > const &V, Matrix< N,TM > M){
	//std::cout << "***************************\tV*M" << std::endl;
	M.transpose();
	return (M*V);
}


#endif
