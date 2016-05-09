#include "lhafunc.h"

Matrix< 2, complex<double> > Pauli(int i){
	Matrix<2,complex<double> > M(0.);
	switch (i) {
	case 1:	M[0][1]=1.;
			M[1][0]=1.;
			break;
	case 2:	M[0][1]=complex<double>(0,-1.);
			M[1][0]=complex<double>(0,+1.);
			break;
	case 3:	M[0][0]=+1.;
			M[1][1]=-1.;
			break;
	default:
		break;
	}
	return M;
}

std::ostream& operator<<( std::ostream& os, VectorPhi const &V )
{
	os << "(";
	for(e_parton i=0; i<e_g; ++i) os << "\t" << V[i] ;
	os << " )" << endl
	   << "(\t\t\t" << V[e_g]  << "\t\t\t)" << endl
	   << "(";
	for(e_parton i=e_g+1; i<N_PARTONS; ++i) os << "\t" << V[i] ;
	os << " )";
	return os;
}

