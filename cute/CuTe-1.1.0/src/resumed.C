#include "resumed.h"

std::ostream& operator<<( std::ostream& os, const K_Integral* K ){
	K->print(os);
	return os;
}

std::ostream& operator<<( std::ostream& os, const Epol& E ){
	os << E[0] << endl << E[1];
	return os;
}

