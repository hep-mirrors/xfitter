#include "param.h"

std::ostream& operator<<( std::ostream& os, CS_params const &p){
	p.print(os);
	return os;
}

std::ostream& operator<<( std::ostream& os,CS_params const *p ){
	p->print(os);
	return os;
}

