#include "ReactionDummy.h"


extern "C" ReactionDummy* create() {
  return new ReactionDummy();
}

void ReactionDummy::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){

  unsigned termID = td->id;
  size_t N=val.size();
  for ( std::size_t i=0; i<N; i++) val[i]=1.0;

}
