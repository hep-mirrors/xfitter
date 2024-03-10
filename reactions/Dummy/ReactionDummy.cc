#include "ReactionDummy.h"


extern "C" ReactionDummy* create() {
  return new ReactionDummy();
}

void ReactionDummy::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){

  unsigned termID = td->id;
  val = 1.0;

}
