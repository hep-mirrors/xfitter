// pointer to the chi2 calculator
typedef double (*Fcn)(int*, double*, double*, double*, int*, int*);
Fcn fcn;

// c++ functions
void   genetic_search(Fcn fcn_in);

// wrapper to make the main search function callable from FORTRAN
// the construction here seems to be fragile
// fix it if possible
extern "C" {
  void genetic_search_(Fcn fcn_in) {
    return genetic_search(fcn_in);
  };
}
