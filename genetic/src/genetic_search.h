#define PARDIM 200	// size of all parameters array (par1 - Ag, ... , par100 - alpha, ...)
#define NDIM 85		// search parameters for optimization in the first NDIM (Ag, ..., Xx only)

#define NRANDPHASE 1000	// number of iterations for the random phase
#define NMINPHASE 1000	// number of iterations for minimization cooling down
#define POPULATION 80	// size of the population 
#define NBEST 10	// number of best minima survive after the random phase 
#define SURVIVE 40	// number of elemenets go to the next round

#define INVT 3.		// cooling parameter


// access to initial conditions (read already by MINUIT)
extern struct init_name {

  char CPNAM[PARDIM][10];

} mn7nam_;

extern struct init_coord {

  double U[PARDIM];
  double ALIM[PARDIM];
  double BLIM[PARDIM];

} mn7ext_;

extern struct npar {

  int maxint;
  int NPAR;
  int maxext;
  int nu;

} mn7npr_;

extern struct debug {

  bool lprint;

} doprint_;

// description of a single point
typedef struct Creature {

  double* position;
  int age;
  double chi2;

} Creature;

// the zoo operations
inline void swap(int i, int j);
inline int  partition(int left, int right, int pivot);
inline void qsort(int left, int right);

// pointer to the chi2 calculator
typedef double (*Fcn)(int*, double*, double*, double*, int*, int*);
Fcn fcn;

// c++ functions
void   genetic_search(Fcn fcn_in);
void   ini_zoo();
void   eliminate();
void   generate_new_creature(int position);
void   generate_offspring(int parent, int child, int generation);
void   print_candidate_parameters(int ipoint, int icand);
double call_fcn(double* pars, int iflag = 4);
void   run_std_fcn();
void   copy_file(const char* filename, const char* dirname);

// wrapper to make the main search function callable from FORTRAN
// the construction here seems to be fragile
// fix it if possible
extern "C" {
  void genetic_search_(Fcn fcn_in) {
    return genetic_search(fcn_in);
  };
}

// make singleton to keep all the global info
// let it be class to allow further development possible
class Zoo {

  public:

    static Zoo* Instance();

    static Creature* creatures[POPULATION];
    static Creature* best_random[NBEST];
    // number of non-zero members in mn7npr_.U[]
    static int NOptPars;
    // array of pointers to non-zero members in mn7npr_.U[]
    static int OptParPointers[NDIM];
    // up and down limits for the parameters variations
    static double down_bound[NDIM];
    static double up_bound[NDIM];
    static double delta_bound[NDIM];
    // names of the parameters
    static char par_names[NDIM][13];

  private:

    // constructors are here to prohibit calls from outside
    Zoo() {};
    Zoo(Zoo const&) {};
    Zoo& operator=(Zoo const&) {return *this;};

    static Zoo* pInstance;

};

Creature* Zoo::creatures[POPULATION];
Creature* Zoo::best_random[NBEST];
int Zoo::NOptPars;
int Zoo::OptParPointers[NDIM];
double Zoo::down_bound[NDIM];
double Zoo::up_bound[NDIM];
double Zoo::delta_bound[NDIM];
char Zoo::par_names[NDIM][13];
Zoo* Zoo::pInstance = NULL;
 
