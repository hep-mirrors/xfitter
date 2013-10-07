#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "genetic_search.h"
#include "mixmax_wrapper.h"



//------------------------------------------------------------------------------
void genetic_search(Fcn fcn_in)
{
  // the main function

  setbuf(stdout, NULL); // don't buffer the output
  printf("\n");

  // associate the chi2 calculator with the global name
  fcn = fcn_in;

  // run first time with std parameters from the steering
  // to avoid problems with pre-calculation (see RT FAST scheme)
  run_std_fcn();

  // initialize a random numbers generator
  ini_rnd();

  // initialize the global arrays of the points
  (void)Zoo::Instance();
  ini_zoo();

  // first make random guesses
  // then walk around good tries

  // do loop over generations
  for (int igen = 0; igen < NRANDPHASE; ++igen) {
    printf("GENETIC: generation %i (random phase)\n", igen);

    for (int icr = NBEST; icr < POPULATION; ++icr) {
      generate_new_creature(icr);
    }
    eliminate();	// sort creatures to get the best on the top of the array

  }


  // copy the best creatures after the first phase to the separate zoo
  for (int icr = 0; icr < NBEST; ++icr) {
    Creature* creature = (Creature*)malloc(sizeof(Creature));
    memcpy(creature, Zoo::creatures[icr], sizeof(Creature));
    Zoo::best_random[icr] = creature;
  }

  // print the results of the RandomPhase
  printf("GENETIC: The results (random phase)\n");
  for (int i = 0; i < NBEST; ++i) {
    printf("%i. age = %i   chi2 = %f\n", i, Zoo::best_random[i]->age, Zoo::best_random[i]->chi2);
  }
  printf("GENETIC: start an evolution phase\n");

  // minimize chi2 for each of the best points from
  // the previous step
  for(int ibcr = 0; ibcr < NBEST; ++ibcr) {
    // fill the zoo with large chi2 solutions
    // not to stuck with the best point for
    // all the other good points
    ini_zoo();
    // flood the first (SURVIVE) positions with
    // the given good creature from RP
    for (int icr = 0; icr < SURVIVE; ++icr) {
      memcpy(Zoo::creatures[icr], Zoo::best_random[ibcr], sizeof(Creature));
    }

    // do minimization search
    // (actually any search algorithm could go here)
    for (int igen = 0; igen < NMINPHASE; ++igen) {
      printf("GENETIC: generation %i (evolution phase)\n", igen);

      for (int icr = 0; icr < SURVIVE; ++icr) {
	int parent = icr;
	int child  = icr + SURVIVE;
	generate_offspring(parent, child, igen);
      }
      eliminate();

      // print the current results of the evolution phase 
      if (igen % 500 == 0) {
	printf("GENETIC: The results (evolution phase)\n");
	for (int i = 0; i < NBEST; ++i) {
	  printf("%i. age = %i   chi2 = %f\n", i, Zoo::best_random[i]->age, Zoo::creatures[i]->chi2);
	}
      }

    }

    // print the results of the mimimization for the random phase candidate
    print_candidate_parameters(0, ibcr);

  }

}


//------------------------------------------------------------------------------
void print_candidate_parameters(int ipoint, int icand)
{
  // print the parameters obtained by the genetic fit
  // in MINUIT-like form

  // open file
  char scand[10];
  sprintf(scand, "%d", icand);
  char dirname[100];
  strcpy(dirname, "output/genetic.");
  strcat(dirname, scand);
  struct stat st = {0};
  if (stat(dirname, &st) == -1) {
    mkdir(dirname, 0777);
  }

  char filename[100];
  strcpy(filename, dirname);
  strcat(filename, "/");
  strcat(filename, "minuit_genetic.out.");
  strcat(filename, scand);
  strcat(filename, ".txt");    
  printf("filename %s", filename);
  FILE* fout = fopen(filename, "w");

  // dump the parameters
  fprintf(fout, "C   The fit parameters for the candidate n%i of age %i with chi2 = %f\n\n", icand, Zoo::creatures[ipoint]->age, Zoo::creatures[ipoint]->chi2);
  for (int ipos = 0; ipos < NDIM; ++ipos) {
    if (Zoo::par_names[ipos][1] != ')') {
      char thename[20];
      for (int j = 0; j < 20; ++j) thename[j] = ' ';
      thename[0] = '\0';
      thename[19] = '\0';
      strcpy(thename, Zoo::par_names[ipos]);
      thename[strlen(Zoo::par_names[ipos])] = ' ';
      float step  = (mn7ext_.BLIM[ipos] - mn7ext_.ALIM[ipos]) / 100.; // set here some step for a further minimization
      fprintf(fout, "  %i    %s%f   %f   %f   %f\n", ipos+1, thename, Zoo::creatures[ipoint]->position[ipos], step, mn7ext_.ALIM[ipos], mn7ext_.BLIM[ipos]);
    }
  }

  
  // call the same fcn() with the special flag to prepare the text files
  call_fcn(Zoo::creatures[ipoint]->position, 3);
  copy_file("pdfs_q2val_01.txt", dirname);
  copy_file("pdfs_q2val_02.txt", dirname);
  copy_file("pdfs_q2val_03.txt", dirname);
  copy_file("pdfs_q2val_04.txt", dirname);
  copy_file("pdfs_q2val_05.txt", dirname);
  copy_file("pdfs_q2val_06.txt", dirname);
  copy_file("lhapdf.block.txt",  dirname);
  copy_file("fittedresults.txt", dirname);
  copy_file("parsout_0",	 dirname);
  copy_file("Results.txt",	 dirname);


  // close file
  fclose(fout);
  printf("GENETIC: Best fit parameters for the candidate %i were dumped to %s\n", icand, filename);

}


//------------------------------------------------------------------------------
void copy_file(const char* filename, const char* dirname)
{
  char infilename[100];
  strcpy(infilename, "output/");
  strcat(infilename, filename);

  char outfilename[100];
  strcpy(outfilename, dirname);
  strcat(outfilename, "/");
  strcat(outfilename, filename);

  FILE* infile  = fopen(infilename,  "r");
  FILE* outfile = fopen(outfilename, "w");

  if (infile != NULL) {
    char ch;
    while( ( ch = fgetc(infile) ) != EOF ) fputc(ch, outfile);
    fclose(infile);
  }

  fclose(outfile);

}

//------------------------------------------------------------------------------
void ini_zoo()
{
  // fill the initial objects

  for (int icr = 0; icr < POPULATION; ++icr) {
    Creature* creature = (Creature*)malloc(sizeof(Creature));
    double* pos = (double*)malloc(PARDIM*sizeof(double));
    creature->position = pos;
    // copy all the coordinates like non-optimized pars, alpha-s, ...
    for (int ipos = 0; ipos < PARDIM; ++ipos) creature->position[ipos] = mn7ext_.U[ipos];
    Zoo::creatures[icr] = creature;
    // randomize creature not to start from the local minimum
    generate_new_creature(icr);	
  }

}


//------------------------------------------------------------------------------
void eliminate()
{
  // sort the creatures according to their chi2,
  // kill some of them, etc.

  // sort over chi2; use qsort
  qsort(0, POPULATION-1);

}


//------------------------------------------------------------------------------
// ALL IS HERE
//------------------------------------------------------------------------------
void generate_new_creature(int position)
{
  // generate a child from the parent

  // the function _does_not_touch_ parameters should not be optimized, so do it manually if needed

  const int Npar = Zoo::NOptPars;
  for (int ipos = 0; ipos < Npar; ++ipos) {
    double coord = get_uniform_double(Zoo::down_bound[ipos], Zoo::up_bound[ipos]);  
    Zoo::creatures[position]->position[Zoo::OptParPointers[ipos]] = coord;
  }
  Zoo::creatures[position]->chi2 = call_fcn(Zoo::creatures[position]->position);
  Zoo::creatures[position]->age = 0;

}


//------------------------------------------------------------------------------
void generate_offspring(int parent, int child, int generation)
{
  // generate a child from the parent

  Zoo::creatures[parent]->age += 1;

  const double gmax = NMINPHASE * 100./99.;

  const int Npar = Zoo::NOptPars;
  for (int ipos = 0; ipos < Npar; ++ipos) {
    int pos_pointer = Zoo::OptParPointers[ipos];
    double parent_pos = Zoo::creatures[parent]->position[pos_pointer];
    double sigma = Zoo::delta_bound[ipos]/2./INVT * (1. - generation / gmax);
    double delta = get_gaus_double(0., sigma);
    Zoo::creatures[child]->position[pos_pointer] = parent_pos + delta;
  }
  Zoo::creatures[child]->chi2 = call_fcn(Zoo::creatures[child]->position);
  Zoo::creatures[child]->age = 0;

}
//------------------------------------------------------------------------------
// END OF BLOCK
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
Zoo* Zoo::Instance()
{
  // initialization of the global variables
  if (Zoo::pInstance != NULL) {

    return Zoo::pInstance;

  } else {
    
    // suppress the output at each fcn calculation
    doprint_.lprint = false;

    // search for the non-zero parameters and
    // fill the arrays of limits, pointers to
    // non-zero parameters, etc.
    NOptPars = 0;
    for (int ipos = 0; ipos < NDIM; ++ipos) {
      // if the parameter is under optimization, fill limits, etc. 
      if (mn7ext_.BLIM[ipos] - mn7ext_.ALIM[ipos] > 0.) {
	Zoo::OptParPointers[NOptPars] = ipos;
	Zoo::down_bound[NOptPars] = mn7ext_.ALIM[ipos];
	Zoo::up_bound[NOptPars] = mn7ext_.BLIM[ipos];
	Zoo::delta_bound[NOptPars] = Zoo::up_bound[NOptPars] - Zoo::down_bound[NOptPars];
	++NOptPars;
      }
      // copy param name
      Zoo::par_names[ipos][0] = '\'';
      int ic = 1;
      for (; ic < 11; ++ic) {
	char c = mn7nam_.CPNAM[ipos][ic-1];
	if (c >= '!') {
	  Zoo::par_names[ipos][ic] = c;
	} else {
	  break;
	}
      }
      Zoo::par_names[ipos][ic] = '\'';
      Zoo::par_names[ipos][ic+1] = '\0';
    }
    Zoo::NOptPars = NOptPars;

    pInstance = new Zoo();

    printf("GENETIC: optimization goes over %i parameters\n", NOptPars);

    return pInstance;
  }
}

//------------------------------------------------------------------------------
double call_fcn(double* pars, int iflag)
{
  // call fcn function

  static int     npar    = NDIM;
  static double* g_dummy = (double*)calloc(NDIM, sizeof(double)); // fill array with zeroes with the calloc()
  static int     futil   = 0;					  // to feed int instead of function pointer is not a good idea indeed
  static double  chi2fcn;

  double chi2 = fcn(&npar, g_dummy, &chi2fcn, pars, &iflag, &futil);

  return chi2;
}


//------------------------------------------------------------------------------
void run_std_fcn()
{
  // run fcn for the first time to avoid negative "memory effects"
  printf("GENETIC: calling fcn with the initial parameters\n");
  double chi2 = call_fcn(mn7ext_.U);
  printf("GENETIC: the first fcn call on the initial parameters gives chi2 = %f\n", chi2);

}

//------------------------------------------------------------------------------
void swap(int i, int j)
{
  // swap two elements of the array

  if (i == j) return;

  Creature* tmp_cr = Zoo::creatures[i];
  Zoo::creatures[i] = Zoo::creatures[j];
  Zoo::creatures[j] = tmp_cr;
}


//------------------------------------------------------------------------------
void qsort(int left, int right)
{
  // qsort implementation

  int lengthm1 = right-left;
  if (lengthm1 < 1) {

    return;

  } else if (lengthm1 == 1) {  // just for a small acceleration

    if (Zoo::creatures[left]->chi2 > Zoo::creatures[right]->chi2) swap(left, right);
    return;

  } else {

    // choose a pivot value; use randomization?
    int pivot = (right + left)/2;
    int new_pivot = partition(left, right, pivot);

    qsort(left, new_pivot-1);
    qsort(new_pivot+1, right);

  }
}

int partition(int left, int right, int pivot)
{
  // function is used by qsort

  double pivotValue = Zoo::creatures[pivot]->chi2;
  swap(pivot, right);
  int new_pivot = left;

  for (int i = left; i < right; ++i) {
    if (Zoo::creatures[i]->chi2 < pivotValue) {
      swap(i, new_pivot);
      ++new_pivot;
    }
  }
  swap(new_pivot, right);

  return new_pivot;
}

