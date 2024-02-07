/*
   @file ReactionRT_DISNC.cc
   @date 2017-04-10
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-10
*/

#include <iostream>
#include "ReactionRT_DISNC.h"
#include "xfitter_cpp_base.h"
#include "ext_pdfs.h"

#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include "xfitter_steer.h"

// the class factories
extern "C" ReactionRT_DISNC *create()
{
  return new ReactionRT_DISNC();
}

struct stf
{
  int idx;
  double f2;
  double fl;
  double f2c;
  double flc;
  double f2b;
  double flb;
};

// RT wrappers from RT/src/mstw2008_wrap.f:
extern "C"
{
  void mstwnc_wrap_(const double &x, const double &q2, const int &ipn,
                    double &f2, double &f2c, double &f2b, double &fl, double &flc, double &flb,
                    const int &iflag, const int &index, const double &f2QCDNUM, const double &flQCDNUM,
                    const int &usekfactors = 0);
  void rt_setalphas_(const double &alphaSzero);
  void rt_set_input_(const double *varin, const double &mCharmin, const double &mBottomin, const double &alphaSQ0in,
                     const double &alphaSMZin, const int &alphaSorderin, const double &alphaSnfmaxin, const int &iordin);
  void wate96_();
}

// Initialize at the start of the computation
void ReactionRT_DISNC::atStart()
{
  Super::atStart();
}

void ReactionRT_DISNC::initTerm(TermData *td)
{
  Super::initTerm(td);

  // Allocate internal arrays:
  unsigned termID = td->id;
  _f2rt[termID].resize(GetNpoint(termID));
  _flrt[termID].resize(GetNpoint(termID));
}

//
void ReactionRT_DISNC::atIteration()
{
  Super::atIteration();
  // Flag for internal arrays
  for (auto ds : _dsIDs)
  {
    (_f2rt[ds])[0] = -100.;
    (_flrt[ds])[0] = -100.;
  }
}

//
void ReactionRT_DISNC::compute(TermData *td, valarray<double> &val, map<string, valarray<double>> &err)
{
  // First init, then call base class:
  td->actualizeWrappers();
  vector<double> varin = {*td->getParamD("varin0"), *td->getParamD("varin1"), *td->getParamD("varin2"), *td->getParamD("varin3")}; // {0.0, 1.0, -2./3., 1.0};
  const double mc = *td->getParamD("mch");
  const double mb = *td->getParamD("mbt");
  const double mZ = *td->getParamD("Mz");
  const double qs0 = 1.0;
  const double as_q0 = alphas_wrapper_(sqrt(qs0));
  const double as_MZ = alphas_wrapper_(mZ);

  const string order = td->getParamS("Order");

  const int iord = OrderMap(order) - 1;
  const int asOrederIn = 0; // ???
  const int alphaSnfmaxin = 3;

  // set PDFs, alphaS functions:
  rt_set_pdfs_alphaS(pdf_xfxq_wrapper_, alphas_wrapper_);

  rt_set_input_(&varin[0], mc, mb, as_q0, as_MZ, asOrederIn, alphaSnfmaxin, iord);
  wate96_();

  Super::compute(td, val, err);
}

//
void ReactionRT_DISNC::F2 BASE_PARS
{
  unsigned termID = td->id;

  valarray<double> f2base, f2gamma_base;
  valarray<double> f2gamma_RT(GetNpoint(termID));

  // Get RT F2gamma
  F2gamma_RT(td, f2gamma_RT, err);

  if (GetDataFlav(termID) == dataFlav::incl)
  {
    // Get ZMVFNs F2s:
    Super::F2gamma(td, f2gamma_base, err);
    Super::F2(td, f2base, err);

    // Re-scale F2:
    if (td->getParamI("additive") == 1)
      val = f2gamma_RT + f2base - f2gamma_base;
    else
      val = f2base * f2gamma_RT / f2gamma_base;
  }
  else
    val = f2gamma_RT;
}

void ReactionRT_DISNC::FL BASE_PARS
{
  unsigned termID = td->id;
  valarray<double> flbase, flgamma_base;
  valarray<double> flgamma_RT(GetNpoint(termID));

  // Get RT F2gamma
  FLgamma_RT(td, flgamma_RT, err);

  // OZ 19.10.2017 TODO: in dis_sigma.f there is no rescaling for FL at order = 1, should it be here?
  if (GetDataFlav(termID) == dataFlav::incl)
  {
    // Get ZMVFNs F2s:
    Super::FLgamma(td, flgamma_base, err);
    Super::FL(td, flbase, err);

    // Re-scale FL:
    if (td->getParamI("additive") == 1)
      val = flgamma_RT + flbase - flgamma_base;
    else
      val = flbase * flgamma_RT / flgamma_base;
  }
  else
    val = flgamma_RT;
}

void ReactionRT_DISNC::F2gamma_RT BASE_PARS
{
  auto termID = td->id;
  calcF2FL(td);
  val = _f2rt[termID];
}

void ReactionRT_DISNC::FLgamma_RT BASE_PARS
{
  auto termID = td->id;
  calcF2FL(td);
  val = _flrt[termID];
}

// Place calculations in one function, to optimize calls.
void ReactionRT_DISNC::calcF2FL(TermData *td)
{
  unsigned termID = td->id;
  if (!(_f2rt[termID][0] < -99.))
    return;
  
  // compute
  // Get x,Q2 arrays:
  auto *q2p = GetBinValues(td, "Q2"), *xp = GetBinValues(td, "x");
  auto q2 = *q2p, x = *xp;

  const size_t Np = GetNpoint(termID);
  int iflag = 1;

  double f2(0), f2b(0), f2c(0), fl(0), flc(0), flb(0);

  int threads = xfitter::xf_ncpu( td->getParamI("threads") );
  std::cout << "threads " << threads << "\n";
  if (threads < 2)
    {
      for (size_t i = 0; i < Np; i++)
	{
	  if (q2[i] > 1.0)
	    {

	      mstwnc_wrap_(x[i], q2[i], 1,
			   f2, f2c, f2b, fl, flc, flb,
			   iflag, i + 1, 1., 0.1, 0);
	    }

	  switch (GetDataFlav(termID))
	    {
	    case dataFlav::incl:
	      _f2rt[termID][i] = f2;
	      _flrt[termID][i] = fl;
	      break;
	    case dataFlav::c:
	      _f2rt[termID][i] = f2c;
	      _flrt[termID][i] = flc;
	      break;
	    case dataFlav::b:
	      _f2rt[termID][i] = f2b;
	      _flrt[termID][i] = flb;
	      break;
	    }
	}
      return;
    }
      
  //fork wait parallelisation
  int fd[2];
  if (pipe(fd) < 0)
    {
      std::cout << "Error in fork/wait: could not create pipe" << std::endl;
      exit(-1);
    }

  size_t Npr = Np/threads+1;
  //std::cout << " Np " << Np << " Npr " << Npr << std::endl;
  for (int P = 0; P < threads; P++)
    {
      pid_t id = xfitter::xf_fork(threads);
      if (id == 0)
	{
	  close(fd[0]);
	  for (size_t i = P*Npr; i < std::min(Np,(P+1)*Npr); i++)
	    {
	      if (!(q2[i] > 1.0))
		continue;
	      stf fs;

	      //std::cout << P << "  " << i << std::endl;
	      mstwnc_wrap_(x[i], q2[i], 1,
			   f2, f2c, f2b, fl, flc, flb,
			   iflag, i + 1, 1., 0.1, 0);
	      fs.f2 = f2;
	      fs.fl = fl;
	      fs.f2c = f2c;
	      fs.flc = flc;
	      fs.f2b = f2b;
	      fs.flb = flb;
	      fs.idx = i;
	  
	      int status;
	      status = write(fd[1], &fs, sizeof fs);
	    }
	  exit(0);
	}
      else if (id < 0)
	{
	  std::cout << "Error: failed to fork" << std::endl;
	  exit (-1);
	}
      
    }
  //wait for all children to finish
  int status;
  pid_t wpid;
  while ((wpid = wait(&status)) > 0)
    if (status < 0)
      {
	std::cout << "Process " << wpid << " terminated with status " << status << std::endl;
	exit(-1);
      }
  
  //Read out buffer
  close(fd[1]);
  for (size_t i = 0; i < Np; i++)
    {
      if (!(q2[i] > 1.0))
	continue;
      stf fs;
      int nbytes = read(fd[0], &fs, sizeof fs);
      if(!(nbytes > 0))
	{
	  string message = "E: Error in fork/wait: nothing on the pipe.";
	  hf_errlog_(22082501, message.c_str(), message.size());
	}
	
      switch (GetDataFlav(termID))
	{
	case dataFlav::incl:
	  _f2rt[termID][fs.idx] = fs.f2;
	  _flrt[termID][fs.idx] = fs.fl;
	  break;
	case dataFlav::c:
	  _f2rt[termID][fs.idx] = fs.f2c;
	  _flrt[termID][fs.idx] = fs.flc;
	  break;
	case dataFlav::b:
	  _f2rt[termID][fs.idx] = fs.f2b;
	  _flrt[termID][fs.idx] = fs.flb;
	  break;
	}
    }
  close(fd[0]);
}
