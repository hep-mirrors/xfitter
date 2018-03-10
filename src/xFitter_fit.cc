#include <boost/python.hpp>
#include <iostream>
#include <list>
#include <typeinfo>

#include "../include/action.h"
#include "../include/actions.h"
#include "../include/xfitter_cpp.h"

using boost::python::object;


// Fortran interface
extern "C" {
  void hfbanner_();
  void read_steer_();
  void read_data_();
  void init_theory_modules_();
  void do_fit_();
  void close_theor_eval_();
  void hf_errsum_(const int& io);
  void init_func_map_();
  int read_reactions_();
  void parse_params_();
  void init_rnd_seeds_();
  void init_func_map_();
}


// Some c connection:
void fit() {
  do_fit_();
  actions_finalize_();
  close_theor_eval_();
  int io = 6;
  hf_errsum_(io);
}

// some more connections

void init_pars() {
  init_rnd_seeds_();
  parse_params_();
  int i = read_reactions_();
}

void init_theo() {
  init_func_map_();
  init_theory_modules_();
}

// Full connection
void run() {
  hfbanner_();
  read_steer_();
  init_pars();
  read_data_();
  init_theory_modules_();
  fit();
}

///////////   Methods to access information

// print theory
double theo(int ibin) {
  return c_theo_.theo[ibin];
}

double getChi2() {
  return chi2fit_.chi2_tot;
}

int ndata() {
  return cndatapoints_.npoints;
}

// Python wrap for the actions class:

class ActionWrap : public Action, public boost::python::wrapper<Action>
{
public:
  virtual void Initialize() override final
  {
    if (boost::python::override Initialize = this->get_override("Initialize")) {
      Initialize();
      return; // *note*
    }
    else {
      Action::Initialize();
      return;
    }
  }
  
  virtual void AtIteration() override final
  {
    if (boost::python::override AtIteration = this->get_override("AtIteration")) {
      AtIteration();
      return; // *note*
    }
    else {
      Action::AtIteration();
      return;
    }
  }

  virtual void Finalize() override final
  {
    if (boost::python::override Finalize = this->get_override("Finalize")) {
      Finalize();
      return; // *note*
    }
    else {
      Action::Finalize();
      return;
    }
  }
  void default_Initialize() { return this->Action::Initialize(); }
  void default_AtIteration() { return this->Action::AtIteration(); }
  void default_Finalize() { return this->Action::Finalize(); }
};

///////////////////////////////////////////////////


void calls_action_iteration(Action& a) {
  a.AtIteration();
}



// Python interface
BOOST_PYTHON_MODULE(libxfitter_fit)
{
    using namespace boost::python;
    def("logo", hfbanner_);
    def("read_steer",read_steer_);
    def("read_data" ,read_data_);
    def("init_theory",init_theory_modules_);
    def("init_pars",init_pars);
    def("fit",fit);
    def("run",&run);    
    def("theo",theo);
    def("ndata",ndata);
    def("chi2",getChi2);
    
    // Actions
    class_<ActionWrap,boost::noncopyable>("Action")
      .def("Initialize", &Action::Initialize, &ActionWrap::default_Initialize)
      .def("AtIteration", &Action::AtIteration, &ActionWrap::default_AtIteration)
      .def("Finalize", &Action::Finalize, &ActionWrap::default_Finalize)
      ;
    def("addAction",&addAction);

    // For debugging
    def("actionIteration",calls_action_iteration);
    def("actionIterate",&actions_at_iteration_);

    //    register_ptr_to_python<ActionPtr>();
}
