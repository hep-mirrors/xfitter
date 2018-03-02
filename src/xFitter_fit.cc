#include <boost/python.hpp>
#include <iostream>
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
}
