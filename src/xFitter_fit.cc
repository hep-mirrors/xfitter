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
}


// Some c connection:
void fit() {
  do_fit_();
  close_theor_eval_();
  int io = 6;
  hf_errsum_(io);
}


// Test to move an array around

float sumvals(int array_length, boost::python::numeric::array in) {
  float value = 0.0;
  for (unsigned int i = 0; i<array_length; i++) {
    value += boost::python::extract<float>(in[i]);
  }
  return value;
}


// Python interface
BOOST_PYTHON_MODULE(libxfitter_fit)
{
    using namespace boost::python;

    boost::python::numeric::array::set_module_and_type("numpy","ndarray");
    def("sumvals",sumvals);

    def("logo", hfbanner_);
    def("read_steer",read_steer_);
    def("read_data" ,read_data_);
    def("init_theory",init_theory_modules_);
    def("fit",fit);
}
