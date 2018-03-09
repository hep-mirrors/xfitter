#include <boost/python.hpp>
#include <iostream>
#include <list>
#include <typeinfo>

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

// print theory

double theo(int ibin) {
  return c_theo_.theo[ibin];
}

int ndata() {
  return cndatapoints_.npoints;
}

// some test object


struct Base
{
    virtual ~Base() {}
    virtual int f() { return 0; }
};

struct BaseWrap : Base, boost::python::wrapper<Base>
{
    int f()
    {
        if (boost::python::override f = this->get_override("f"))
            return f(); // *note*
        return Base::f();
    }

    int default_f() { return this->Base::f(); }
};

int calls_f(Base& b) { return b.f(); }

////////////////////////////////////////////////////////////////////////////////

struct Action
{
  Action(){}
  virtual ~Action() {}
  virtual void Initialize() { }
  virtual void AtIteration() {std::cout << "At Iter " << std::endl;}
};

struct ActionWrap : Action, boost::python::wrapper<Action>
{
    void Initialize()
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
  
    void AtIteration()
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

    void default_Initialize() { return this->Action::Initialize(); }
    void default_AtIteration() { return this->Action::AtIteration(); }
};

void calls_action_iteration(Action& a) {
  a.AtIteration();
  std::cout << typeid(a).name() << std::endl;
}

// a list to store actions
std::list<Action> _actionList_;


void actionIterate() {
  for ( auto  &action : _actionList_) {
    action.AtIteration();
    std::cout << typeid(action).name() << std::endl;
  }
}

void addAction(Action& a) {
  std::cout << typeid(a).name() << std::endl;
  _actionList_.push_back(a);
  actionIterate();
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
    def("theo",theo);
    def("ndata",ndata);


    class_<BaseWrap, boost::noncopyable>("Base")
      .def("f", &Base::f, &BaseWrap::default_f)
      ;
    def("calls_f",calls_f);


    class_<ActionWrap, boost::noncopyable>("Action")
      .def("Initialize", &Action::Initialize, &ActionWrap::default_Initialize)
      .def("AtIteration", &Action::AtIteration, &ActionWrap::default_AtIteration)
      ;
    def("addAction",addAction);
    def("actionIteration",calls_action_iteration);
    def("actionIterate",actionIterate);
}
