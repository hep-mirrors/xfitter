#include "MSHT_Decomposition.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include <iostream>
#include <cmath>

namespace xfitter
{
  // Динамічне завантаження модуля
  extern "C" MSHT_Decomposition* create(const char* name) {
    return new MSHT_Decomposition(name);
  }

  MSHT_Decomposition::MSHT_Decomposition(const char* name) : BasePdfDecomposition{name} {}
  const char* MSHT_Decomposition::getClassName() const { return "MSHT_Decomposition"; }

  // Зчитування параметрів з parameters.yaml
  void MSHT_Decomposition::atStart() {
    using namespace std;
    const YAML::Node node = XFITTER_PARS::getDecompositionNode(_name);
    
    par_xuv      = getParameterisation(node["xuv"].as<string>());
    par_xdv      = getParameterisation(node["xdv"].as<string>());
    par_S        = getParameterisation(node["xS"].as<string>());
    par_splus    = getParameterisation(node["xsplus"].as<string>());
    par_sminus   = getParameterisation(node["xsminus"].as<string>());
    par_d_over_u = getParameterisation(node["xd_over_u"].as<string>());
    par_xg       = getParameterisation(node["xg"].as<string>());
  }

  // Правила сум (виконується на кожній ітерації фіту)
  void MSHT_Decomposition::atIteration() {
    // Фіксуємо кількість валентних кварків (2 up, 1 down)
    par_xuv->setMoment(-1, 2.0);
    par_xdv->setMoment(-1, 1.0);
    
    // Правило сум імпульсу: кварки
    // Оскільки S = 2(ubar + dbar) + s_plus, а s_plus = s + sbar,
    // загальний імпульс - це просто сума валентних кварків та загального моря S
    double xsumq = 0.0;
    xsumq += par_xuv->moment(0);
    xsumq += par_xdv->moment(0);
    xsumq += par_S->moment(0);
    
    // Залишок імпульсу віддаємо глюону
    par_xg->setMoment(0, 1.0 - xsumq);
  }

  // Обчислення кінцевих PDG частинок для даного x
  std::map<int,double> MSHT_Decomposition::xfxMap(double x) const {
    // 1. Отримуємо значення з наших параметризацій
    double xuv    = (*par_xuv)(x);
    double xdv    = (*par_xdv)(x);
    double S      = (*par_S)(x);
    double splus  = (*par_splus)(x);
    double sminus = (*par_sminus)(x);
    double ratio  = (*par_d_over_u)(x); // Це наше d_bar / u_bar
    double g      = (*par_xg)(x);

    // 2. Розв'язуємо алгебру для дивних кварків
    double sbar = (splus - sminus) / 2.0;
    double s    = (splus + sminus) / 2.0;

    // 3. Розв'язуємо алгебру для u_bar та d_bar
    // S = 2(ubar + dbar) + splus  =>  2(ubar + ratio*ubar) = S - splus
    double ubar = (S - splus) / (2.0 * (1.0 + ratio));
    double dbar = ubar * ratio;

    // 4. Формуємо повні u та d кварки
    double u = xuv + ubar;
    double d = xdv + dbar;

    // 5. Повертаємо мапу у форматі LHAPDF
    std::map<int,double> out = {
      {-6, 0.0},
      {-5, 0.0},
      {-4, 0.0},
      {-3, sbar},
      {-2, ubar},
      {-1, dbar},
      { 1, d},
      { 2, u},
      { 3, s},
      { 4, 0.0},
      { 5, 0.0},
      { 6, 0.0},
      {21, g}
    };
    
    return out;
  }
}