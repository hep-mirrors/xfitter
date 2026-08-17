#include "ChebyshevFlex_PdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <iostream>

namespace xfitter {

extern "C" ChebyshevFlex_PdfParam* create(const char* name) {
  return new ChebyshevFlex_PdfParam(name);
}

void ChebyshevFlex_PdfParam::atStart() {
  using namespace std;
  BasePdfParam::atStart();
  const size_t n = getNPar();
  
  if (n < 4) {
    cerr << "[ERROR] Too few parameters given to ChebyshevFlex parameterisation \"" 
         << _name << "\", expected at least 4, got " << n << endl;
    hf_errlog(18120702, "F: Wrong number of parameters for ChebyshevFlex, see stderr");
  }
}

double ChebyshevFlex_PdfParam::operator()(double x) const {
  const unsigned int npar = getNPar();
  
  // Витягуємо параметр k (тепер він на нульовому індексі)
  double k = *pars[0];

  // Базовий множник тепер використовує індекси 1, 2 та 3
  double power = (*pars[1]) * pow(x, (*pars[2])) * pow((1.0 - x), (*pars[3]));
  
  // Якщо поліномів немає (передали лише 4 параметри: k, A, B, C)
  if (npar == 4) return power;

  // Нелінійне перетворення координати з використанням параметра k
  double y = 1.0 - 2.0 * pow(x, k);

  // Ініціалізуємо суму одиницею (зафіксований T_0)
  double cheb_sum = 1.0; 
  
  // Рекурентні змінні
  double t_prev = 1.0; // T_0(y)
  double t_curr = y;   // T_1(y)
  
  // pars[4] залишається коефіцієнтом для T_1
  cheb_sum += (*pars[4]) * t_curr;
  
  // Додаємо всі наступні члени починаючи з T_2
  for (unsigned int i = 5; i < npar; i++) {
    // Швидке обчислення наступного полінома Чебишова
    double t_next = 2.0 * y * t_curr - t_prev; 
    
    cheb_sum += (*pars[i]) * t_next;
    
    // Зсув для наступної ітерації циклу
    t_prev = t_curr;
    t_curr = t_next;
  }

  return power * cheb_sum;
}

} // namespace xfitter