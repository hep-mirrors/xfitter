#include "ChebyshevGluon_PdfParam.h"
#include <cmath>
#include <iostream>

namespace xfitter {

extern "C" ChebyshevGluon* create(const char* name) {
  return new ChebyshevGluon(name);
}

ChebyshevGluon::ChebyshevGluon(const char* name) : BasePdfParam(name) {}

const char* ChebyshevGluon::getClassName() const {
  return "ChebyshevGluon";
}

double ChebyshevGluon::operator()(double x) const {
  // Використовуємо вбудований метод xFitter для отримання розміру
  const unsigned int npar = getNPar();
  
  // Захист від помилок: перевіряємо, чи передано мінімум 7 базових параметрів
  if (npar < 7) {
    std::cerr << "[ERROR] Too few parameters given to ChebyshevGluon parameterisation \"" 
              << _name << "\", expected at least 7, got " << npar << std::endl;
    return 0.0; 
  }

  // Отримуємо значення, розіменовуючи вказівники (*pars[i])
  // 0. Кінематичний параметр
  double k = *pars[0];

  // 1. Окремий доданок (без Чебишова)
  double Ag_minus     = *pars[1];
  double etag_minus   = *pars[2];
  double deltag_minus = *pars[3];
  double term_separate = Ag_minus * std::pow(1.0 - x, etag_minus) * std::pow(x, deltag_minus);

  // 2. Базові параметри для Чебишова
  double Ag     = *pars[4];
  double etag   = *pars[5];
  double deltag = *pars[6];

  // 3. Рахуємо кінематику та суму поліномів
  double y = 1.0 - 2.0 * std::pow(x, k);
  double cheb_sum = 1.0; 

  // Якщо у конфізі є коефіцієнти Чебишова (тобто npar > 7)
  if (npar > 7) {
    double t_prev = 1.0; // T_0(y)
    double t_curr = y;   // T_1(y)
    
    // pars[7] - це коефіцієнт для T_1
    cheb_sum += (*pars[7]) * t_curr;
    
    // Додаємо всі наступні члени починаючи з T_2 (індекси 8 і далі)
    for (unsigned int i = 8; i < npar; i++) {
      double t_next = 2.0 * y * t_curr - t_prev; 
      
      cheb_sum += (*pars[i]) * t_next;
      
      // Зсув для наступної ітерації
      t_prev = t_curr;
      t_curr = t_next;
    }
  }

  // 4. Збираємо перший терм (база * суму Чебишова)
  double term_cheb = Ag * std::pow(1.0 - x, etag) * std::pow(x, deltag) * cheb_sum;

  // 5. Повертаємо повне значення функції
  return term_cheb + term_separate;
}

} // namespace xfitter