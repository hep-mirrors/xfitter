#include "Hathor.h"
#include "SgTop.h"
#include "cstring"
#include "xfitter_cpp.h"

extern "C" {
  double sd2_(double* acc, double (*f)(double*), void (*r)(int, double*, double*, double*));
}

class IHathorGenericIntegrator {
public:
  virtual ~IHathorGenericIntegrator() = default;
  virtual double getXsection(double m, double mur, double muf, const std::string& integrator) = 0;
  virtual void setPartonicEnergy(const double x[]) = 0;
  virtual void evaluatePDFs(const double h1[], const double h2[], const double h2left[] = 0, const double h2right[] = 0) = 0;
  virtual void evaluateScalingFunctions() = 0;
  virtual double evaluateIntegral(double as, double wgt) = 0;
  virtual void setSwq(double x) = 0;
  virtual void setAlpha(double x) = 0;
  virtual void setCkmMatrix(const double ckm[3][3]) = 0;
  virtual void getCkmMatrix(double ckm[3][3]) = 0;
  virtual void PrintCkmMatrix() = 0;
  virtual void setParticle(SgTop::PARTICLE p) = 0;
  virtual double getAlphas(double mur) = 0;
  virtual void setColliderType(AbstractHathor::COLLIDERTYPE type) = 0;
  virtual void setScheme(unsigned int newscheme) = 0;
  virtual void sethc2(double hcq_) = 0;
  virtual void setSqrtShad(double sqrts) = 0;
  virtual void setPrecision(int n) = 0;
  virtual void getResult(int pdfset, double &integral, double &err) = 0;
  virtual void getResult(int pdfset, double &integral, double &err, double & chi2a) = 0;
};

class HathorLikeSgTop : public Hathor {
  public:
    HathorLikeSgTop(Pdf & pdf_) : Hathor(pdf_) {;};
    // dummy staff which is never called for Hathor
    void setCkmMatrix(const double ckm[3][3]) {;};
    void getCkmMatrix(double ckm[3][3]) {;};
    void PrintCkmMatrix() {;};
    void setParticle(SgTop::PARTICLE p) {;};
};

template <class Integrand> 
class HathorGenericIntegrator : public IHathorGenericIntegrator, public Integrand {
  public:
    HathorGenericIntegrator(Pdf & pdf_) : Integrand(pdf_) {
      hathor = this;
    };
    void setPartonicEnergy(const double x[]) override {
      Integrand::setPartonicEnergy(x);
    };
    void evaluatePDFs(const double h1[], const double h2[], 
			    const double h2left[] = 0, 
			    const double h2right[] = 0) override {
      Integrand::evaluatePDFs(h1, h2, h2left, h2right);
    };
    void evaluateScalingFunctions() override {
      Integrand::evaluateScalingFunctions();
    };
    double evaluateIntegral(double as, double wgt) override {
      return Integrand::evaluateIntegral(as, wgt);
    };
    void setSwq(double x) override {
      Integrand::setSwq(x);
    };
    void setAlpha(double x) override {
      Integrand::setAlpha(x);
    };
    void setCkmMatrix(const double ckm[3][3]) override {
      Integrand::setCkmMatrix(ckm);
    };
    void getCkmMatrix(double ckm[3][3]) override {
      Integrand::getCkmMatrix(ckm);
    };
    void PrintCkmMatrix() override {
      Integrand::PrintCkmMatrix();
    };
    void setParticle(SgTop::PARTICLE p) override {
      Integrand::setParticle(p);
    };
    double getAlphas(double mur) override {
      return Integrand::getAlphas(mur);
    }
    void setColliderType(AbstractHathor::COLLIDERTYPE type) override {
      Integrand::setColliderType(type);
    }
    void setScheme(unsigned int newscheme) override {
      Integrand::setScheme(newscheme);
    }
    void sethc2(double hcq_) override {
      Integrand::sethc2(hcq_);
    }
    void setSqrtShad(double sqrts) override {
      Integrand::setSqrtShad(sqrts);
    }
    void setPrecision(int n) override {
      Integrand::setPrecision(n);
    }
    void getResult(int pdfset, double &integral, double &err) override {
      Integrand::getResult(pdfset, integral, err);
    }
    void getResult(int pdfset, double &integral, double &err, double & chi2a) override {
      Integrand::getResult(pdfset, integral, err, chi2a);
    }

    double getXsection(double m, double mur, double muf, const std::string& integrator);
    // number of cross section evaluation calls
    unsigned long _ncalls;
  private:
    static HathorGenericIntegrator* hathor;
    double getXsectionVegas(double m, double mur, double muf);
    double getXsectionSd2(double m, double mur, double muf);
    static double callback(double x[]) {
      const double wgt = 1.0;
      double res[1];
      double ret = hathor->f(x, wgt, res);
      return ret;
    }
    static void callback_region(int ll, double* xx, double* aa, double* bb)
    {
      const double del = 1e-7;
      aa[0] = del;
      bb[0] = 1.0 - del;
      aa[1] = del;
      bb[1] = 1.0 - del;
    }
    // private methods taken from AbstractHathor.h
    double f(const double x[], const double wgt, double res[]);
    static inline void ChargeConjugation(double h[13]) {
      // Charge conjugate proton pdf's to obtain anti-proton pdf's 
      double tmp;
      for (int i = hathor->DOWN; i < hathor->ADOWN; i++) {
        tmp = h[i];
        h[i] = h[i + hathor->ADOWN-hathor->DOWN];
        h[i + hathor->ADOWN-hathor->DOWN] = tmp;
      }
    }
};

template <class Integrand>
HathorGenericIntegrator<Integrand>* HathorGenericIntegrator<Integrand>::hathor = nullptr;

template <class Integrand>
double HathorGenericIntegrator<Integrand>::getXsection(double m_, double mur_, double muf_, const std::string& integrator) {
  _ncalls = 0;
  double result = 0.;
  if (integrator == "vegas") {
    result = getXsectionVegas(m_, mur_, muf_);
  }
  else if (integrator == "sd2") {
    result = getXsectionSd2(m_, mur_, muf_);
  }
  else {
    char str[256];
    sprintf(str, "F: unknow integrator %s for HATHOR", integrator.c_str());
    hf_errlog_(26060301, str, strlen(str));
  }
  return result;
}

template <class Integrand>
double HathorGenericIntegrator<Integrand>::getXsectionVegas(double m_, double mur_, double muf_) {
  hathor->m = m_; hathor->mur = mur_; hathor->muf = muf_;
  hathor->update();
  
  Vegas<AbstractHathor> integrator(*this, hathor->pdfmax);
  //integrator.setnprn(1);
  integrator.vegas();
  
  for (int i = 0; i < hathor->pdfmax; i++) {
    hathor->integral[i] = integrator.avgi_[i];
    hathor->error[i] = integrator.sd_[i];
    hathor->chi2a[i] = integrator.chi2a_[i];
  }
  
  // reset pdf member
  if (hathor->pdfmax > 1)
    hathor->pdf.InitMember(0);
  
  return integrator.avgi;
}

template <class Integrand>
double HathorGenericIntegrator<Integrand>::getXsectionSd2(double m_, double mur_, double muf_) {
  hathor->m = m_; hathor->mur = mur_; hathor->muf = muf_;
  hathor->update();
  
  double wgt;
  double x[MAXDIM];
  double res[FDIMMAX];
  double acc = 0.;
  double s0=sd2_(&acc, callback, callback_region);
  acc = s0 * 1. / this->calls;
  double s1=sd2_(&acc, callback, callback_region);
  //printf("getXsection _ncalls, s0, s1 = %ld, %f, %f\n", _ncalls,s0,s1);
  
  for (int i = 0; i < hathor->pdfmax; i++) {
    hathor->integral[i] = s1;//integrator.avgi_[i];
    hathor->error[i] = 0.;//integrator.sd_[i];
    hathor->chi2a[i] = 0.;//integrator.chi2a_[i];
  }
  
  // reset pdf member
  if (hathor->pdfmax > 1)
    hathor->pdf.InitMember(0);
  
  return s1;
}

template <class Integrand>
double HathorGenericIntegrator<Integrand>::f(const double x[], const double wgt, double res[]) {
   _ncalls++;

   /*
   *  Differential cross section including the parton luminosities.
   */

  
  // calculate kinematic variables
  hathor->setPartonicEnergy(x);
  
  if ((hathor->delta > 0) && (hathor->x2 + hathor->delta > 1.))
    return 0.;
  
  // evaluate partonic cross sections
  hathor->evaluateScalingFunctions();
  
  double hi[13], hj[13], hjleft[13], hjright[13];
  // loop over pdf set members
  for(int ipdf = 0; ipdf < hathor->pdfmax; ipdf++) {
    if (hathor->pdfmax > 1)
      hathor->pdf.InitMember(ipdf);

    hathor->pdf.GetPdf(hathor->x1, hathor->muf, hi);
    hathor->pdf.GetPdf(hathor->x2, hathor->muf, hj);

    if (hathor->delta > 0) {
      hathor->pdf.GetPdf(hathor->x2-hathor->delta, hathor->muf, hjleft);
      hathor->pdf.GetPdf(hathor->x2+hathor->delta, hathor->muf, hjright);
    }
    
    if (hathor->collider == hathor->PPBAR) {
      ChargeConjugation(hj);
      if (hathor->delta > 0) {
	ChargeConjugation(hjleft);
	ChargeConjugation(hjright);
      }
    }
    
    // evaluate parton fluxes
    if (hathor->delta > 0) {
      hathor->evaluatePDFs(hi, hj, hjleft, hjright);
    } else {
      hathor->evaluatePDFs(hi, hj);
    }
    
    // calculate full integral
    res[ipdf] = hathor->evaluateIntegral(hathor->as_pdf[ipdf], wgt);
  }
  
  return res[0];
}
