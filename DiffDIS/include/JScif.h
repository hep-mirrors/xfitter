/*
  F2_heavy via
  gamma* parton -> Q + Qbar + X folded with the parton densities
*/


#ifdef __cplusplus
extern "C" {
#endif

void JSsetFact(double qfac, double mfac);
void JSsetContr(int born, int gluon, int quark);
void JSsetSF(int Ftype, int q_id, double mq2);
void JSsetParams(int Ftype, int q_id, double mq, double qfac, double mfac);
double JackSmith(double x, double QQ);

#ifdef __cplusplus
}
#endif
