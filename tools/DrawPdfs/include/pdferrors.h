#ifndef pdferrors_h
#define pdferrors_h

#include <vector>

using namespace std;

//Monte Carlo PDF errors
extern double mean(vector <double> xi);
extern double rms(vector <double> xi);
extern double median(vector <double> xi);
extern double cl(int sigma);
extern double delta(vector <double> xi, double central, double ConfLevel);
extern void deltaasym(vector <double> xi, double central, double& delta_p, double& delta_m, double ConfLevel);
//Hessian PDF errors
extern double ahessdelta(vector <double> xi);
extern void ahessdeltaasym(vector <double> xi, double& delta_p, double& delta_m);
extern double shessdelta(vector <double> xi);
//Experimental+Model+Parametrisation PDF errors
extern void vardeltaasym(vector <double> xi, int npar, double& delta_p, double& delta_m);
extern void deltaenvelope(vector <double> xi, double& delta_p, double& delta_m);

#endif
