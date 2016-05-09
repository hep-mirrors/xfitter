#include "CuteInterface.h"

#include <xfitter_cpp.h>

#include <vector>
#include <string.h>
#include <algorithm>

using namespace std;

//extern "C" void hf_errlog_(const int &id, const char *text, int); 

double quadchi2(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double accuracy)
{
  return 0.;
}

//Constructor
Cute::Cute(string file): infile(file)
{
  gsl_set_error_handler_off();
  plot.init(file.c_str());
  //store nominal accuracy
  accuracy = plot.p.gen;
}

void Cute::SetBins(vector <double> lowedge, vector <double> upedge, int bd)
{
  values.resize(lowedge.size());

  //points-per-bin
  bindensity = bd;

  //create a number of qt points corresponding to bindensity * nbins
  vector<double>::iterator itl = lowedge.begin();
  vector<double>::iterator itu = upedge.begin();
  for (; itl != lowedge.end(); itl++, itu++)
    {
      binwidth.push_back(*itu - *itl);
      for (int i = 0; i < bindensity; i++)
	points[*itl + double(i) * (*itu - *itl)/double(bindensity)] = 0;
    }
  points[*(upedge.end()-1)] = 0;

  plot.q.clear();
  for (map<double, double>::iterator it = points.begin(); it != points.end(); it++)
    if (it->first > 0.05) //lower limit for cross section is 0.05 GeV, set to 0 below such limit
      plot.q.push_back(it->first);
  plot.PP = new PlotPoint(plot.p,plot.y.size(),plot.matchingScheme,plot.toPlot);
}

void Cute::Calculate(const double mur, const double muf)
{
  //plot.PP = new PlotPoint(plot.p,plot.y.size(),plot.matchingScheme,plot.toPlot);
  plot.p.Lambda=max(0., c_lambda_.lambdanp_);
  plot.PP->init(plot.p);

  //recalculate qstar at current value of alphas, use accuracy < 1e-05
  plot.p.gen = min(1e-05, plot.p.gen);
  plot.PP->init(plot.p);
  plot.qstar = plot.p.q_star();

  plot.p.gen = accuracy;
  plot.PP->init(plot.p);

  cout << "qstar: " << plot.qstar << " alphas: " << appl_fnalphas_(91.1876) << " lambda: " << plot.p.Lambda << endl;

  //loop on qt points
  for(uint i = 0; i < plot.q.size(); ++i)
    {
      //calculate cross section at qt point plot.q[i]
      plot.p.qT=plot.q[i];

      //set fixed scale
      if (plot.mu.size() > 0)
	plot.p.mu=mur*plot.mu[0];
      else //set qt dependent renormalisation scale
	plot.p.mu=mur*(plot.p.qT+plot.qstar);

      //calculate cross section
      //loop on y points
      if(plot.y.size() > 0)
	for(uint j = 0; j < plot.y.size(); ++j)
	  {
	    plot.PP->values[e_y]=plot.p.y=plot.y[j];
	    plot.p.update(plot.p.qT,plot.p.mu,plot.p.y);
	    plot.PP->update();
	  }
      else
	{
	  plot.p.update(plot.p.qT,plot.p.mu,plot.p.y);
	  plot.PP->update();
	}

      //cout << "qT = " << plot.q[i] << "; mu = " <<  plot.p.mu << "; xs = " <<  plot.PP->values[e_erg] << endl;

      //Store cross section
      points[plot.q[i]] = plot.PP->values[e_erg];
    }

  //free pointer
  //  delete plot.PP;

  /*
  for (map<double, double>::iterator pit = points.begin(); pit != points.end(); pit++)
    cout << pit->first << "  " << pit->second << endl;
  */

  //Check for jumps in the cross section and recalculate with better accuracy


  //Check for maximum one changes of sign in the first derivative of the positive and negative cross sections
  bool checkderivative = false;
  int iteration = 0;
  while (plot.p.gen > 1e-07 && !checkderivative)
    {
      checkderivative = true;
      iteration++;
      plot.p.gen = plot.p.gen / 10.;
      plot.PP->init(plot.p);
      vector <int> poslocextr;
      vector <int> neglocextr;
      vector <int> locextr;
      //loop on qt points, add a leading and a trailing 0 to catch issues in the first and last points
      for(uint i = -1; (i+2) <= plot.q.size(); ++i)
	{
	  double ap, bp, cp;
	  if (i == -1)
	    ap = 0;
	  else
	    ap = max(0., points[plot.q[i]]);
	  bp = max(0., points[plot.q[i+1]]);
	  if ((i+2) == plot.q.size())
	    cp = 0;
	  else
	    cp = max(0., points[plot.q[i+2]]);

	  double am, bm, cm;
	  if (i == -1)
	    am = 0;
	  else
	    am = min(0., points[plot.q[i]]);
	  bm = min(0., points[plot.q[i+1]]);
	  if ((i+2) == plot.q.size())
	    cm = 0;
	  else
	    cm = min(0., points[plot.q[i+2]]);

	  //check if b is a local maximum or minimum
	  if ((cp-bp)*(bp-ap) < 0)
	    {
	      poslocextr.push_back(i+1);
	      locextr.push_back(i+1);
	    }
	  if ((cm-bm)*(bm-am) < 0)
	    {
	      neglocextr.push_back(i+1);
	      locextr.push_back(i+1);
	    }

	}

      if (poslocextr.size() > 1 || neglocextr.size() > 1)
	{
	  cout << "failed local extremes test, positive: " << poslocextr.size() << "; negative: " << neglocextr.size() << endl;
	  checkderivative = false;
	}
      
      if (!checkderivative)
	{
	  int id = 0505201401;
	  char text[] = "W: Cute accuracy warning";
	  int textlen = strlen(text);
	  hf_errlog_(id, text, textlen);
	  cout << "CuteInterface, recalculate points with better accuracy: " << plot.p.gen << endl;
	}

      if (iteration <= 1)
	checkderivative = false;

      //revaluate extremal points
      if (!checkderivative)
	{
	  sort (locextr.begin(), locextr.end());
	  vector<int>::iterator it = unique(locextr.begin(), locextr.end());
	  locextr.resize(distance(locextr.begin(), it));
	  for (vector<int>::iterator it = locextr.begin(); it != locextr.end(); it++)
	    {
	      plot.p.qT=plot.q[*it];
	      if (plot.mu.size() > 0)
		plot.p.mu=mur*plot.mu[0];
	      else
		plot.p.mu=mur*(plot.p.qT+plot.qstar);
	      plot.p.update(plot.p.qT,plot.p.mu,plot.p.y);
	      plot.PP->update();
	      double before = points[plot.q[*it]];
	      points[plot.q[*it]] = plot.PP->values[e_erg];
	      cout << "qT = " << plot.q[*it] << "; mu = " <<  plot.p.mu << "; old xs = " << before << "; new xs = " <<  plot.PP->values[e_erg] << endl;
	    }
	}
    }
  if (plot.p.gen <= 1e-07)
    {
      int id = 0607201401;
      char text[] = "E: Cute numerical error";
      int textlen = strlen(text);
      hf_errlog_(id, text, textlen);
      cout << "Error in CuteInterface, could not reach accuracy: " << plot.p.gen << endl;
      cout << "lambda = " <<  plot.p.Lambda << "; alphas = " <<  appl_fnalphas_(91.1876) << endl;
    }
  plot.p.gen = accuracy;

  /*
  //Change the following to a chi2 test of a quadratic fit to 4 points
  bool checkchi2 = false;
  while (plot.p.gen > 1e-07 && !checkchi2)
    {
      checkchi2 = true;
      plot.p.gen = plot.p.gen / 10.;
      plot.PP->init(plot.p);
      for(uint i = 0; (i+3) < plot.q.size(); ++i)
	{
	    double check = true;
	    double a = points[plot.q[i]];
	    double b = points[plot.q[i+1]];
	    double c = points[plot.q[i+2]];
	    double d = points[plot.q[i+4]];

	    //quadratic regression test
	    double chi2 = quadchi2(plot.q[i], a, plot.q[i+1], b,  plot.q[i+2], c, plot.q[i+3], d, accuracy);
	    if (chi2 > 5)
	      {
		cout << "failed quadratic regression test: " << chi2 << endl;
		checkchi2 = false;
	      }
	      
	    //revaluate the four points
	    if (!checkchi2)
	      {
		int id = 0505201401;
		char text[] = "W: Cute accuracy warning";
		int textlen = strlen(text);
		hf_errlog_(id, text, textlen);
		cout << "CuteInterface, recalculate points with better accuracy: " << plot.p.gen << endl;
		for(uint j = i; j <= i+3; ++j)
		  {
		    plot.p.qT=plot.q[j];
		    if (plot.mu.size() > 0)
		      plot.p.mu=mur*plot.mu[0];
		    else
		      plot.p.mu=mur*(plot.p.qT+plot.qstar);
		    plot.p.update(plot.p.qT,plot.p.mu,plot.p.y);
		    plot.PP->update();
		    double before = points[plot.q[j]];
		    points[plot.q[j]] = plot.PP->values[e_erg];
		    cout << "qT = " << plot.q[j] << "; mu = " <<  plot.p.mu << "; old xs = " << before << "; new xs = " <<  plot.PP->values[e_erg] << endl;
		  }
	      }
	}
    }
  if (plot.p.gen <= 1e-07)
    {
      int id = 0607201401;
      char text[] = "E: Cute numerical error";
      int textlen = strlen(text);
      hf_errlog_(id, text, textlen);
      cout << "Error in CuteInterface, could not reach accuracy: " << plot.p.gen << endl;
      cout << "lambda = " <<  plot.p.Lambda << "; alphas = " <<  appl_fnalphas_(91.1876) << endl;
    }
  plot.p.gen = accuracy;
  */

  //convert points into bins
  map<double, double>::iterator pit = points.begin();
  vector<double>::iterator itbw = binwidth.begin();
  for (vector<double>::iterator it = values.begin(); it != values.end(); it++, itbw++)
    {
      *it = pit->second/2.; 
      pit++;
      for (int i = 1; i <= (bindensity - 1); i++)
	{
	  *it += pit->second;
	  pit++;
	}
      *it += pit->second/2.; 

      *it *= 1./double(bindensity);
    }
  return;
}
