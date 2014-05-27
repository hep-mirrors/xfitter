#include "Dataset.h"

#include "CommandParser.h"

#include <TObjArray.h>
#include <TObjString.h>

#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

Dataset::Dataset(int dataindex, string dataname) : index(dataindex), name(dataname)
{
  string trimmed = name;
  trimmed.erase(trimmed.find_last_not_of(" ")+1, string::npos);
  if (trimmed.substr(0, 1) == " ")
    trimmed.erase(0, trimmed.find_first_not_of(" "));
  name = trimmed;
 }

void Subplot::AddPoint(map <string, float> fline)
{
  bins1.push_back(fline["lbin"]);
  bins2.push_back(fline["ubin"]);
  valx.push_back(fline["x"]);
  data.push_back(fline["data"]);
  uncorerr.push_back(fline["uncor"]);
  toterr.push_back(fline["tot"]);
  theory.push_back(fline["thorig"]);
  theoryshifted.push_back(fline["thmod"]);
  therrup.push_back(fline["therr+"]);
  therrdown.push_back(fline["therr-"]);
  pulls.push_back(fline["pull"]);
}

Subplot::Subplot(string plotoptions) :  xmin(0), xmax(0), yminr(0), ymaxr(0), extralabel(""), experiment(""), title(""), xlabel(""), ylabel(""), logx(false), logy(false)
{
  valid = false;
  //parse plot options
  //Rewrite in standard C parsing, avoid Root string parsing
  TString popt(plotoptions.c_str());
  TObjArray* array = popt.Tokenize("@");
  
  // first loop to detect x axis
  for (int i=0; i<array->GetEntries(); i++) 
    {
      TString str( ((TObjString*)array->At(i))->GetString().Data());
      if(str.BeginsWith("Xmin:")) 
	{
	  str.ReplaceAll("Xmin:","");
	  xmin = str.Atof();
	}
      else if(str.BeginsWith("Xmax:")) 
	{
	  str.ReplaceAll("Xmax:","");
	  xmax = str.Atof();
	}
    }
  for(int i=0; i<array->GetEntries(); i++) 
    {
      TString str( ((TObjString*)array->At(i))->GetString().Data());
      if(str.BeginsWith("ExtraLabel:")) 
	{
	  str.ReplaceAll("ExtraLabel:","");
	  extralabel = str.Data();
	}
      if(str.BeginsWith("Experiment:")) 
	{
	  str.ReplaceAll("Experiment:","");
	  experiment = str.Data();
	}
      if(str.BeginsWith("Title:")) 
	{
	  str.ReplaceAll("Title:","");
	  title = str.Data();
	}
      else if(str.BeginsWith("XTitle:")) 
	{
	  str.ReplaceAll("XTitle:","");
	  xlabel = str.Data();
	}
      else if(str.BeginsWith("YTitle:")) 
	{
	  str.ReplaceAll("YTitle:","");
	  ylabel = str.Data();
	}
      else if(str.BeginsWith("YminR:")) 
	{
	  str.ReplaceAll("YminR:","");
	  yminr = str.Atof();
	}
      else if(str.BeginsWith("YmaxR:")) 
	{
	  str.ReplaceAll("YmaxR:","");
	  ymaxr = str.Atof();
	}
      else if(str.BeginsWith("Xlog"))
	logx = true;
      else if(str.BeginsWith("Ylog"))
	logy = true;
    }
  delete array;
}

void Subplot::Init(string label, int dataindex, int subplotindex)
{
  //If no bin edges information is available, need to use TGraph for plotting
  maketgraph = false;
  for (vector <float>::iterator it = valx.begin(); it != valx.end(); it++)
    if (*it != 0)
      maketgraph = true;

  //check bins sanity
  vector <float> b1 = bins1;
  vector <float> b2 = bins2;
  vector<float>::iterator it1 = b1.begin();
  vector<float>::iterator it2 = b2.begin();
  if (b1.size() < 1)
      return;   //      cout << "zero bins for " << label << ", dataset: " << dataindex << ", subplot: " << subplotindex << ". skipping..." << endl;

  bool skip = false;
  for (; (it1+1) != b1.end(); it1++, it2++)
    if (*(it1+1) < *it2 || *it1 >= *(it1+1))
      skip = true;
  if (skip)
    {
      if (!maketgraph)
	{
	  cout << "bin inconsistency for " << label << ", dataset: " << dataindex << ", subplot: " << subplotindex << ". skipping..." << endl;
	  cout << "Cannot plot data, skipping" << endl;
	  return;
	}
    }


  //bins array
  float bin[bins1.size() + 1];
  vector <double> bins;

  if (maketgraph)
    {
      //make arbitrary bin edges
      vector <double> temp;
      for (vector <float>::iterator it = valx.begin(); it != valx.end(); it++)
	temp.push_back(*it);
      
      sort(temp.begin(), temp.end());

      //adjust x axis
      double xaxmin = *(temp.begin());
      double xaxmax = *(temp.end()-1);
      if (valx.size() > 1 && (*(temp.end()-1) - *(temp.begin()) > 0))
	{
	  if (logx)
	    {
	      double axislength = *(temp.end()-1) / *(temp.begin());
	      xaxmin = xaxmin / pow(axislength,1./20);
	      xaxmax = xaxmax * pow(axislength,1./20);
	    }
	  else
	    {
	      double axislength = *(temp.end()-1) - *(temp.begin());
	      xaxmin = xaxmin - axislength/20.;
	      xaxmax = xaxmax + axislength/20.;
	    }
	}
      else
	{
	  xaxmin *= 0.9999;
	  xaxmax *= 1.0001;
	}
      bins.push_back(xaxmin);
      if (valx.size() > 1)
	for (vector<double>::iterator it = temp.begin()+1; it != temp.end(); it++)
	  bins.push_back((*it + *(it-1))/2);
      bins.push_back(xaxmax);


      for (vector<double>::iterator it = bins.begin(); it != bins.end(); it++)
	bin[it-bins.begin()] = *it;
    }
  else
    {
      //fill empty gap
      int pos = 0;
      float bmin, bmax;
      while (pos != -1)
	{
	  pos = -1;
	  vector<float>::iterator it1 = bins1.begin();
	  vector<float>::iterator it2 = bins2.begin();
	  for (; (it1+1) != bins1.end(); it1++, it2++)
	    if (*(it1+1) != *it2 && *it2 < *(it1+1))
	      {
		pos = (it1 - bins1.begin()) + 1;
		bmin = *it2;
		bmax = *(it1+1);
	      }
	  if (pos != -1)
	    {
	      bins1.insert(bins1.begin()+pos, bmin);
	      bins2.insert(bins2.begin()+pos, bmax);
	      data.insert(data.begin()+pos, 0);
	      uncorerr.insert(uncorerr.begin()+pos, 0); 
	      toterr.insert(toterr.begin()+pos, 0); 
	      theory.insert(theory.begin()+pos, 0); 
	      theoryshifted.insert(theoryshifted.begin() +pos, 0); 
	      therrup.insert(therrup.begin()+pos, 0); 
	      therrdown.insert(therrdown.begin()+pos, 0); 
	      pulls.insert(pulls.begin()+pos, 0);
	    }
	}

      //make bins array
      int i = 0;
      for (vector<float>::iterator it = bins1.begin(); it != bins1.end(); it++)
	{
	  bin[i] = *it;
	  i++;
	}
      bin[i] = *(bins2.end()-1);
    }
  char hnm[300];
  sprintf (hnm, "data_%s_%d-%d", label.c_str(), dataindex, subplotindex);
  string hname = hnm;
  hdata = new TH1F((hname +"_data").c_str(), "", bins1.size(),  bin);
  hdatatot = new TH1F((hname + "_datatot").c_str(), "", bins1.size(),  bin);
  hth = new TH1F((hname + "_th").c_str(), "", bins1.size(),  bin);
  hthshift = new TH1F((hname + "_thshift").c_str(), "", bins1.size(),  bin);
  htherr = new TH1F((hname + "_therr").c_str(), "", bins1.size(),  bin);
  htherrup = new TH1F((hname + "_therrup").c_str(), "", bins1.size(),  bin);
  htherrdown = new TH1F((hname + "_therrdown").c_str(), "", bins1.size(),  bin);
  hpull = new TH1F((hname + "_pull").c_str(), "", bins1.size(),  bin);

  if (xmin != 0 && xmax != 0)
    {
      hdata->SetAxisRange(xmin, xmax);
      hdatatot->SetAxisRange(xmin, xmax);
      hth->SetAxisRange(xmin, xmax);
      hthshift->SetAxisRange(xmin, xmax);
      htherr->SetAxisRange(xmin, xmax);
      htherrup->SetAxisRange(xmin, xmax);
      htherrdown->SetAxisRange(xmin, xmax);
      hpull->SetAxisRange(xmin, xmax);
    }
  else
    {
      xmin = hdata->GetXaxis()->GetBinLowEdge(hdata->GetXaxis()->GetFirst());
      xmax = hdata->GetXaxis()->GetBinUpEdge(hdata->GetXaxis()->GetLast() - 1);
    }
      
  hdata->SetXTitle(xlabel.c_str());
  hdata->SetYTitle(ylabel.c_str());

  //set default labels
  if (xlabel == "" && ylabel == "")
    {
      hdata->SetXTitle("(Set XTitle:<label>)");
      hpull->SetXTitle("(Set XTitle:<label>)");
      hdata->SetYTitle("(Set YTitle:<label>)");
    }

  for (unsigned int b = 0; b < data.size(); b++)
    {
      hdata->SetBinContent(b + 1, data[b]);
      hdata->SetBinError(b + 1, uncorerr[b]);
      hdatatot->SetBinContent(b + 1, data[b]);
      hdatatot->SetBinError(b + 1, toterr[b]);
      hth->SetBinContent(b + 1, theory[b]);
      hthshift->SetBinContent(b + 1, theoryshifted[b]);
      htherr->SetBinContent(b + 1, theory[b] + (therrup[b] - therrdown[b]) / 2);
      htherr->SetBinError(b + 1, (therrup[b] + therrdown[b]) / 2);
      htherrup->SetBinContent(b + 1, theory[b] + therrup[b]);
      htherrdown->SetBinContent(b + 1, theory[b] - therrdown[b]);
      //invert pulls -> (theory - data)
      if (!opts.ratiototheory)
	hpull->SetBinContent(b + 1, -pulls[b]);
      else
	hpull->SetBinContent(b + 1, pulls[b]);
      hpull->SetBinError(b + 1, 0);
    }

  //Prepare ratio histograms
  r_th = (TH1F*)hth->Clone();
  r_thshift = (TH1F*)hthshift->Clone();
  r_therr = (TH1F*)htherr->Clone();
  r_therrup = (TH1F*)htherrup->Clone();
  r_therrdown = (TH1F*)htherrdown->Clone();

  //Select reference, data or theory
  //  TH1F * refdata;
  if (opts.ratiototheory)
    href = (TH1F*)hth->Clone();
  else
    href = (TH1F*)hdata->Clone();

  //Prevent double counting of errors in ratio
  for (int b = 1; b <= href->GetNbinsX(); b++)
    href->SetBinError(b, 0);

  /*
  r_th->Divide(refdata);
  r_thshift->Divide(refdata);
  r_therr->Divide(refdata);
  r_therrup->Divide(refdata);
  r_therrdown->Divide(refdata);

  for (int b = 1; b <= r_th->GetNbinsX(); b++)
    r_th->SetBinError(b, 0);
  for (int b = 1; b <= r_thshift->GetNbinsX(); b++)
    r_thshift->SetBinError(b, 0);
  for (int b = 1; b <= r_therr->GetNbinsX(); b++)
    r_therr->SetBinError(b, (r_therrup->GetBinContent(b) - r_therrdown->GetBinContent(b)) / 2 );
  for (int b = 1; b <= r_therrup->GetNbinsX(); b++)
    r_therrup->SetBinError(b, 0);
  for (int b = 1; b <= r_therrdown->GetNbinsX(); b++)
    r_therrdown->SetBinError(b, 0);
  */
  valid = true;
}

Data::Data(string dirname, string label)
{
  //Open datasets file (fittedresults.txt)
  char filename[300];
  sprintf (filename, "%s/fittedresults.txt", dirname.c_str());
  ifstream infile(filename);
  if (!infile.is_open()) //fittedresults.txt not found
    return;

  //Read datasets
  string line;
  int dtindex, nextdtindex;
  string name;
  float buffer;
  double bin1, bin2, data, uncorrerr, toterr, theory, theory_mod, therr_up, therr_down, pull;
  map <int, string> coltag;

  int Ndatasets;
  getline(infile, line);
  istringstream issnd(line);
  issnd >> Ndatasets; // Total number of data sets

  if (infile.eof() || Ndatasets == 0)
    return; //empty file, or no datasets found

  getline(infile, line);
  istringstream issdi(line);
  issdi >> nextdtindex;  //Dataset index

  //Loop on datasets
  while(!infile.eof()) 
    {
      dtindex = nextdtindex;

      //Read dataset name
      getline(infile, name);
      
      //Initialise new dataset
      Dataset dtset(dtindex, name);

      getline(infile, line);
      //Read plot options
      while(!infile.eof())
	{
	  if (line.find("Plot") == string::npos)
	    break;

	  //Read subplot index
	  TString temp(line.c_str());
	  TObjArray* array = temp.Tokenize("@");
	  temp.Form(((TObjString*)  array->At(0))->GetString().Data());
	  delete array;
	  temp.ReplaceAll("Plot","");
	  int iplot = temp.Atoi();
	  //End of reading subplotindex

	  dtset.subplots[iplot] = Subplot(line);
	  getline(infile, line);
	}

      //Read columns tags
      string col;
      int i = 0;
      istringstream iss(line);

      while (iss >> col)
	{
	  //	  coltag[i] = col;
	  //	  i++;
	}

      //hard coded patch, until the fittedresults.txt format is improved
      coltag[0] = "lbin";
      coltag[1] = "ubin";
      coltag[2] = "dummy";
      coltag[3] = "data";
      coltag[4] = "uncor";
      coltag[5] = "tot";
      coltag[6] = "thorig";
      coltag[7] = "thmod";
      coltag[8] = "therr+";
      coltag[9] = "therr-";
      coltag[10] = "pull";
      coltag[11] = "iset";
      coltag[12] = "iplot";
      coltag[13] = "x";
      //end of patch

      getline(infile, line);
      //Loop on data points
      while(!infile.eof())
	{
	  istringstream iss(line);
	  
	  //Read a line of data
	  int iplot;
	  map <string, float> fline;
	  int i = 0;
	  while (iss >> buffer)
	    {
	      if (coltag[i] == "iplot")
		iplot = (int) buffer;
	      else
		fline[coltag[i]] = buffer;
	      i++;
	      //patch, to be cleaned up when the fittedresults.txt format is improved
	      if (i == 12)
		{
		  string s;
		  iss >> s;		  
		  TString str(s.c_str());
		  TObjArray* array = str.Tokenize("/");
		  iplot = ((TObjString*) array->At(0))->GetString().Atoi();
		  fline["x"] = ((TObjString*) array->At(1))->GetString().Atof();
		  break;
		}
	      //end of patch
	    }
	  if (i == 1) //New dataset
	    {
	      //Set dataset index for the incoming dataset
	      nextdtindex = (int) buffer;
	      break; 
	    }
	  //Add point to subplot
	  dtset.subplots[iplot].AddPoint(fline);
	  getline(infile, line);
	}//End loop on data points
      for (map <int, Subplot>::iterator sit = dtset.subplots.begin(); sit != dtset.subplots.end(); sit++)
	{
	  (*sit).second.Init(label, dtindex, (*sit).first);
	}
      datamap[dtindex] = dtset;
    }//End loop on datasets
}
