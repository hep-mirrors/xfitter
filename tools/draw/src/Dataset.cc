#include "Dataset.h"
#include "Outdir.h"
#include "CommandParser.h"
#include "pdferrors.h"

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

Subplot::Subplot(string plotoptions) :  xmin(0), xmax(0), yminr(0), ymaxr(0), ymin(0), ymax(0), extralabel(""), experiment(""), title(""), xlabel(""), ylabel(""), logx(false), logy(false)
{
  hastherr = false;
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
      if(str.BeginsWith("Lumi:"))
        {
          str.ReplaceAll("Lumi:","");
          lumilabel = str.Data();
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
      else if(str.BeginsWith("Ymin:"))
        {
          str.ReplaceAll("Ymin:","");
          ymin = str.Atof();
        }
      else if(str.BeginsWith("Ymax:"))
        {
          str.ReplaceAll("Ymax:","");
          ymax = str.Atof();
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
  hdata      = new TH1F((hname + "_data"     ).c_str(), "", bins1.size(), bin);
  hdatatot   = new TH1F((hname + "_datatot"  ).c_str(), "", bins1.size(), bin);
  hth        = new TH1F((hname + "_th"       ).c_str(), "", bins1.size(), bin);
  hthshift   = new TH1F((hname + "_thshift"  ).c_str(), "", bins1.size(), bin);
  htherr     = new TH1F((hname + "_therr"    ).c_str(), "", bins1.size(), bin);
  htherrup   = new TH1F((hname + "_therrup"  ).c_str(), "", bins1.size(), bin);
  htherrdown = new TH1F((hname + "_therrdown").c_str(), "", bins1.size(), bin);
  hpull      = new TH1F((hname + "_pull"     ).c_str(), "", bins1.size(), bin);

  if (xmin == 0 && xmax == 0)
    {
      xmin = hdata->GetXaxis()->GetBinLowEdge(hdata->GetXaxis()->GetFirst());
      xmax = hdata->GetXaxis()->GetBinUpEdge(hdata->GetXaxis()->GetLast());
    }

  //Compute range for histograms
  lowrange = max(1, hdata->GetXaxis()->FindFixBin(xmin));
  uprange = min(hdata->GetNbinsX(), hdata->GetXaxis()->FindFixBin(xmax));
  hdata->GetXaxis()->SetRange(lowrange, uprange);
  hdatatot->GetXaxis()->SetRange(lowrange, uprange);
  hth->GetXaxis()->SetRange(lowrange, uprange);
  htherr->GetXaxis()->SetRange(lowrange, uprange);
  htherrup->GetXaxis()->SetRange(lowrange, uprange);
  htherrdown->GetXaxis()->SetRange(lowrange, uprange);
  hpull->GetXaxis()->SetRange(lowrange, uprange);

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

  for (unsigned int b = 0; b < data.size(); b++)
    if (therrup[b] != 0 || therrdown[b] != 0)
      hastherr = true;

  //Prepare ratio histograms
  r_th = (TH1F*)hth->Clone();
  r_thshift = (TH1F*)hthshift->Clone();
  r_therr = (TH1F*)htherr->Clone();
  r_therrup = (TH1F*)htherrup->Clone();
  r_therrdown = (TH1F*)htherrdown->Clone();

  //Select reference, data or theory
  if (opts.ratiototheory)
    href = (TH1F*)hth->Clone();
  else
    href = (TH1F*)hdata->Clone();

  //Prevent double counting of errors in ratio
  for (int b = 1; b <= href->GetNbinsX(); b++)
    href->SetBinError(b, 0);

  //Ratio histograms are evaluated in plotting in order to use a common reference

  valid = true;
}

vector <string> readLineFiles (vector <ifstream*> infiles)
{
  vector <string> lines;
  for (int i=0; i<infiles.size(); i++)
    {
      string line;
      getline(*infiles[i],line);
      lines.push_back(line);
    }
  return lines;
}

void getTheoryShift (vector<pdfshift> pdfshifts, vector <vector <double> > cor_matrix, pdferr err, vector <string> lines, bool scale68, double& corSum, double& errplus, double& errminus)
{
  // Decode string thing
  int N = ( err == AsymHess ) ? 2*pdfshifts.size()+1 : pdfshifts.size()+1;

  if (N == 1) return;

  vector <double> val;

  for (int i=0; i<N; i++) {
    istringstream iss(lines[i]);
    // Hardwire !!! //
    double buffer;
    for (int j =0; j<7; j++)
      iss >> buffer;
    val.push_back(buffer);
  }

  vector <double> xi;
  corSum = 0;
  double cent = val[0];
  xi.push_back(cent);
  if ( err == AsymHess )
    {
      for ( int i = 0; i<pdfshifts.size(); i++ )
        {
          double plus, minus;
          if (scale68)
            {
              plus  = (val[i*2+1] - cent) / 1.645;
              minus = (val[i*2+2] - cent) / 1.645;
            }
          else
            {
              plus  = val[i*2+1] - cent;
              minus = val[i*2+2] - cent;
            }
          double valShift = pdfshifts[i].val;
          double errShift = pdfshifts[i].err;
          //compute shifted central value
          double cor = 0.5*(plus - minus)*valShift + 0.5*(plus+minus)*valShift*valShift;
          corSum += cor;
          //compute reduced uncertainties
          xi.push_back(plus*errShift+cent);
          xi.push_back(minus*errShift+cent);
        }
      ahessdeltaasym(xi, errplus, errminus, cor_matrix);
    }
  else if ( err == SymHess )
    {
      for ( int i = 0; i < pdfshifts.size(); i++ )
        {
          double plus  = val[i+1] - cent;
          double valShift = pdfshifts[i].val;
          double errShift = pdfshifts[i].err;
          //compute shifted central value
          double cor = plus*valShift;
          corSum += cor;
          //compute reduced uncertainties
          xi.push_back(plus*errShift+cent);
        }
      errplus = errminus = shessdelta(xi, cor_matrix);
    }
}

void getTheoryReweight (vector<double> weights, vector <string> lines, double& val, double& err)
{
  val = 0;
  err = 0;

  // Decode string thing
  int N = lines.size();

  if (N == 0)
    return;

  vector <double> xi;

  for (int i = 0; i < N; i++) {
    istringstream iss(lines[i]);
    // Hardwire !!! //
    double buffer;
    for (int j =0; j<7; j++)
      iss >> buffer;
    xi.push_back(buffer);
  }

  val = mean(xi, weights);
  err = rms(xi, weights);
}

Data::Data(string dirname, string label)
{
  //Open datasets file (fittedresults.txt)
  char filename[300];
  if (outdirs[label].IsMCreplica())
    sprintf (filename, "%s/fittedresults.txt", outdirs[label].dirlist.begin()->c_str());
  else
    sprintf (filename, "%s/fittedresults.txt", dirname.c_str());

  ifstream infile(filename);
  if (!infile.is_open())
    {
      cout << "Error " << filename << " not found " << endl;
      return;
    }

  //open files for shifted or reweighted fittedresults.txt
  int nfiles = 0;
  pdferr err = pdfmap[label].err;
  if (outdirs[label].IsProfiled())
    {
      if (err == AsymHess)
        nfiles = pdfmap[label].pdfshifts.size()*2+1;
      else if (err == SymHess)
        nfiles = pdfmap[label].pdfshifts.size()+1;
    }
  if (outdirs[label].IsReweighted())
    if (err == MC)
      nfiles = pdfmap[label].mcw.size();

  // reset if shifts is empty.
  if (nfiles == 1) {nfiles = 0;}

  vector <ifstream*> infiles;
  for ( int i = 0; i < nfiles; i++)
    {
      sprintf (filename, "%s/fittedresults.txt_set_%04i", dirname.c_str(), i);
      infiles.push_back(new ifstream(filename));
      if (!infiles[i]->is_open())
        {
          cout << "Error " << filename << " not found " << endl;
          return;
        }
    }

  if (outdirs[label].IsMCreplica())
    for (vector <string>::iterator it = outdirs[label].dirlist.begin(); it != outdirs[label].dirlist.end(); it++)
      {
        string fname = (*it) + "/fittedresults.txt";
        infiles.push_back(new ifstream(fname.c_str()));
        if (!infiles[infiles.size()-1]->is_open())
          {
            cout << "Error " << filename << " not found " << endl;
            return;
          }
      }

  //Read datasets
  string line;
  vector <string> lines;
  lines.reserve(nfiles);

  int dtindex, nextdtindex;
  string name;
  float buffer;
  double bin1, bin2, data, uncorrerr, toterr, theory, theory_mod, therr_up, therr_down, pull;
  map <int, string> coltag;

  int Ndatasets;
  getline(infile, line); lines = readLineFiles(infiles);

  istringstream issnd(line);
  issnd >> Ndatasets; // Total number of data sets

  if (infile.eof() || Ndatasets == 0)
    return; //empty file, or no datasets found

  getline(infile, line); lines = readLineFiles(infiles);
  istringstream issdi(line);
  issdi >> nextdtindex;  //Dataset index

  //Loop on datasets
  while(!infile.eof())
    {
      dtindex = nextdtindex;

      //Read dataset name
      getline(infile, name); lines = readLineFiles(infiles);

      //Initialise new dataset
      Dataset dtset(dtindex, name);

      getline(infile, line); lines = readLineFiles(infiles);
      //Read plot options
      while(!infile.eof())
        {
          if (line.find("Plot") == string::npos)
            break;

          //Read subplot index
          TString temp(line.c_str());
          TObjArray* array = temp.Tokenize("@");
	  temp.Form(((TObjString*)  array->At(0))->GetString().Data());
	  //snprintf(temp,sizeof(temp),"%s",((TObjString*)  array->At(0))->GetString().Data());
	  delete array;
          temp.ReplaceAll("Plot","");
          int iplot = temp.Atoi();
          //End of reading subplotindex

          dtset.subplots[iplot] = Subplot(line);
          getline(infile, line); lines = readLineFiles(infiles);
        }

      //Read columns tags
      string col;
      int i = 0;
      istringstream iss(line);

      while (iss >> col)
        {
          //      coltag[i] = col;
          //      i++;
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

      getline(infile, line);  lines = readLineFiles(infiles);
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

          //Plain 90cl -> 68cl scaling
          if (outdirs[label].Scale68())
            {
              fline["therr+"] = fline["therr+"]/1.645;
              fline["therr-"] = fline["therr-"]/1.645;
            }

          //Hessian profile the theory prediction
          if (outdirs[label].IsProfiled())
            {
              double cor, eplus, eminus;
              getTheoryShift(pdfmap[label].pdfshifts, pdfmap[label].cor_matrix, err, lines, outdirs[label].Scale68(), cor, eplus, eminus);
              fline["thorig"] += cor;
              fline["therr+"] = eplus;
              fline["therr-"] = eminus;
            }

          //Bayesian reweight the theory prediction
          if (outdirs[label].IsReweighted())
            {
              double value, error;
              getTheoryReweight(pdfmap[label].mcw, lines, value, error);
              fline["thorig"] = value;
              fline["therr+"] = error;
              fline["therr-"] = error;
            }

          //Cumulative theory predictions for MC replica runs
          if (outdirs[label].IsMCreplica())
            {
              double value, error;
              vector <double> w;
              getTheoryReweight(w, lines, value, error);
              fline["thorig"] = value;
              fline["therr+"] = error;
              fline["therr-"] = error;
              fline["thmod"] = value;
              fline["pull"] = 0.;
            }

          //Add point to subplot
          dtset.subplots[iplot].AddPoint(fline);
          getline(infile, line); lines = readLineFiles(infiles);

        }//End loop on data points
      for (map <int, Subplot>::iterator sit = dtset.subplots.begin(); sit != dtset.subplots.end(); sit++)
        (*sit).second.Init(label, dtindex, (*sit).first);

      datamap[dtindex] = dtset;
    }//End loop on datasets

  //Close files to avoid memory leak
  for (int i = 0; i < infiles.size(); i++)
    if (infiles[i]->is_open())
      infiles[i]->close();
}
