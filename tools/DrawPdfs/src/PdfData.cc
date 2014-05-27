#include "PdfData.h"

#include "Outdir.h"
#include "CommandParser.h"
#include "Par.h"

#include <TMath.h>

#include <fstream>
#include <stdlib.h>
#include <math.h>

//pdf type
pdftype pdfts[] = {g, uv, dv, ubar, dbar, s, U, D, Sea, Ubar, Dbar, c, b, dbarubar, sdbar};
//pdf labels
string pdflab[] = {"g", "uv", "dv", "#bar{u}", "#bar{d}", "s", "U", "D", "S", "#bar{U}", "#bar{D}", "c", "b", "#bar{d}-#bar{u}", "(s - #bar{d})/(s + #bar{d})"};
//pdf filenames
string pdffil[] = {"g", "uv", "dv", "ubar", "dbar", "s", "U", "D", "sea", "UBar", "DBar", "c", "b", "dbar-ubar", "sdbar"};

vector <pdftype> pdfs(pdfts, pdfts + sizeof(pdfts) / sizeof(pdftype));
vector <string> pdflabels(pdflab, pdflab + sizeof(pdflab) / sizeof(string));
vector <string> pdffiles(pdffil, pdffil + sizeof(pdffil) / sizeof(string));


Pdf::Pdf(string filename) : Q2value(0), NxValues(0), NPdfs(0), Xmin(0), Xmax(0)
{
  //Open PDF file
  ifstream infile(filename.c_str());
  if (!infile.is_open()) //PDF file not found
    return;

  // Read the file header
  infile >> Q2value >> NxValues >> NPdfs >> Xmin >> Xmax;

  // Read the column names
  vector <pdftype> PdfTypes;
  for (int i = 0; i <= NPdfs; i++)
    {
      string var;
      infile >> var;
      pdftype ipdf;
      if (var ==  "g") ipdf = g;
      else if (var ==  "U") ipdf = U;
      else if (var ==  "D") ipdf = D;
      else if (var ==  "u_val") ipdf = uv;
      else if (var ==  "d_val") ipdf = dv;
      else if (var ==  "Ubar") ipdf = Ubar;
      else if (var ==  "Dbar") ipdf = Dbar;
      else if (var ==  "sea") ipdf = Sea;
      else if (var ==  "u_sea") ipdf = ubar;
      else if (var ==  "d_sea") ipdf = dbar;
      else if (var ==  "str") ipdf = s;
      else if (var ==  "chm") ipdf = c;
      else if (var ==  "bot") ipdf = b;
      else if (var != "x")
	{
	  cout << "Error: pdf " << var << " not recognised" << endl;
	  exit(1);
	}
      PdfTypes.push_back(ipdf);
    }

  // Read the table
  double val;
  for (int ix = 0; ix < NxValues; ix++)
    {
      //read x value
      infile >> val;
      xpoints.push_back(val);
      for (int j = 1; j <= NPdfs; j++)
	{
	  infile >> val;
	  tablemap[PdfTypes[j]].push_back(val);
	}
    }

  //custom pdf types
  PdfTypes.push_back(dbarubar);  NPdfs++;
  for (int ix = 0; ix < NxValues; ix++)
    tablemap[dbarubar].push_back(tablemap[dbar][ix] - tablemap[ubar][ix]);

  PdfTypes.push_back(sdbar);  NPdfs++;
  for (int ix = 0; ix < NxValues; ix++)
    if ((tablemap[dbar][ix] + tablemap[s][ix]) != 0)
      tablemap[sdbar].push_back((tablemap[s][ix] - tablemap[dbar][ix])/(tablemap[s][ix] + tablemap[dbar][ix]));
    else
      tablemap[sdbar].push_back(0);
}

TGraphAsymmErrors* Pdf::GetPdf(pdftype ipdf)
{
  if (GetNx() == 0)
    return 0;
  TGraphAsymmErrors* temp = new TGraphAsymmErrors(GetNx());
  for (int i = 0; i < GetNx();  i++)
    temp->SetPoint(i, xpoints[i],tablemap[ipdf][i]);
  if (tablemapup[ipdf].size() > 0 && tablemapdn[ipdf].size() > 0)
    for (int i = 0; i < GetNx();  i++)
      {
	temp->SetPointEYhigh(i, tablemapup[ipdf][i]);
	temp->SetPointEYlow(i, tablemapdn[ipdf][i]);
      }
  return temp;
}

PdfData::PdfData(string dirname, string label)
{
  if (outdirs[label].IsMCreplica())
    {
      err = MC;
      //Loop on Q2 values, read PDF tables
      int iq2 = 0;
      while (true)
	{
	  iq2++;
	  //Load PDF MC replica sets
	  char fname[300];
	  int iband = 0;
	  Pdf temppdf;
	  float q2 = 0;
	  for (vector <string>::iterator it = outdirs[label].dirlist.begin(); it != outdirs[label].dirlist.end(); it++)
	    {
	      iband++;

	      sprintf(fname, "%s/pdfs_q2val_%02d.txt", (*it).c_str(), iq2);
	      Pdf temperr(fname);
	      temppdf = temperr;
	      q2 = temperr.GetQ2();
	      if (temperr.GetQ2() == 0)
		break;

	      if (temperr.GetNx() > 0)
		Errors[temperr.GetQ2()].push_back(temperr);
	    }
	  if (q2 == 0)
	    break;
	  //build central PDF afterwards
	  Central[temppdf.GetQ2()] = temppdf;
	}
    }
  else
    {
      //check if errors are symmetric, asymmetric hessian, Monte Carlo
      char filename[300];
      //asymmetric hessian
      sprintf (filename, "%s/pdfs_q2val_s%02dm_%02d.txt", dirname.c_str(), 1, 1);
      ifstream asfile(filename);
      if (asfile.is_open())
	{
	  err = AsymHess;
	  asfile.close();
	}

      //symmetric hessian errors
      sprintf (filename, "%s/pdfs_q2val_s%02ds_%02d.txt", dirname.c_str(), 1, 1);
      ifstream sfile(filename);
      if (sfile.is_open())
	{
	  err = SymHess;
	  sfile.close();
	}

      //MC errors
      sprintf (filename, "%s/pdfs_q2val_mc%03ds_%02d.txt",dirname.c_str(), 1, 1);
      ifstream mcfile(filename);
      if (mcfile.is_open())
	{
	  err = MC;
	  mcfile.close();
	}

      //Loop on Q2 values, read PDF tables
      int iq2 = 0;
      while (true)
	{
	  iq2++;

	  //Load central PDF
	  char fname[300];
	  sprintf(fname, "%s/pdfs_q2val_%02d.txt", dirname.c_str(), iq2);
	  Pdf temppdf(fname);
	  if (temppdf.GetQ2() == 0)
	    break; //if central file for this q2 index not found, break q2 loop

	  Central[temppdf.GetQ2()] = temppdf;

	  //Get Pdf errors if requested
	  if (!opts.dobands)
	    continue;
  
	  //Load PDF error sets
	  int iband = 0;
	  if (err == MC)
	    while (true)
		{
		  iband++;
		  sprintf (fname, "%s/pdfs_q2val_mc%03ds_%02d.txt", dirname.c_str(), iband, iq2);
		  Pdf temperr(fname);
		  if (temperr.GetQ2() == 0)
		    break;
		  if (temperr.GetNx() > 0)
		    Errors[temperr.GetQ2()].push_back(temperr);
		}
	  else if (err == SymHess)
	    while (true)
	      {
		iband++;
		sprintf (fname, "%s/pdfs_q2val_s%02ds_%02d.txt", dirname.c_str(), iband, iq2);
		Pdf temperr(fname);
		if (temperr.GetQ2() == 0)
		  break;
		if (temperr.GetNx() > 0)
		  Errors[temperr.GetQ2()].push_back(temperr);
	      }
	  else if (err == AsymHess)
	    while (true)
	      {
		iband++;
		//positive variation
		sprintf (fname, "%s/pdfs_q2val_s%02dm_%02d.txt", dirname.c_str(), iband, iq2);
		Pdf temperrplus(fname);
		if (temperrplus.GetQ2() == 0)
		  break;
		if (temperrplus.GetNx() > 0)
		  Errors[temperrplus.GetQ2()].push_back(temperrplus);

		//negative variation
		sprintf(fname, "%s/pdfs_q2val_s%02dp_%02d.txt", dirname.c_str(), iband, iq2);
		Pdf temperrminus(fname);
		if ((temperrminus.GetQ2() != temperrplus.GetQ2()) || temperrminus.GetNx() == 0)
		  {
		    cout << "Error, Asymmetric hessian PDF uncertainties, positive variation found, but cannot find down variation: " << fname << endl;
		    break;
		  }
		Errors[temperrminus.GetQ2()].push_back(temperrminus);
	      }
	}
    }

  //Remake central PDF
  if (err == MC)
    for (map<float, Pdf>::iterator pdfit = Central.begin(); pdfit != Central.end(); pdfit++) //Loop on q2 values
      {
	float q2 = pdfit->first;
	for (vector <pdftype>::iterator pit = pdfs.begin(); pit != pdfs.end(); pit++) //loop on pdf types
	  for (int ix = 0; ix < Central[q2].GetNx(); ix++) //Loop on x points
	    {
	      vector <double> xi;
	      for (vector <Pdf>::iterator eit = Errors[q2].begin(); eit != Errors[q2].end(); eit++)
		xi.push_back((*eit).GetTable(*pit)[ix]);
	      double val;
	      if (outdirs[label].IsMedian())
		val = Median(xi);
	      else
		val = TMath::Mean(xi.begin(), xi.end());
	      pdfit->second.SetPoint(*pit, ix, val);
	    }
      }

  //Compute PDF uncertainty bands
  if (!opts.dobands)
    return;

  //Loop on q2 values
  for (map<float, Pdf>::iterator pdfit = Central.begin(); pdfit != Central.end(); pdfit++)
    {
      float q2 = pdfit->first;
      Pdf Cent = pdfit->second;
      Up[q2] = Cent;
      Down[q2] = Cent;

      //loop on pdf types
      for (vector <pdftype>::iterator pit = pdfs.begin(); pit != pdfs.end(); pit++)
	{
	  //Loop on x points
	  for (int ix = 0; ix < Cent.GetNx(); ix++)
	    {
	      double val = Cent.GetTable(*pit)[ix];
	      double eminus = 0;
	      double eplus = 0;
	    
	      //MC replica errors
	      if (err == MC)
		{
		  vector <double> xi;
		  for (vector <Pdf>::iterator eit = Errors[q2].begin(); eit != Errors[q2].end(); eit++)
		    xi.push_back((*eit).GetTable(*pit)[ix]);
		  if (outdirs[label].IsMedian())
		    val = Median(xi);
		  else
		    val = TMath::Mean(xi.begin(), xi.end());

		  if (outdirs[label].Is68cl())
		    {
		      if (outdirs[label].IsAsym())
			deltaasym(xi, val, eplus, eminus, cl(1));
		      else
			eplus = eminus = delta(xi, val, cl(1));
		    }
		  else if (outdirs[label].Is90cl())
		    {
		      if (outdirs[label].IsAsym())
			deltaasym(xi, val, eplus, eminus, 0.90);
		      else
			eplus = eminus = delta(xi, val, 0.90);
		    }
		  else
		    eminus = eplus = TMath::RMS(xi.begin(), xi.end());
		}
	      else if (err == AsymHess)
		{
		  //Loop on error sets
		  for (vector <Pdf>::iterator eit = Errors[q2].begin(); eit != Errors[q2].end(); eit+=2)
		    {
		      double vm = (*eit).GetTable(*pit)[ix];
		      double vp = (*(eit+1)).GetTable(*pit)[ix];
		      if (!outdirs[label].IsAsym()) //symmetrise errors
			{
			  double err = 0.5*(vp-vm);
			  eminus += err*err;
			  eplus  += err*err;
			}
		      else //asymmetric errors
			{
			  // down variation:
			  double d1 = val - vm;
			  double d2 = val - vp;
			  double ed = ( d2>d1) ? d2 : d1;
			
			  if (ed<0) { ed = 0;}
			
			  // up variation
			  d1 = -d1;
			  d2 = -d2;
			  double ep = (d2>d1) ? d2 : d1;    
			  if (ep<0) {ep = 0;}
			  eminus += ed*ed;
			  eplus  += ep*ep;
			}
		    }
		  eminus = sqrt(eminus);
		  eplus = sqrt(eplus);
		}
	      else if (err == SymHess)
		{
		  double maxp = 0;
		  double maxm = 0;
		  //Loop on error sets
		  for (vector <Pdf>::iterator eit = Errors[q2].begin(); eit != Errors[q2].end(); eit+=2)
		    {
		      double v = (*eit).GetTable(*pit)[ix];
		      double err = v - val;
		      eminus += err*err;
		      eplus  += err*err;
		    }
		  eminus = sqrt(eminus);
		  eplus = sqrt(eplus);
		}
	      Up[q2].SetPoint(*pit, ix, val+eplus);
	      Down[q2].SetPoint(*pit, ix, val-eminus);
	      pdfit->second.SetPoint(*pit, ix, val);
	      pdfit->second.SetErrUp(*pit, ix, eplus);
	      pdfit->second.SetErrDn(*pit, ix, eminus);
	    } //End loop on x points
	}//End loop on pdf types
    }//End loop on q2 values
}
