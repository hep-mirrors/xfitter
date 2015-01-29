#include "PdfData.h"

#include "Outdir.h"
#include "CommandParser.h"
#include "pdferrors.h"
#include "Par.h"

#include <TMath.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "FileOpener.h"

//pdf type
pdftype pdfts[] = {uv, dv, g, Sea, ubar, dbar, s, Rs, c, b, dbarubar, uvdv, U, D, Ubar, Dbar};
//pdf labels
string pdflab[] = {"u_{V}", "d_{V}", "g", "#Sigma", "#bar{u}", "#bar{d}", "s", "(s+#bar{s})/(#bar{u}+#bar{d})", "c", "b", "#bar{d}-#bar{u}", "u_{V}-d_{V}", "U", "D", "#bar{U}", "#bar{D}"};
//pdf filenames
string pdffil[] = {"uv", "dv", "g", "Sea", "ubar", "dbar", "s", "Rs", "c", "b", "dbar-ubar", "uvdv", "U", "D", "UBar", "DBar"};

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

  PdfTypes.push_back(Rs);  NPdfs++;
  for (int ix = 0; ix < NxValues; ix++)
    if (tablemap[dbar][ix] != 0)
      tablemap[Rs].push_back((2*tablemap[s][ix])/(tablemap[ubar][ix]+tablemap[dbar][ix]));
    else
      tablemap[Rs].push_back(0);

  PdfTypes.push_back(uvdv);  NPdfs++;
  for (int ix = 0; ix < NxValues; ix++)
    tablemap[uvdv].push_back(tablemap[uv][ix] - tablemap[dv][ix]);
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

PdfData::PdfData(string dirname, string label) : model(false), par(false)
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
          Pdf temppdf;
          float q2 = 0;
          for (vector <string>::iterator it = outdirs[label].dirlist.begin(); it != outdirs[label].dirlist.end(); it++)
            {
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
	  if (!opts.dobands && !outdirs[label].IsProfiled() && !outdirs[label].IsRotated())
	    continue;
  
          //Load PDF error sets
          int iband = 1;
          if (err == MC)
            while (true)
                {
                  sprintf (fname, "%s/pdfs_q2val_mc%03ds_%02d.txt", dirname.c_str(), iband, iq2);
                  Pdf temperr(fname);
                  if (temperr.GetQ2() == 0)
                    break;
                  if (temperr.GetNx() > 0)
                    Errors[temperr.GetQ2()].push_back(temperr);
                  iband++;
                }
          else if (err == SymHess)
            while (true)
              {
                sprintf (fname, "%s/pdfs_q2val_s%02ds_%02d.txt", dirname.c_str(), iband, iq2);
                Pdf temperr(fname);
                if (temperr.GetQ2() == 0)
                  break;
                if (temperr.GetNx() > 0)
                  Errors[temperr.GetQ2()].push_back(temperr);
                iband++;
              }
          else if (err == AsymHess)
            while (true)
              {
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
                iband++;
              }

          //check for model errors
          sprintf (filename, "%s/pdfs_q2val_m%02dm_%02d.txt",dirname.c_str(), iband+1, 1);
          ifstream modelfile(filename);
          if (modelfile.is_open())
            {
              model = true;
              modelfile.close();
            }
          if (model)
            while (true)
              {
                //positive variation
                sprintf (fname, "%s/pdfs_q2val_m%02dm_%02d.txt", dirname.c_str(), iband, iq2);
                Pdf temperrplus(fname);
                if (temperrplus.GetQ2() == 0)
                  break;
                if (temperrplus.GetNx() > 0)
                  ModelErrors[temperrplus.GetQ2()].push_back(temperrplus);

                //negative variation
                sprintf(fname, "%s/pdfs_q2val_m%02dp_%02d.txt", dirname.c_str(), iband, iq2);
                Pdf temperrminus(fname);
                if ((temperrminus.GetQ2() != temperrplus.GetQ2()) || temperrminus.GetNx() == 0)
                  {
                    cout << "Error, Model PDF uncertainties, positive variation found, but cannot find down variation: " << fname << endl;
                    break;
                  }
                ModelErrors[temperrminus.GetQ2()].push_back(temperrminus);
                iband++;
              }

          //check for parametrisation errors
          sprintf (filename, "%s/pdfs_q2val_p%02ds_%02d.txt",dirname.c_str(), iband+1, 1);
          ifstream parfile(filename);
          if (parfile.is_open())
            {
              par = true;
              parfile.close();
            }
          if (par)
            while (true)
              {
                sprintf (fname, "%s/pdfs_q2val_p%02ds_%02d.txt", dirname.c_str(), iband, iq2);
                Pdf temperr(fname);
                if (temperr.GetQ2() == 0)
                  break;
                if (temperr.GetNx() > 0)
                  ParErrors[temperr.GetQ2()].push_back(temperr);
                iband++;
              }
        }
    }

  //Remake central PDF
  if (err == MC && opts.dobands)
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
                val = median(xi);
              else
                val = mean(xi);
              pdfit->second.SetPoint(*pit, ix, val);
            }
      }

  //Compute PDF uncertainty bands
  if (!opts.dobands && !outdirs[label].IsProfiled() &&!outdirs[label].IsRotated())
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
                    val = median(xi);
                  else
                    val = mean(xi);

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
                    eminus = eplus = rms(xi);
                }
              else if (err == AsymHess)
                {
                  vector <double> xi;
                  xi.push_back(val);
                  for (vector <Pdf>::iterator eit = Errors[q2].begin(); eit != Errors[q2].end(); eit++)
                    xi.push_back((*eit).GetTable(*pit)[ix]);

		  if (!outdirs[label].IsAsym()) //symmetrise errors
		    eplus = eminus = ahessdelta(xi);
		  else //asymmetric errors
		    ahessdeltaasym(xi, eplus, eminus);
		  if (outdirs[label].Scale68())
		    {
		      eplus = eplus/1.645;
		      eminus = eminus/1.645;
		    }
		}
              else if (err == SymHess)
                {
                  vector <double> xi;
                  xi.push_back(val);
                  for (vector <Pdf>::iterator eit = Errors[q2].begin(); eit != Errors[q2].end(); eit++)
		    xi.push_back((*eit).GetTable(*pit)[ix]);

		  eplus = eminus = shessdelta(xi);
		  if (outdirs[label].Scale68())
		    {
		      eplus = eplus/1.645;
		      eminus = eminus/1.645;
		    }
		}

              //Add model and parametrisation uncertainties
              if (model)
                {
                  vector <double> xi;
                  xi.push_back(val);

                  for (vector <Pdf>::iterator eit = ModelErrors[q2].begin(); eit != ModelErrors[q2].end(); eit++)
                    xi.push_back((*eit).GetTable(*pit)[ix]);
                  
                  double modeplus, modeminus;
                  if (!outdirs[label].IsAsym()) //symmetrise errors
                    modeplus = modeminus = ahessdelta(xi);
                  else //asymmetric errors
                    ahessdeltaasym(xi, modeplus, modeminus);
                  eplus = sqrt(pow(eplus, 2) + pow(modeplus, 2));
                  eminus = sqrt(pow(eminus, 2) + pow(modeminus, 2));
                }
              if (par)
                {
                  vector <double> xi;
                  xi.push_back(val);

                  for (vector <Pdf>::iterator eit = ParErrors[q2].begin(); eit != ParErrors[q2].end(); eit++)
                    xi.push_back((*eit).GetTable(*pit)[ix]);

                  double pareplus, pareminus;
                  deltaenvelope(xi, pareplus, pareminus);
                  if (!outdirs[label].IsAsym()) //symmetrise errors
                    {
                      double parerr = (pareminus + pareplus)/2;
                      pareplus = pareminus = parerr;
                    }
                  eplus = sqrt(pow(eplus, 2) + pow(pareplus, 2));
                  eminus = sqrt(pow(eminus, 2) + pow(pareminus, 2));
                }

              Up[q2].SetPoint(*pit, ix, val+eplus);
              Down[q2].SetPoint(*pit, ix, val-eminus);
              pdfit->second.SetPoint(*pit, ix, val);
              pdfit->second.SetErrUp(*pit, ix, eplus);
              pdfit->second.SetErrDn(*pit, ix, eminus);
            } //End loop on x points
        }//End loop on pdf types
    }//End loop on q2 values

  if (outdirs[label].IsProfiled() || opts.profiled)
    profile(dirname, label);
  if (outdirs[label].IsRotated() )
    pdfRotate(dirname, label);
  
}


void PdfData::pdfRotate(string dirname, string label)
{ 
  // Extra rotations from rot.dat file
  string fname = dirname + "/pdf_rot.dat";
  ifstream f(fname.c_str());
  if ( ! f.good() ) {
    cout << "File " << fname << " is empty (or io error)" << endl;
    exit(1);
  }

  vector< vector<double> > rotation;
  string line;
  getline (f,line);
  std::cout<<line <<std::endl;
  getline (f, line);
  istringstream iss(line);
  int N;
  iss >> N;
  int idx1 = 0;
  while ( getline (f,line) ) 
    {
	vector <double> aline;
	istringstream iss(line);
	int idx2;
	iss >> idx2;
	for ( int i = 0; i<N; i++) {
	  double val;
	  iss >> val;
	  aline.push_back(val);
	}
	rotation.push_back(aline);
    }
  f.close();

  int iRotation = outdirs[label].rSet()-1;
  //  for (int id=0; id<N; id++) {
  //  std::cout << iRotation << " " << rotation[iRotation][id]<<"\n";
  // }

  for ( map<float, Pdf>::iterator pdfit = Central.begin(); pdfit != Central.end(); pdfit++) {
    float q2 = pdfit->first;
    Pdf Cent = pdfit->second;
    

    // loop over pdf types
    for (vector <pdftype>::iterator pit = pdfs.begin(); pit != pdfs.end(); pit++) {
      //Loop on x points
      for (int ix = 0; ix < Cent.GetNx(); ix++)
	{
	  double val = Cent.GetTable(*pit)[ix];
	  double corsum = 0;
	  double eminus = 0; // also  errors
	  double eplus = 0;  
	
	  // For now CT10 only:
	  for ( int id=0; id<N; id++) {
	    Pdf Up = Errors[q2].at(2*(id));
	    Pdf Dn = Errors[q2].at(2*(id)+1);
	    double plus  = Up.GetTable(*pit)[ix] - val;
	    double minus = Dn.GetTable(*pit)[ix] - val;


	    corsum += 0.5*(plus-minus)*rotation[iRotation][id];

	    //	    corsum += 0.5*(plus-minus)*rotation[iRotation][id];
	  }
	  
	  Cent.SetPoint(*pit, ix, val+corsum);
	  Cent.SetErrUp(*pit, ix, eplus);
	  Cent.SetErrDn(*pit, ix, eminus);
	  

	  Up[q2].SetPoint(*pit, ix, val+corsum+eplus);
	  Down[q2].SetPoint(*pit, ix, val+corsum-eminus);
	}
    }
    pdfit->second = Cent;
  }
  
}

void PdfData::profile(string dirname, string label)
{

  // Extract PDF shifts from Results.txt:

  // string fname = dirname + "/Results.txt";
  // ifstream f(fname.c_str());
  // if ( ! f.good() ) {
    // cout << "File " << fname << " is empty (or io error)" << endl;
    // return;
  // }
  
  InFileOpener_t fo;
  fo.Add(dirname + "/Results.txt");
  fo.Add(dirname + "/Results_0.txt");
  // fo.Clear();
  if(fo.Open()) return;
  ifstream &f = fo.GetStream();

  string line;
  string buffer = "";
  while (buffer != "Name")
    {
      getline(f, line);
      istringstream iss(line);
      iss >> buffer; 
    }
  string systlabel, dummy;
  float systindex, value, error;
  int counter = 0;
  while (getline(f, line))
    {
      istringstream iss(line);
      iss >> systindex >> systlabel  >> value  >> dummy  >> error;       
      if ( systlabel.find("PDF_nuisance_param") == 0 ) {
        ++counter;
        pdfshift shift;
        shift.val = value;
        shift.err = error;
        shift.id  = counter;
        pdfshifts.push_back(shift);
      }
    }
  f.close();

  // Read also PDF correlation matrix:
  string ffname = dirname + "/pdf_vector_cor.dat";
  ifstream ff(ffname.c_str());
  vector< vector<double> > cor_matrix;
  if ( ! ff.good() ) {
    cout << "File " << ffname << " is empty (or io error). Use diagonal approximation for the PDF nuisance parameters." << endl;
  }
  else {
    getline (ff,line);
    istringstream iss(line);
    int N;
    iss >> N;
    int idx1 = 0;
    while ( getline (ff,line) ) 
      {
        vector <double> aline;
        istringstream iss(line);
        int idx2;
        iss >> idx2;
        for ( int i = 0; i<N; i++) {
          double val;
          iss >> val;
          aline.push_back(val);
        }
        cor_matrix.push_back(aline);
      }
    ff.close();
  }
  // over all Q2 values
  for ( map<float, Pdf>::iterator pdfit = Central.begin(); pdfit != Central.end(); pdfit++) {
    float q2 = pdfit->first;
    Pdf Cent = pdfit->second;
    

    // loop over pdf types
    for (vector <pdftype>::iterator pit = pdfs.begin(); pit != pdfs.end(); pit++) {
      //Loop on x points
      for (int ix = 0; ix < Cent.GetNx(); ix++)
        {
          double val = Cent.GetTable(*pit)[ix];
          double t1 = Up[q2].GetTable(*pit)[ix];
          double t2 = Down[q2].GetTable(*pit)[ix];
          double corsum = 0;
          double eminus = 0; // also  errors
          double eplus = 0;  
          vector <double> xi;
          xi.push_back(val);

          for ( vector<pdfshift>::iterator shift = pdfshifts.begin(); shift != pdfshifts.end(); shift++) {      

            int id = shift->id;
            double valShift = shift->val;
            double errShift = shift->err;


            if (err == AsymHess) {
              Pdf Up = Errors[q2].at(2*(id-1));
              Pdf Dn = Errors[q2].at(2*(id-1)+1);

              double plus  = Up.GetTable(*pit)[ix] - val;
              double minus = Dn.GetTable(*pit)[ix] - val;
              

              double cor = 0.5*(plus - minus)*valShift   - 0.5*(plus+minus)*valShift*valShift;
              
              xi.push_back(plus*errShift+val);
              xi.push_back(minus*errShift+val);

              corsum += cor;
            }
            else if (err = SymHess) {
              Pdf Up = Errors[q2].at(id-1);
              double plus =  Up.GetTable(*pit)[ix] - val;
              double cor = - plus*valShift;              
              xi.push_back(plus*errShift+val);

              corsum += cor;
            }                
          }

          if ( err == AsymHess ) {

            if (!outdirs[label].IsAsym()) //symmetrise errors
              eplus = eminus = ahessdelta(xi, cor_matrix);
            else //asymmetric errors
              ahessdeltaasym(xi, eplus, eminus, cor_matrix);            
          }

          else if (err == SymHess) {
            eplus = eminus = shessdelta(xi, cor_matrix );
          }                
          
          Cent.SetPoint(*pit, ix, val+corsum);
          Cent.SetErrUp(*pit, ix, eplus);
          Cent.SetErrDn(*pit, ix, eminus);


          Up[q2].SetPoint(*pit, ix, val+corsum+eplus);
          Down[q2].SetPoint(*pit, ix, val+corsum-eminus);
        }
    }
    pdfit->second = Cent;
  }
}

