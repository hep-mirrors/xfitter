#include "Outdir.h"

#include <dirent.h>
#include <stdlib.h>
#include <algorithm>

map <string, PdfData> pdfmap;
map <string, Data> datamap;
map <string, Chi2> chi2map;
map <string, Par> parmap;

//Constructor, load all the directory data
Outdir::Outdir(string dir, string lab) : dirname(dir), label(lab)
{
  //check directory exists
  DIR *pDir;
  pDir = opendir (dirname.c_str());
  if (pDir != NULL)
    closedir (pDir);
  else
    {
      cout << "Error: directory " << dirname << " not found" << endl;
      exit(1);
    }

  //Load PDF data
  PdfData pdf(dirname, label);
  pdfmap[label] = pdf;

  //Load datasets
  Data data(dirname, label);
  datamap[label] = data;

  //Load chi2
  Chi2 chi2(dirname, label);
  chi2map[label] = chi2;

  //Load parameters
  Par par(dirname, label);
  parmap[label] = par;
}

vector <float> q2list()
{
  vector <float> q2;
  for (map <string, PdfData>::iterator pit = pdfmap.begin(); pit != pdfmap.end(); pit++)
    for (map <float, Pdf>::iterator qit = pit->second.Central.begin(); qit != pit->second.Central.end(); qit++)
      q2.push_back(qit->first);

  sort (q2.begin(), q2.end());
  vector<float>::iterator it = unique (q2.begin(), q2.end());
  q2.resize(distance(q2.begin(), it));

  return q2;
}

vector <int> datalist()
{
  vector <int> dl;
  for (map <string, Data>::iterator dit = datamap.begin(); dit != datamap.end(); dit++) //loop on directories
    for (map <int, Dataset>::iterator iit = dit->second.datamap.begin(); iit != dit->second.datamap.end(); iit++) //loop on dataset index
      for (map <int, Subplot>::iterator sit = iit->second.subplots.begin(); sit != iit->second.subplots.end(); sit++) //loop on subplots
	dl.push_back(sit->first + iit->first * 100);

  sort (dl.begin(), dl.end());
  vector<int>::iterator it = unique (dl.begin(), dl.end());
  dl.resize(distance(dl.begin(), it));

  return dl;
}

vector <int> chi2list()
{
  vector <int> chi2;
  for (map <string, Chi2>::iterator cit = chi2map.begin(); cit != chi2map.end(); cit++) //loop on directories
    for (map <int, chi2type>::iterator it = cit->second.chi2list.begin(); it != cit->second.chi2list.end(); it++) //loop on dataset index
      chi2.push_back(it->first);

  sort (chi2.begin(), chi2.end());
  vector<int>::iterator it = unique (chi2.begin(), chi2.end());
  chi2.resize(distance(chi2.begin(), it));

  return chi2;
}

vector <int> parlist()
{
  vector <int> par;
  for (map <string, Par>::iterator dit = parmap.begin(); dit != parmap.end(); dit++) //loop on directories
    for (map <int, partype>::iterator it = dit->second.parlist.begin(); it != dit->second.parlist.end(); it++) //loop on parameters
      par.push_back(it->first);

  sort (par.begin(), par.end());
  vector<int>::iterator it = unique (par.begin(), par.end());
  par.resize(distance(par.begin(), it));

  return par;
}

int finddataindex(string name)
{
  int index = -1;

  for (map <string, Data>::iterator dit = datamap.begin(); dit != datamap.end() && index == -1; dit++) //loop on directories
    for (map <int, Dataset>::iterator iit = dit->second.datamap.begin(); iit != dit->second.datamap.end() && index == -1; iit++) //loop on dataset index
      if (iit->second.GetName() == name)
	index = iit->second.GetIndex();

  return index;
}
string finddataname(int index)
{
  string name = "";

  for (map <string, Data>::iterator dit = datamap.begin(); dit != datamap.end() && name == ""; dit++) //loop on directories
    for (map <int, Dataset>::iterator iit = dit->second.datamap.begin(); iit != dit->second.datamap.end() && name == ""; iit++) //loop on dataset index
      if (iit->second.GetIndex() == index)
	name = iit->second.GetName();

  return name;
}
string findparname(int index)
{
  string name = "";

  for (map <string, Par>::iterator dit = parmap.begin(); dit != parmap.end() && name == ""; dit++) //loop on directories
    for (map <int, partype>::iterator iit = dit->second.parlist.begin(); iit != dit->second.parlist.end() && name == ""; iit++) //loop on parameters index
      if (iit->first == index)
	name = iit->second.name;

  return name;
}
