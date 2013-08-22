#include <iostream>
#include <string>
#include <stdlib.h>
#include <TError.h>

#include <Output.h>
#include <DataPainter.h>

using std::cout;
using std::cerr;
using std::endl;

using namespace std;

int main(int argc, char **argv) 
{

  string OutputPath("output/");
  
  gErrorIgnoreLevel=1001;

  int ibase1 = 0;
  int ibase2 = 0;
  int ibase3 = 0;

  if ( argc == 1 ) 
    {
      printf("program usage:\n DrawDataPulls dir1 [dir2] [dir3]\n");
      return -1;
    }
  
  if ( argc > 1 ) 
    for (int iar=1; iar<argc; iar++)
      if (ibase1 == 0) 
	ibase1 = iar;
      else if (ibase2 == 0)
	ibase2 = iar;
      else if (ibase3 == 0)
	ibase3 = iar;
      else
	{
	  printf("program usage:\n DrawDataPulls dir1 [dir2] [dir3]\n");
	  return -1;
	}

  //initialize datasets
  if (ibase2 ==0){
    OutputPath = argv[ibase1];
    cout << "plots are stored in: " << (argv[ibase1]) << endl;
  }
  else{
    OutputPath = "datapulls/";
    cout << "plots are stored in: datapulls/" << endl;
  }

  vector <Output*> info_output;
  info_output.push_back(new Output(argv[ibase1]));
  if (ibase2 != 0)
    info_output.push_back(new Output(argv[ibase2]));
  if (ibase3 != 0)
    info_output.push_back(new Output(argv[ibase3]));

  map <int, vector<dataseth> > datamap;
  for (unsigned int o = 0; o < info_output.size(); o++)
    {
      info_output[o]->Prepare(false);
      for (unsigned int d = 0; d < info_output[o]->GetNsets(); d++)
	{
	  if (info_output[o]->GetSet(d)->GetNSubPlots() == 1)
	    {
	      int id = info_output[o]->GetSet(d)->GetSetId();
	      //check bins sanity
	      vector <float> b1 = info_output[o]->GetSet(d)->getbins1(0);
	      vector <float> b2 = info_output[o]->GetSet(d)->getbins2(0);
	      vector<float>::iterator it1 = b1.begin();
	      vector<float>::iterator it2 = b2.begin();
	      bool skip = false;
	      for (; (it1+1) != b1.end(); it1++, it2++)
		if (*(it1+1) != *it2 || *it1 >= *(it1+1))
		  skip = true;
	      if (skip)
		{
		    cout << "bin inconsistency for dataset: " << info_output[o]->GetSet(d)->GetName() << endl;
		    cout << "Cannot plot, skipping" << endl;
		    continue;
		}

	      dataseth dt = dataseth(info_output[o]->GetSet(d)->GetName(),
				     info_output[o]->GetName()->Data(),
				     info_output[o]->GetSet(d)->getbins1(0),
				     info_output[o]->GetSet(d)->getbins2(0),
				     info_output[o]->GetSet(d)->getdata(0),
				     info_output[o]->GetSet(d)->getuncor(0),
				     info_output[o]->GetSet(d)->gettoterr(0),
				     info_output[o]->GetSet(d)->gettheory(0),
				     info_output[o]->GetSet(d)->gettheoryshifted(0),
				     info_output[o]->GetSet(d)->getpulls(0));
	      datamap[id].push_back(dt);
	    }
	}
    }

  vector <TCanvas*> canvaslist;
  for (map <int, vector <dataseth> >::iterator it = datamap.begin(); it != datamap.end(); it++)
    canvaslist.push_back(DataPainter((*it).first, (*it).second));

  if (OutputPath.rfind("/") != OutputPath.size() - 1)
    OutputPath.append("/");
  system(((string)"mkdir -p " + OutputPath).c_str());
  for (vector <TCanvas*>::iterator it = canvaslist.begin(); it != canvaslist.end(); it++)
    {
      string out = OutputPath + (*it)->GetName() + ".eps";
      (*it)->Print(out.c_str());
    }

  return 0;
}
