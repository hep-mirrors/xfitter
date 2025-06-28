#include "metamorphCollection.h"

using namespace std;

double vernum = 1.0;
char timeBuffer[20];

//========================================================================

void getTime()
//Get system time
{
  // Set all elements to '\0'
  memset(timeBuffer, '\0', sizeof(timeBuffer));
    
  // Write the value (with a timestamp)
  // Get current time
  time_t now = time(0);
  tm* localTime = localtime(&now);
  
  // Format the time as YYYY-MM-DD HH:MM:SS
  strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", localTime);
} //getTime ----------------------------------------------

template<typename T>
bool isNumber(T x){
  //Returns true if argument is a number. 
  string s;
  stringstream ss; 
  ss << x;
  ss >>s;
  if(s.empty() || isspace(s[0]) || isalpha(s[0])) return false ;
  char * p ;
  strtod(s.c_str(), &p) ;
  return (*p == 0) ;
} //isNumber( T x) ----------------------------------------------

metamorphCollection::metamorphCollection()
{
  MetaVector.reserve(maxMet);
} // metamorphCollection constructor ----------------------------------------------

void metamorphCollection::PushMember()
{
  string flagcheck;
  vector <int> newpositions;
  
  vector<double> Xsvec, fmvec, fpvec, Scmvec; //temporary vectors used to
                                  //initialize values of new control points

  //Read new carrier parameters 
  for (int isc = 0; isc < maxSc; isc++){ // loop copies Sc parameters read
               //from ReadCard() to Scm0[NMeta] array; set DSwitches to 1
               //to update Sc parameters with deltas
    Scm0[NMeta][isc] = stod(strScm[NMeta][isc]);
    DSwitch[NMeta][isc]=1.0;
  }

  //Parse control point parameters
  // For each input string strScm[i], the loop checks if strScm[i] contains
  // a number. If yes, then the input parameters for the corresponding
  // control point will be pushed into a vector to eliminate empty entries
  // caused by non-number CPs. The vector values will be copied by the
  // arrays used to calculate the metamorph objects for each flavor.
  // If strScm is found to be not a number, it will perform additional parsing
  // depending on the type of this special CP
  for (int i = 0; i < k; i++) {
    
    //The string value for the respective CP
    flagcheck = strScm[NMeta][i+maxSc];
    
    if (isNumber(flagcheck)){
      //Numerical CP: fill in the initial value Scm0 and set DSwitch=1
      //to update CP value Sm with de
      Scm0[NMeta][i+maxSc] = stod(flagcheck); // converts string to value
      DSwitch[NMeta][i+maxSc]=1.0; //turn on deltas for this parameter
      l++;
    }
    else if (flagcheck == fixflag){
      //FIX flag: for a FIXED CP, its respective CP value is always zero, 
      //turn off delta changes for this CP
      Scm0[NMeta][i+maxSc] = 0.0; 
      DSwitch[NMeta][i+maxSc]=0.0;
      l++;
    }
    else if (flagcheck == newflag){
      //NEW flag: create a NEW CP from the modulator for other CPs,
      //push its index into the roster of NEW positions;
      // its initial value is computed below; set its DSwitch=1.0 to update 
      DSwitch[NMeta][i+maxSc]=0.0;
      newpositions.push_back(i);
      l++;
    }
    else if (flagcheck == calcflag){
      //CALC flag: do not add the control point, just add the value at
      //at the respective Xs0; do nothing her
    }
    else {
      cerr << "ERROR in parsing the input steering card, metamorph "<<NMeta<<
	",\ncontrol point "<<Xs0[NMeta][i]<<"  "<< flagcheck << endl;
      cerr << "Replace "<<flagcheck << "by a number or a known flag (CALC/FIX/NEW)" << endl;
      exit(1);
    }//else

    if (flagcheck==calcflag || flagcheck==fixflag ||flagcheck==newflag) 
      if (VerbosityLevel > 0) {
	cout << endl << "Special value " << flagcheck << " detected for control point x = " << Xs0[NMeta][i] << endl;
	cout << "Using metamorph::f(x) = Carrier(x) for control point." << endl;
      }
    
    
    //Push known CP values into vectors for creating new control points
    if (isNumber(flagcheck) || flagcheck == fixflag){
      Xsvec.push_back(Xs0[NMeta][i]);
      Scmvec.push_back(Scm0[NMeta][i+maxSc]);
      fpvec.push_back(fp0[NMeta][i]);
      fmvec.push_back(fm0[NMeta][i]);
    } // if (isNumber(flagcheck) == true)
    
  }// for (int i = 0; i < k; i++)

  //If the zeroth or last CPs are fixed, set the xstretching parameters to their x values
  vector <double> vstretch;
  if(strScm[NMeta][maxSc] == fixflag)
    vstretch.push_back(Xs0[NMeta][0]);
  if(strScm[NMeta][maxSc+k-1] == fixflag)
    vstretch.push_back(Xs0[NMeta][k-1]);
		       
  if (! vstretch.empty() && VerbosityLevel > 0){
    cout << "Updating stretching parameters for metamorph "<<NMeta << " to ";
    for (int ipos =0; ipos < (int)vstretch.size(); ipos++)
      cout << vstretch[ipos] << " ";
    cout << endl;
  }
  
  //Compute values for NEW control points
  int Nnew=newpositions.size();
  if (Nnew !=0){
    //First, check that the counts of numerical, fixed, and new control
    //points add up to Nm[Meta] +1

    if ((int)(Xsvec.size() + Nnew) != Nm[NMeta]+1){
      cerr <<"Something is wrong in metamorphCollection::PushMember()"<<endl;
      cerr << "Counts of numerical, fixed, new CPs do not add up to Nm+1"<< endl;
      exit(3);
    }

    int Nmtmp = Xsvec.size()-1;   // degree of the temporary metamorph
                                  //constructed from old and fixed CPs
    double Xstmp[maxctrlpts] = {0}, Scmtmp[maxScm] = {0}, fmtmp[maxctrlpts] = {0}, fptmp[maxctrlpts] = {0};

    for (int i = 0; i < maxSc; i++) //Carrier
      Scmtmp[i] = Scm0[NMeta][i];
    
    for (int i = 0; i < Nmtmp+1; i++){
      Xstmp[i] = Xsvec[i];
      Scmtmp[i+maxSc] = Scmvec[i];
      fmtmp[i] = fmvec[i];
      fptmp[i] = fpvec[i];
    }
    
    metamorph metatmp(Nmtmp, Xstmp, (Scmtmp+maxSc), Scmtmp, xPower[NMeta],vstretch);
    metatmp.SetBoundary(MappingMode[NMeta], fmtmp, fptmp);
    metatmp.UpdateModulator();
    
    //Fill in the values of new points with the values from temporary metamorph
    for (int ipos = 0; ipos < Nnew; ipos++){
      int inew=newpositions[ipos];
      double xtmp = Xs0[NMeta][inew];
      Scm0[NMeta][inew+maxSc] = metatmp.Modulator(xtmp)-1; // sets each Sm parameter to Pi for all control points
    }
    
  }//if (Nnew !=0)
  
  //Construct the final metamorph and push into metaCollection
  for (int i = 0; i < maxSc; i++)
    Scm[NMeta][i] = Scm0[NMeta][i];
  
  for (int i = 0; i < k; i++)
    {
      Xs[NMeta][i] = Xs0[NMeta][i];
      Scm[NMeta][i+maxSc] = Scm0[NMeta][i+maxSc];
      fm[NMeta][i] = fm0[NMeta][i];
      fp[NMeta][i] = fp0[NMeta][i];
    }
  
  MetaVector.emplace_back(Nm[NMeta], Xs[NMeta], Scm[NMeta]+maxSc, Scm[NMeta], xPower[NMeta],vstretch);
  MetaVector[NMeta].SetBoundary(MappingMode[NMeta], fm[NMeta], fp[NMeta]);
  MetaRoster.insert(pair<int, metamorph*>(iflavor[NMeta], &(MetaVector[NMeta])));

  //Set a unique ID for the created metamorph
  MetaVector[NMeta].ID=iflavor[NMeta];
  
  //check the condition number for T
  double condnum=MetaVector[NMeta].GetConditionNumber();
  if (condnum > 10000){
    cerr << "WARNING: a high condition number ="<<condnum
	 << "in metamorph " << MetaVector[NMeta].ID << endl;
    cerr << "Check x spacing of its control points" << endl;
  }
  
} // metamorphCollection::PushMember ----------------------------------------------

void metamorphCollection::ReadCard(const string& inputcard)
{
  ifstream fantosteerin(inputcard);
  if (!fantosteerin.is_open()) {
    cerr << "Unable to open Fantomas input file: " << inputcard << endl;
    exit(3);
  }
  
  // Read and check version of the steering card, then read the timestamp
  string versionheader, timestamp;
  if (!getline(fantosteerin, versionheader) || !getline(fantosteerin, timestamp)) {
    cerr << "Error reading version header from input card." << endl;
    exit(3);
  }
  
  istringstream iss(versionheader);
  string word;
  double vernum_in = 0.0;
  
  while (iss >> word) {
    try {
      vernum_in = stod(word);
    }
    catch (const invalid_argument&) {
      // Skip non-numeric words
    }
  }
  
  if (vernum_in == vernum) {
    if (VerbosityLevel > 0)
      cout << "Reading the Fantomas steering card version " << vernum_in << endl;
  } else {
    cerr << "Fantomas steering card version does not match. Expected " << vernum
	 << ", but found " << vernum_in << "." << endl;
    exit(3);
  }
  
  //Read and parse metamorph blocks until the end of file  
  while (true){
    string line;
    
    // Read first comment line (Flavor header)
    if (!getline(fantosteerin, line))
      break; // end of file cleanly
    if (line.empty()) continue; // skip blank lines
    
    flvcomment[NMeta][1] = line;
    
    // Read second comment line 
    getline(fantosteerin, flvcomment[NMeta][2]);
    
    // Read parameter line
    fantosteerin >> iflavor[NMeta] >> Nm[NMeta] >> MappingMode[NMeta] >> xPower[NMeta];

    if (MappingMode[NMeta] !=0 ){
            cerr << "STOP: a wrong mapping mode "<<MappingMode[NMeta]<<
	      ") for metamorph "<< NMeta << " found in the steering card" << endl;
      cerr << "Only MappingMode==0 currently implemented. Correct it." << endl;
      exit(1);
    }
      
    if (Nm[NMeta] > maxNm){
      cerr << "STOP: the order of Bezier polynomial ("<<Nm[NMeta]<<
	") for metamorph "<< NMeta<<" exceeds\nthe maximal value maxNm="<<maxNm<< endl;
      cerr << "Increase maxNm in metamorphCollection.h" << endl;
      exit(1);
    }

    if (VerbosityLevel > 0) {
      cout << "iflavor = " << iflavor[NMeta] << ", Nm = " << Nm[NMeta]
	   << ", MappingMode = " << MappingMode[NMeta] << ", xPower = " << xPower[NMeta] << endl;
    }
    
    // Read Sc(0), Sc(1), Sc(2)
    for (int i = 0; i < maxSc; ++i) {
      fantosteerin >> strScm[NMeta][i];
      if (VerbosityLevel > 0)
	cout << "Scm[" << i << "] = " << strScm[NMeta][i] << endl;
    }
    
    fantosteerin.ignore(numeric_limits<streamsize>::max(), '\n');
    
    // Read third comment line
    getline(fantosteerin, flvcomment[NMeta][3]);
    
    // Read Nm+1 lines of Xs, Sm, fp, fm
    k = l = m = 0;
    for (int i = 0; i <= Nm[NMeta]; ++i)
      {
	if (!getline(fantosteerin, line)) {
	  cerr << "Error: Unexpected end of file reading control points." << endl;
	  exit(3);
	}
	
	istringstream linestream(line);
	
	linestream >> Xs0[NMeta][i] >> strScm[NMeta][i + maxSc];
	
	if (linestream.eof()) {
	  fp0[NMeta][i] = INFINITY;
	  fm0[NMeta][i] = -1;
	} else {
	  linestream >> fp0[NMeta][i];
	  if (linestream.eof())
	    fm0[NMeta][i] = -1;
	  else
	    linestream >> fm0[NMeta][i];
	}
	++k;
      }
    
    PushMember();
    
    if (l != Nm[NMeta] + 1) {
      cerr << "Error: Nm mismatch for iflavor " << iflavor[NMeta] << endl;
      exit(3);
    }
    
    PositionRoster[iflavor[NMeta]] = NMeta;
    iPts[NMeta] = k;
    ++NMeta;

    if (NMeta > maxMet){
      cerr << "STOP: the number of metamorphs ("<<NMeta<<
	") in the steering card exceeds\nthe maximal allowed number maxMet="<<maxMet<< endl;
      cerr << "Increase maxMet in metamorphCollection.h" << endl;
      exit(1);
    }
     
  }//while(true)
  
  // IMPORTANT: After parameters were read or changed in another way, ALWAYS
  // call UpdateMetamorphs or UpdateModulator[ifl] to update modulator functions
  this->UpdateMetamorphs();
  
  fantosteerin.close();
} // metamorphCollection::ReadCard() ----------------------------------------------

void metamorphCollection::WriteCard(const string& outputcard)
{
  getTime();
    
  ofstream fantosteerout;
  fantosteerout.open(outputcard, ofstream::out);
  
  ios oldState(nullptr); //save the old i/o format such as significant figures
  oldState.copyfmt(fantosteerout);
  
  if (fantosteerout.is_open())
  {
    fantosteerout << std::fixed << std::setprecision(1) <<
      "# Fantomas steering card v. " << vernum << endl;
    fantosteerout << "# " << timeBuffer << endl;
    fantosteerout.copyfmt(oldState); //restore the old i/o format
    
    for (int i = 0; i < NMeta; i++)
    // output fantomas parameters into out card and loop over each input flavor
    {
      fantosteerout << flvcomment[i][1] << endl << flvcomment[i][2] << endl;
      fantosteerout << iflavor[i] << "\t" << Nm[i] << "\t" << MappingMode[i] << "\t" << xPower[i] << "\t";
      // lk24 changed output card to no longer have tab after last Scm value
      for (int j = 0; j < maxSc - 1; j++)            // loop that writes out all Sc parameter values
        fantosteerout << Scm[i][j] << "\t";
      fantosteerout << Scm[i][maxSc - 1] << endl;
      fantosteerout << flvcomment[i][3] << endl; // header for metamorph parameters
      for (int j = 0; j < iPts[i]; j++)          // loop that writes out all metamorph parameter values
      {
        // lk22 added routine to print out all control point PDF values.
        fantosteerout << Xs0[i][j];

        if (strScm[i][j+maxSc] == calcflag)
          fantosteerout << "\t" << (MetaRoster[iflavor[i]]->Modulator(Xs0[i][j])) - 1;
        else if (strScm[i][j+maxSc] == fixflag)
          fantosteerout << "\t" << fixflag;
        else
          fantosteerout << "\t" << Scm[i][j+maxSc];

        if ((fp[i][j] == INFINITY && fm[i][j] == -1) || (fp[i][j] == 0 && fm[i][j] == 0))
          fantosteerout << endl;
        else if (fp[i][j] != INFINITY)
          fantosteerout << "\t" << fp[i][j] << "\t" << fm[i][j] << endl;
           
      } // for (int j = 0; j < iPts[i]; j++)
    } // for (int i = 0; i < NMeta; i++)
       
  } // if (fantosteerout is_open())

  if (!fantosteerout.is_open())
  {
    cout << "Unable to open Fantomas output file" << endl;
    exit(3); // exit program if fantomas steering output card cannot be created
  }

  fantosteerout.close();

} // Writecard() ----------------------------------------------

int metamorphCollection::GetMetamorphCount()
{
  return NMeta;
} //GetMetamorphCount->

void metamorphCollection::UpdateMetamorphs()
//Update modulators for all metamorph members
{
  for (int ifl =0; ifl < NMeta; ifl++)
    MetaVector[ifl].UpdateModulator();
} //UpdateMetamorphs ----------------------------------------------

void metamorphCollection::UpdateParameters(const int ifl, double *deltas)
{
  map<int, int>::iterator it;
  it = PositionRoster.find(ifl);  // iterator used to locate if there is an entry with iflavor = ifl

  if (it == PositionRoster.end()) // returns error if entry for ifl is not found
  {
    cout << "Metamorph not found for iflavor = " << ifl << "\n"
              << "Enter initial parameters in Fantomas steering card before continuing" << endl;
    exit(3);
  }

  paraiMet = PositionRoster[ifl];

  // Update normalization Sc[0] if the delta switch DSwitch = 1
  Scm[paraiMet][0] = DSwitch[paraiMet][0]*deltas[0];

  //Update all other parameters of the carrier and modulator
  for (int i = 1; i < maxSc+Nm[paraiMet] + 1; i++)
  {
    Scm[paraiMet][i] = Scm0[paraiMet][i] + DSwitch[paraiMet][i]*deltas[i];
    if (VerbosityLevel > 0)
      cout << fixed << "updated Scm[" << paraiMet << "][" << i << "] = " << Scm[paraiMet][i]
	   << ", deltas["<<i<<"] = "<< deltas[i] << endl; 
  }
      
  //Scm[paraiMet][0] = -1+exp(Scm0[paraiMet][0] + deltas[0]);
  //Recompute the metamorph using the updated parameters
  if (Nm[paraiMet] != 0)
    MetaVector[paraiMet].UpdateModulator();

} // Metamorph::UpdateParameters ----------------------------------------------

double metamorphCollection::f(const int ifl, const double x)
{
  double ftmp;

  map<int, metamorph*>::iterator it;
  it = MetaRoster.find(ifl);  // iterator used to locate if there is an entry with iflavor = ifl

  if (it == MetaRoster.end()) // returns error if entry for ifl is not found
  {
    cout << "Metamorph not found for iflavor = " << ifl << "\n"
              << "Enter initial parameters in Fantomas steering card before continuing" << endl;
    exit(3);
  }

  ftmp = MetaRoster[ifl]->f(x); // metamorph::f(x) acting on specified metamorph object
  return ftmp;
} // metamorphCollection::f -----------------------------------------------------------

double metamorphCollection::MellinMoment(int ifl, double MellinPower, int npts)
{
  double momenttmp = MetaRoster[ifl]->GetMellinMoment(MellinPower, npts);
  return momenttmp;
} // metamorphCollection::MellinMoment

double metamorphCollection::GetConditionNumber(int ifl)
{
  double condnumtmp = MetaRoster[ifl]->GetConditionNumber();
  return condnumtmp;
}

metamorphCollection::~metamorphCollection()
{
  PositionRoster.clear();
  MetaRoster.clear();
  MetaVector.clear();
}//metamorphCollection::~metamorphCollection ----------------------------------------------
