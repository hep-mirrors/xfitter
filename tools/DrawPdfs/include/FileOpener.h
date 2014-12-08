/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2014
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 1.0.0
  
_____________________________________________________________*/

#ifndef FILE_OPENER_H_
#define FILE_OPENER_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

// -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
class InFileOpener_t {
  int Verbose;
  ifstream f;
  int ind;
    
  public:

  typedef vector<string> flist_t;
  flist_t Flist;

  // ==================================
  InFileOpener_t() {
    Verbose = 1;
    // mode = opmode;
  }
  
  // // ==================================
  // FileOpener_t(const flist_t& fl) : Flist(fl) {}
  
  // ==================================
  ~InFileOpener_t() {
    Close();
    // cout << "Bye! " << f.is_open() << endl;
  }
  
  // ==================================
  void Close() {
    // if(f.is_open()) // --- not needed
    f.close();
  }
  
  // ==================================
  void SetVerbose(int verb) {Verbose = verb;}
  
  // ==================================
  ifstream& GetStream() {return f;}
  
  // ==================================
  int GetIndex() const {return ind;}
  
  // ==================================
  string GetPath() const {return ind < 0 ? "" : Flist[ind];}
  
  // ==================================
  int Add(const string& fname) {
    Flist.push_back(fname);
  }

  // ==================================
  void Clear() {
    Flist.clear();
  }

  // ==================================
  /**
    \brief Open first available file from \c Flist.
    
    Return  0 on success
  */
  int Open() {
    ind = -1;
    if(Flist.empty()) {
      if(Verbose) cout <<"InFileOpener_t: No files to open." << endl;
      return 2;
    }
    Close();
    vector<string>::iterator it;
    for(it = Flist.begin(); it != Flist.end(); it++) {
      if(Verbose > 1) cout <<"Trying to open '" << *it <<"' ... " << flush;
      // f.open(it->c_str(), mode);
      f.open(it->c_str());
      if(!f) {
        if(Verbose > 1) cout <<" Failed." << endl;
      } else {
        if(Verbose > 1) cout <<" OK." << endl;
        ind = it - Flist.begin();
        return 0;
      }
    }
    
    if(Verbose) {
      cout << "InFileOpener_t: Could not open any of\n";
      Show();
    }
    return 1;
  }

  // ==================================
  void Show() {
    vector<string>::iterator it;
    for(it = Flist.begin(); it != Flist.end(); it++) {
      cout <<"  "<< *it << '\n';
    }
  }

};

#endif
