/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/// \cond DataTable_t
/**
  \file DataTable.h
  \mainpage
  \author W. Slominski
  \version 2.02
  \date 2011-11
  
  Tools for numerical data organized in columns.
*/
/// \endcond

#ifndef _TXTTBL_H_
#define _TXTTBL_H_

/** @addtogroup gr_EDTool
 *  @{
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
// #include <valarray>
#include <ctype.h>

using namespace std;
// #include "emsg.h"
#include "genut.h"
#include "Coca.h"

#include "dbgtools.h"

typedef vector<double> DblVector;

/**
  \brief Titled column.
  \internal
  \par Important changes wrt. previous version:
    ?
*/
  // \endinternal
// oooooooooooooooooooooooooooooooooooooooo
class Column_t : public DblVector {

public:

  Column_t(size_t n=0) : DblVector(n) {}
  
  //=========================
  void Scale(double s) {
    if(s == 1.) return;
    // for(int j = 0; j < size(); j++) at(j) *= s;
    iterator vp;
    for(vp = begin(); vp != end(); vp++) *vp *= s;
  }

  //=========================
  void Add(double s) {
    if(s == 0.) return;
    // for(int j = 0; j < size(); j++) at(j) *= s;
    iterator vp;
    for(vp = begin(); vp != end(); vp++) *vp += s;
  }

  // ==========================================
  void Add(const Column_t& B, double s=1) {
    for(int j = 0; j < size(); j++) at(j) += s*B[j];
  }
  
  //===========================================
  Column_t& operator-() {
    Scale(-1);
    return *this;
  }

  //===========================================
  Column_t& operator+() {
    return *this;
  }

  //===========================================
  Column_t& operator+=(double s) {
    Add(s);
    return *this;
  }

  //===========================================
  Column_t& operator+=(const Column_t& B) {
    int n = min(size(), B.size());
    for(int j = 0; j < n; j++) at(j) += B[j];
    return *this;
  }

  //================================================================
  const Column_t operator+(double s) const {
    Column_t C(*this);
    C.Add(s);
    return C;
  }
  
  //================================================================
  Column_t operator+(const Column_t& B) const {
    Column_t C(*this);
    return C += B;
  }
  
  //===========================================
  Column_t& operator*=(double s) {
    Scale(s);
    return *this;
  }

  //===========================================
  Column_t& operator*=(const Column_t& B) {
    int n = min(size(), B.size());
    for(int j = 0; j < n; j++) at(j) *= B[j];
    return *this;
  }

  //================================================================
  const Column_t operator*(double s) const {
    // int N = A.size();
    Column_t C(*this);
    C.Scale(s);
    return C;
  }
  
  //================================================================
  // friend Column_t operator* (double s, const Column_t& B);
  // friend Column_t operator* (const Column_t& B, double s);
  // friend Column_t operator* (const Column_t& A, const Column_t& B);
  
  // //================================================================
  // Column_t& operator* (double s, const Column_t& A) {
    // return A*s;
  // }
  
  //================================================================
  Column_t operator* (const Column_t& B) const {
    // int N = A.size();
    Column_t C(*this);
    return C *= B;
  }
  
};


// //================================================================
// Column_t operator* (const Column_t& B, double s) {
  // return s*B;
// }

// //================================================================
// Column_t operator* (double s, const Column_t& B) {
  // Column_t C(B);
  // C.Scale(s);
  // return C;
// }

// //================================================================
// Column_t operator* (Column_t const& A, Column_t const& B) {
  // Column_t C(A);
  // C *= B;
  // return C;
// }



/**
  \brief Numerical table with named columns.
  \internal
  \par Important changes wrt. previous version:
    All indices are 0-based.
*/
  // \endinternal
// ooooooooooooooooooooooooooooo
class DataTable_t {
    // typedef vector<double> Column_t;
    // typedef valarray<double> Column_t;
  public:
    typedef bool (*Filter_t)(const DataTable_t& );
    typedef double (*Evaluator_t)(const DataTable_t& );

  private:
  // --- DATA
  // ____________ 
  
    char CmntChar;
    string Leader;
    TextReader_t InTx;
    vector<double> CurPoint; ///< data for the currently read row
    StringList ColNames;
    int M_precis;
    int M_iDataStart;
    bool M_isReading;
    int M_iCurRow;
    // int M_iCurCol;
    vector<int> Index;
    int M_icSplit;
    
    // bool (*InputFilter)(const DataTable_t& );
    Filter_t InputFilter;
    
  public:
    Coca_t *CC; ///< parameters read from and written to the file Header
    Coca_t *xParams; ///< extra parameters used within data sets
    vector<Column_t> Data; ///< Actual data; \c Data[icol][irow] is the number at icol,irow (0-based indices)  
    vector<bool> Accepted; ///< 'Accepted' flag

    
  // --- METHODS
  // ____________ 
  
  private:
    void ParseParams();
    
    //==============================================
    void init_0() {
      M_iDataStart = 0;
      CC = new Coca_t("Parameters",CmntChar,Leader);
      CC->Create("Columns", "").AddAlias("columns:");
      xParams = new Coca_t("Extra params",CmntChar,Leader);
      xParams->SetMode(Coca_t::Auto);
      M_precis = 6;
      InputFilter = NULL;
      M_isReading = false;
      M_icSplit = -1;
    }

    //==============================================
    void reset_index() {
      int npt = GetNrows();
      Index.resize(npt);
      for(int ir=0; ir < npt; ir++) Index[ir] = ir;
    }

  public:
    /// @name Constructors
    //@{
    //==============================================
    DataTable_t(
        /*const char *name="EDT",*/ 
        const char cc='#', ///< comment char
        const char* ldr="#$"  ///< leader
      ) 
      : M_iDataStart(0)
    {
      CmntChar = cc;
      Leader = ldr;
      init_0();
    }

    /// A copy of DataTable, optionally with only a subset of source columns.
    //==============================================
    DataTable_t(const DataTable_t& A, const string& column_names = "") {
      CmntChar = A.CmntChar;
      Leader = A.Leader;
      init_0();
      *CC = *A.CC;
      *xParams = *A.xParams;
      M_icSplit = A.M_icSplit;

      InTx.M_fname = A.InTx.M_fname;
      InTx.M_text = A.InTx.M_text;
      InTx.LineNumBase = A.InTx.LineNumBase;
      InTx.M_flags = A.InTx.M_flags;

      Xstring col_names(column_names);
      if(col_names.empty()) col_names = A.ColNames.Join(",");
      int npt,ir,ic,icA;
      Create(col_names, npt = A.GetNrows());
      for(ic=0; ic < GetNcols(); ic++) {
        icA = A.FindColumn(ColNames[ic]);
        for(ir=0; ir < npt; ir++) Data[ic][ir] = A.Data[icA][ir];
      }
    }
    
    //==============================================
    DataTable_t& operator=(const DataTable_t& A) {
      CmntChar = A.CmntChar;
      Leader = A.Leader;
      init_0();
      *CC = *A.CC;
      *xParams = *A.xParams;
      M_icSplit = A.M_icSplit;

      InTx.M_fname = A.InTx.M_fname;
      InTx.M_text = A.InTx.M_text;
      InTx.LineNumBase = A.InTx.LineNumBase;
      InTx.M_flags = A.InTx.M_flags;

      int npt,ir,ic;
      Create(A.ColNames.Join(","), npt = A.GetNrows());
      Accepted = A.Accepted;
      Index = A.Index;
      for(ic=0; ic < GetNcols(); ic++) {
        for(ir=0; ir < npt; ir++) Data[ic][ir] = A.Data[ic][ir];
      }
    }
    
    //@}

    /** @name Creation, modification and filling
     Methods to:
     \li create space for values (initialized to 0),
     \li (re)name columns
     \li modify layout
     \li assign values
    */
    //@{
    //=========================================
    void Create(int nCols=0, int nRows=0) {
      if(!nCols) nCols = ColNames.size();
      // if(nCols == GetNcols() && nRows == GetNrows()) return; //--- this leaves old data and old Accepted!
      Clear();
      Data.resize(nCols);
      if(nRows) {
        for(int ic=0; ic < nCols; ic++) Data[ic].resize(nRows);
        Accepted.assign(GetNrows(), true);
        reset_index();
      }
    }

    //==================================================
    void Create(const string& col_names, int nRows=0) {
      Clear(); //--- needed for the final check in SetColNames
      SetColNames(col_names);
      Create(0, nRows);
    }

  private:
    //==================================================
    StringList ExpandColRange(const Xstring& rname) {
      // abc[j:k]xyz --> A,j,k,B
      int j,k;
      StringList cList;
      if(!rname.Contains("[")) {cList.push_back(rname); return cList;}
      if(!rname.Contains("]")) throw Fiasco("Syntax error in '%s'", rname.data());
      if(!rname.Contains(":")) throw Fiasco("Syntax error in '%s'", rname.data());
      Xstring rn(rname);
      Xstring nic(" ");
      rn = nic + rn;
      string A(rn.Token("[").Trim());
      if(A.empty() || !isalpha(A[0])) throw Fiasco("Syntax error in '%s'", rname.data());
      
      rn.erase(0,1);
      rn = nic + rn;
      Xstring sj(rn.Token(":").Trim());
      if(sj.empty()) j = 1; else j = sj.GetInt();
      
      rn.erase(0,1);
      rn = nic + rn;
      Xstring sk(rn.Token("]").Trim());
      if(sk.empty()) k = 0; else k = sk.GetInt();
      
      rn.erase(0,1); //--- B
      // string B(rn.Remainder().erase(0,1));

      // #define SHO(a) cout << #a" = |" << (a) <<"|"<< endl;
      // SHO(A) 
      // SHO(sj) 
      // SHO(sk) 
      // SHO(rn) 
      for(; j <= k; j++) {
        char ts[8];
        sprintf(ts, "%d",j);
        cList.push_back(A+ts+rn);
      }
      return cList;
    }
    
  public:
    /**
      \brief Set column names
      
      Column names cannot contain whitespace, comma, semicolon.\n
      Special are: 
      \arg # as the last char; illegal when \c AsymmErrors not defined. See ExpData_t.
      \arg [\em start:end] denoting a sequence, e.g.\n
      <tt> a[1:23]</tt> translates to \c a1,a2,...,a23
    */
    //=========================
    void SetColNames(const string& ColNamesList = "" ///< Column names separated by: whitespace, comma, semicolon
                    ) {
      // ColNames.clear();
      // CC->AddParam("Columns", xstr.Join(ColNames,",").c_str());
      // TParam t;
      Xstring xs(ColNamesList);
      if(xs.empty()) xs = CC->GetString("Columns");
      else CC->SetVal("Columns", xs);
      
      ColNames = xs.Split(",; \t");
      
      //--- expand ranges "[j:k]"
      
      StringList::iterator cnp;
      if(xs.Contains("[")) {
        StringList newCN;
        for(cnp=ColNames.begin(); cnp != ColNames.end(); cnp++) newCN.Append(ExpandColRange(*cnp));
        ColNames = newCN;
      }
      
      //--- resolve "name#"
      
      if(xs.Contains("#")) {
        int asy = 0;
        if(CC->Exists("AsymmErrors")) asy = CC->GetInt("AsymmErrors");
        else throw "'#' found in column names but 'AsymmErrors' not defined.";
        // asy = CC->GetInt("AsymmErrors");
        if(asy) {
          asy = (asy > 0) ? 1 : 0;
          char pm[] = "-+";
          StringList newCN;
          for(cnp=ColNames.begin(); cnp != ColNames.end(); cnp++) {
            int last = cnp->size()-1;
            if((*cnp)[last] == '#') {
              string s(cnp->substr(0,last));
              newCN.push_back(s+pm[asy]);
              newCN.push_back(s+pm[1-asy]);
            }
            else newCN.push_back(*cnp);
          }
          ColNames = newCN;
        }
      }
      
      // --- check and auto-name
      
      if(GetNcols() > ColNames.size()) {
        int nstart = ColNames.size() +1;
        int nend = GetNcols();
        char ts[8];
        string autostr("@");
        for(int n = nstart; n <= nend; n++) {
          sprintf(ts, "%d", n);
          ColNames.push_back(autostr + ts);
        }
      }
    }

    //=========================
    StringList GetColNames(size_t start=0, int count=-1) const {
      // vector<Xstring> A(ColNames.begin()+start, ColNames.begin()+start+count);
      // if(count >= 0) return StringList(vector<Xstring>(ColNames.begin()+start, ColNames.begin()+start+count));
      // return StringList(ColNames.begin()+start, ColNames.end());
      return ColNames.Range(start, count);
    }
    
    //=========================
    void Rename(const Xstring& cold, const Xstring& cnew) { 
      if(cnew.IsBlank()) throw Fiasco("Rename: column name cannot be blank.");
      ColNames[FindColumn(cold)] = cnew;
      CC->SetVal("Columns", ColNames.Join(","));
    }

    // /// Scale column \c iCol by \c s
    // //=========================
    // void Scale(int iCol, double s) {
      // if(s == 1.) return;
      // // int npt = GetNrows();
      // for(M_iCurRow = GetNrows()-1; M_iCurRow >= 0; M_iCurRow--) Data[iCol][M_iCurRow] *= s;
    // }
    
    // /// Scale column \c name by \c s
    // //=========================
    // void Scale(const string& name, double s) {Scale(FindColumn(name), s);}
    
    //=========================
    void Fill(int iCol, Evaluator_t fcn) {
      // iCol--;
      int npt = GetNrows();
      for(M_iCurRow=0; M_iCurRow < npt; M_iCurRow++) Data[iCol][M_iCurRow] = fcn(*this);
    }
    
    //=========================
    void Fill(const string& name, Evaluator_t fcn) { Fill(FindColumn(name), fcn); }
    
    //=========================
    DataTable_t& AddColumn(const Xstring& cnam, const string& before = Xstring::null_string) {
      if(before.empty()) {
        Data.push_back(Column_t(GetNrows()));
        ColNames.push_back(cnam);
        return *this;
      } 
      else return AddColumn(cnam, FindColumn(before));
    }

    //=========================
    DataTable_t& AddColumn(const Xstring& cnam, int before) {
      //--- DO NOT return column index --> it can change after adding another column
      // if(cnam.IsBlank()) throw Fiasco("AddColumn: column name cannot be blank.");
      if(FindColumn(cnam, 0) >= 0) throw Fiasco("AddColumn: column '%s' already exists.", cnam.c_str());
      Column_t v(GetNrows());
      Data.insert(Data.begin()+before, v);
      ColNames.insert(ColNames.begin()+before, cnam);
      CC->SetVal("Columns", ColNames.Join(","));
      return *this;
    }

    /// Add row filled with zeroes
    //=========================
    size_t AddRow() {
      int Ncols = GetNcols();
      if(!Ncols) throw Fiasco("Cannot add rows to non-existent Data.");
      for(int col=0; col < Ncols; col++) Data[col].push_back(0);
      Accepted.push_back(true);
      return GetNrows()-1;
    }

    //=========================
    size_t AddRow(const double* vals) {
      int Ncols = GetNcols();
      for(int col=0; col < Ncols; col++) Data[col].push_back(vals[col]);
      return GetNrows()-1;
    }

    //=========================
    size_t AddRow(const vector<double>& vals) {
      int Ncols = GetNcols();
      for(int col=0; col < Ncols; col++) Data[col].push_back(vals[col]);
      return GetNrows()-1;
    }

    /**
      \brief Swap two columns (contents and names).\n
      Swaps columns' contents as well as their names.
      Note that the standard \c swap(obj[iCol1],obj[iCol2]) swaps only contents w/o the names.
    */
    //==============================================
    DataTable_t& SwapColumns(int ia, int ib) {
      // swap(df["systot+"], df["systot-"]);
      swap(Data[ia], Data[ib]);
      swap(ColNames[ia], ColNames[ib]);
      return *this;
    }

    //==============================================
    DataTable_t& SwapColumns(const string& cnamA, const string& cnamB) {
      int icA = FindColumn(cnamA);
      int icB = FindColumn(cnamB);
      return SwapColumns(icA, icB);
    }

    // //==============================================
    // DataTable_t& MoveColumn(int icol, int before) {
      // Data.erase(Data.begin()+icol);
      // ColNames.erase(ColNames.begin()+icol);
      // CC->SetVal("Columns", ColNames.Join(","));
      // return *this;
    // }

    //==============================================
    DataTable_t& RemoveColumn(int icol) {
      // Data.at(icol); //--- generate potential out_of_range exception
      Data.erase(Data.begin()+icol);
      ColNames.erase(ColNames.begin()+icol);
      CC->SetVal("Columns", ColNames.Join(","));
      return *this;
    }

    //==============================================
    DataTable_t& RemoveColumn(const string& cnam) {
      return RemoveColumn(FindColumn(cnam));
    }

    // /// Remove \c N columns starting from \c icol
    // //==============================================
    // DataTable_t& RemoveColumn(int icol, int n=1) {
      // if(!n) return *this;
      // // if(n < 0) return;
      // int last = icol+n;
      // // Data.at(icol); //--- generate potential out_of_range exception
      // Data.erase(Data.begin()+icol, n > 0 ? Data.begin()+icol+n);
      // ColNames.erase(ColNames.begin()+icol, ColNames.begin()+last);
      // CC->SetVal("Columns", ColNames.Join(","));
      // return *this;
    // }

    /// Discard all Data, but not parameters.
    //=========================
    void Clear() {
      for(int col=0; col < GetNcols(); col++) Data[col].clear();
      Data.clear();
    }

    /// Fill \c Accepted according to \c filt
    //=========================
    int Select(Filter_t filt) {
      // if(filt) 
      int n=0;
      for(M_iCurRow=0; M_iCurRow < GetNrows(); M_iCurRow++) if((Accepted[M_iCurRow] = filt(*this))) n++;
      // else Accepted.assign(GetNrows(), true);
      return n;
    }
    
    /// Set \c Accepted for all data points to \c all
    //=========================
    void Select(bool all=true) {
      Accepted.assign(GetNrows(), all);
    }
    
    /// Stable sort by rows acc. to given columns. \c Sort() w/o arg. restores the original order.
    //=========================
    void Sort(const string& column_names="" ///< comma-separated column names
    );
    // {
      // M_iCurCol = FindColumn(cname);
      // gVal = &Data[M_iCurCol];
      // int npt = GetNrows(), ir;
      // Index.resize(npt);
      // for(ir=0; ir < npt; ir++) Index[ir] = ir;
      // stable_sort (Index.begin(), Index.end(), cmp_vals);
    // }
    
    //@}

    //=========================
    // bool skipblocks(int bcnt);

    //=========================
    // bool skipcomments();

private:
    //=========================
    void read_text_file(const char *fname) {
      InTx.Open(fname);
      InTx.Read();
      InTx.Close();
    }

    //======================================================
    void ParseHeader(const char *SOHpat=NULL, const char *EOHpat="") {
      // string::npos = (unsigned)(-1)
      int istart=0;
      if(SOHpat) {
        istart = InTx.Find(SOHpat);
        if(istart < 0) throw Fiasco("Start-Of-Header pattern '%s' not found.", SOHpat);
      }
      if(EOHpat) {
        int iend = InTx.Find(EOHpat, istart);
        if(iend < 0) throw Fiasco("End-Of-Header pattern '%s' not found.", EOHpat);
        M_iDataStart = iend;
        CC->ParseList(InTx, istart, iend - istart);
      }
      else CC->ParseList(InTx, istart);
      #ifdef DBGT_ON
      CC->Show();
      #endif
      SetColNames(CC->GetString("Columns")); //--- will be called once more by ReadDataSet but needed here for ExpData_t::ReadFile
    }

    //=========================
    int findparam(const string& pnam, double pval, double acc=1e-6) {
      // string sep(WHSP+"=");
      string sep(Xstring::white_space + "=");
      int N=InTx.NLines();
      // string pat(Leader)
      for(int j=M_iDataStart; j < N; j++) {
        // SHOW(InTx[j])
        if(!InTx[j].BeginsWith(Leader)) continue;
        Xstring coca(InTx[j]);
        coca.erase(0, Leader.size());
        string nam = coca.Token(sep);
        if(nam != pnam) continue;
        coca.TrimRight();
        double v = coca.Token().GetReal(); 
        xParams->Define(pnam, v);
        if(pval ? (fabs(v/pval -1) < acc) : (v == 0.0)) return j;
      }
      return -1;
    }

    //=========================
    int FindDblBlank(int start) {
      // --- returns index of the first blank line
      int N=InTx.NLines()-1;
      for(int j=start; j < N; j++) {
        if(InTx[j].empty() && InTx[j+1].empty()) return j;
      }
      return -1;
    }

    
public:
    /** @name Reading data from files
    */
    //@{
    //=========================
    void SetInputFilter(Filter_t ifilter) {
      InputFilter = ifilter;
    }

    //======================================================
    int ReadFile(const string& fname, const char *SOHpat=NULL, const char *EOHpat="") {
      SetInputFile(fname, SOHpat, EOHpat);
      return ReadDataSet();
    }
    
    //======================================================
    void SetInputFile(const string& fname, const char *SOHpat=NULL, const char *EOHpat="") {
      cout << "Opening file '" << fname <<"'"<< endl;
      read_text_file(fname.c_str());
      ParseHeader(SOHpat, EOHpat);
    }
    
    /**
      \internal
      \c Data array is resized according to the data read.
      \c
    */
    //=====================================
    int ReadDataSet(const string& key="", double val=0, double  acc=1e-6) {
      int Ncols=0;
      M_isReading = true;
      xParams->clear();
      bool isKey = !key.empty();
      int j, j0 = M_iDataStart, j1=InTx.NLines();
      DBG_SHOW(j1)
      if(isKey) {
        j0 = findparam(key, val, acc);
        if(j0 < 0) return 0;
        j0++;
        j1 = FindDblBlank(j0);
        if(j1 < 0) j1=InTx.NLines();
      }
      // markpos();
      // SHOW(pos)  SHOW(lineno)
      Clear();
      for(j=j0; j < j1; j++) {
        Xstring dstr(InTx[j]);
        // cout << j<<":  "<< dstr << endl;
        if(dstr.IsBlank()) continue;
        if(isKey && dstr.BeginsWith(Leader)) xParams->ParseList(InTx, j, 1);
        else if(dstr[0] == CmntChar) continue;
        else {
          StringList vals = dstr.Split();
          if(!Ncols) {
            Ncols = vals.size();
            Data.resize(Ncols);
            CurPoint.resize(Ncols);
          }
          else if(Ncols != vals.size()) throw Fiasco("Different # columns at line %i:\n'%s'", j+1, dstr.data());
          int col;
          // CurPoint.clear();
          // for(col=0; col < Ncols; col++) CurPoint.push_back(vals[col].GetReal());
          for(col=0; col < Ncols; col++) CurPoint[col] = vals[col].GetReal();
          if(!InputFilter || InputFilter(*this))
            // for(col=0; col < Ncols; col++) Data[col].push_back(vals[col].GetReal());
            for(col=0; col < Ncols; col++) Data[col].push_back(CurPoint[col]);
        }
      }
      if(ColNames.size() && ColNames.size() != Ncols) throw Fiasco("Different # columns and Column names: %d != %d", Ncols, ColNames.size());
      Accepted.assign(GetNrows(), true);
      reset_index();
      M_isReading = false;
      SetColNames();
      return GetNrows();
    }
    
    //@}

    /// @name Accessing the data
    //@{
    //================================================================
    Column_t& operator[] (const string& name) {return Data.at(FindColumn(name));}
    Column_t& operator[] (int ic) {return Data.at(ic);}

    /** \brief CurrentPoint value by name.
    
      Provides access to current values in \c Filter_t and \c Evaluator_t functions.
    */
    //=========================
    double Value(const string& name) const {return Value(FindColumn(name));}

    /// CurrentPoint value by index (0-based).
    //=========================
    double Value(int ic ///< 0-based column index
    ) const {
      return M_isReading ? CurPoint.at(ic) : Data.at(ic).at(M_iCurRow);
    }
    
    //=========================
    size_t GetNcols() const {
      return Data.size();
    }

    //=========================
    size_t GetNrows() const {
      return Data[0].size();
    }

    //=========================
    /// returns 0-based index
    int FindColumn(const Xstring& cnam, bool strict=true) const { 
      if(cnam.IsBlank()) throw Fiasco("FindColumn: column name cannot be blank.");
      for(int j=0; j < ColNames.size(); j++) {
        if(ColNames[j] == cnam) return j;
      }
      if(strict) throw Fiasco("Column '%s' not found.", cnam.c_str());
      return -1;
    }

    //@}
    
    
    // --- OUTPUT ---
    // _________________
    
    /** @name Saving (printing)
    */
    //@{
    //================================================================
    // void out_row(ostream &ostr=cout) {
      // for(ic=0; ic < GetNcols(); ic++) ostr << setw(cw) << left << Data[ic][ir];
    // }
    
    //================================================================
    void Show(const string& fname, const string& Out_columns = "") {
      ofstream vf(fname.data());
      Show(vf, Out_columns);
      vf.close();
    }
    
    /**
      \internal
      \todo fix the column width basing on the formatted data and title
    */
    //================================================================
    void Show(ostream &ostr=cout, const string& Out_columns = "") {
    
      if(!GetNcols()) return;
      
      int ir,ic;
      int npt = GetNrows();
      
      // --- Output columns selection
      StringList OutColNames;
      Xstring out_columns(Out_columns);
      if(out_columns.empty()) OutColNames = ColNames;
      else OutColNames = out_columns.Split(",");
      
      // --- split into sets
      bool split = M_icSplit >= 0;
      double lastval;
      ostringstream Svals;
      // int nout=0;
      int nCC = CC->size();
      if(split) {
        Sort(ColNames[M_icSplit]);
        int cnt=0;
        lastval = Data[M_icSplit][Index[0]] -1e99;
        // Svals << (lastval = Data[M_icSplit][Index[0]]);
        for(int irow=0; irow < npt; irow++) {
          ir = Index[irow];
          if(!Accepted[ir]) continue;
          // DBG_SHOW(ir)
          if(lastval != Data[M_icSplit][ir]) {
            if(cnt) Svals << ",";
            Svals << (lastval = Data[M_icSplit][ir]);
            cnt++;
          }
        }
        string vs("Values(");
        CC->Define(vs+ColNames[M_icSplit]+")", Svals.str());
        vs = "Count(";
        CC->Define(vs+ColNames[M_icSplit]+")", cnt);
      }
      
      
      // --- Header
      CC->SetTail("--- end-of-header");
      // CC->SetTail("");
      CC->SetVal("Columns", OutColNames.Join(","));
      CC->Show(ostr);
      // ostr << CmntChar << endl;
      ostr << endl;

      // --- Data header
      int cw = M_precis+7;
      StringList::iterator iter;
      ostringstream tits,nums;
      tits.setf(ios::left);
      nums.setf(ios::left);
      char tmps[32];
      // bool nxt=false;
      // ostr << CmntChar << " ";
      for(ic=1,iter = OutColNames.begin(); iter != OutColNames.end(); iter++,ic++ ) {
        tits << "  " << setw(cw) << *iter;
        nums << "  :" << setw(cw-1) << ic;
        // tits << "  " << setw(cw) << iter->Justify(0,cw);
        // sprintf(tmps,"%d", ic);  nums << "  " << setw(cw) << Xstring(tmps).Justify(0,cw);
      }
      // ostr << endl;
      tits.seekp(0).put(CmntChar);
      nums.seekp(0).put(CmntChar);
      
      CC->SetVal("Columns", ColNames.Join(","));
      
      // xParams->ShowPlain(ostr);
      #define OUTSHD \
      ostr << *xParams;\
      ostr << tits.str() << endl;\
      ostr << nums.str() << endl;
      
      // --- Data
      cw = M_precis+9;
      // ostr.setf(ios::scientific);// | ios::internal);
      ostr.setf(ios::right);
      ostr.precision(M_precis);
      DBG_SHOW(npt)
      
      StringList::iterator cnp;
      // if(split) xParams->Define(ColNames[M_icSplit],lastval = Data[M_icSplit][Index[0]]);
      // OUTSHD
      if(split) lastval = Data[M_icSplit][Index[0]] -1e99;
      else {OUTSHD}
      DBG_HERE("Index?")
      int nout=0;
      for(int irow=0; irow < npt; irow++) {
        ir = Index[irow];
        if(!Accepted[ir]) continue;
        // DBG_SHOW(ir)
        if(split && lastval != Data[M_icSplit][ir]) {
          if(nout) ostr << "\n\n";
          // lastval = Data[M_icSplit][ir];
          xParams->Define(ColNames[M_icSplit],lastval = Data[M_icSplit][ir]);
          OUTSHD
        }
        nout++;
        ostr << "  ";
        for(cnp=OutColNames.begin(); cnp != OutColNames.end(); cnp++)  ostr << setw(cw) << left << Data.at(FindColumn(*cnp))[ir];
        // for(ic=0; ic < GetNcols(); ic++) ostr << setw(cw) << left << Data[ic][ir];
        ostr << endl;
      }
      #undef OUTSHD
      xParams->Delete(ColNames[M_icSplit]);
      if(split) CC->resize(nCC);
    }
    
    /// Sets the precision of numbers printed by \c Show.
    //=========================
    void SetPrecision(int p) { M_precis = p; }
    
    /// Output of \c Show will be splitted into sets of constant values of "name".
    //=========================
    void SplitOutput(const string& name) { M_icSplit = FindColumn(name); }
    
    //@}
    
    friend class ExpData_t;

  
};

// int ReadDataTbl(const char* dname, const char* key, double val,
              // double* &xin, double** fiin, int ncols);

// int ReadDataTbl(char* dname, char* key, double val,
              // double** fiin, vector<int>& cols);

///@}

#endif
