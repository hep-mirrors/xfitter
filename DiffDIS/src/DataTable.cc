/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "DataTable.h"

static vector<double> *gVal;

//=========================================
bool cmp_vals(int i, int j) {return gVal->at(i) < gVal->at(j); }

/**
  Stable sort.\n
  It is also indirect — the sorted row indices are kept in a private table \c Index.
  Nb. <tt> A[icol][irow] </tt> refers to an unsorted row index \c irow.\n
  Call \c Sort() w/o argument to recover the original order.
*/
//=========================================
void DataTable_t::Sort(const string& column_names) {
  if(column_names.empty()) {reset_index(); return;}
  // int M_iCurCol = FindColumn(cname);
  // gVal = &Data[M_iCurCol];
  int npt = GetNrows(), ir;
  if(Index.size() != npt) reset_index();
  
  StringList::reverse_iterator cnp;
  StringList clist = Xstring(column_names).Split(",");
  for(cnp = clist.rbegin(); cnp != clist.rend(); cnp++) {
    gVal = &Data.at(FindColumn(*cnp));
    stable_sort (Index.begin(), Index.end(), cmp_vals);
  }
}
