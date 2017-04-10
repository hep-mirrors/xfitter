#include "xfitter_cpp_base.h"
#include <string>
#include <map>


// Maybe some other place is better for this function:
int OrderMap(std::string ord) {
  std::map <std::string,int> mOrderMap  = { {"NNNLO",4}, {"NNLO",3}, {"NLO",2}, {"LO",1} };
  if ( mOrderMap.find(ord) != mOrderMap.end() ) {
    return mOrderMap[ord];
  }
  else {
    std::string text = "F: Can not convert "+ord+" to computation order";
    hf_errlog_(17032803,text.c_str(),text.size());
    return -1;
  }
}
