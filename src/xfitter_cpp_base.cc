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


// Maybe some other place is better for this function:
void hf_errlog(int id,const std::string& message) {
   hf_errlog_(id,message.c_str(),message.size());
}
std::string stringFromFortran(char*s,size_t size){
  char*p=s+size;
  while(p>s){
    if(*(--p)!=' '){
      ++p;
      break;
    }
  }
  return std::string(s,p-s);
}
