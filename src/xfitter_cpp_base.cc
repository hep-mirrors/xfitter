#include "xfitter_cpp_base.h"
#include <string>
#include <map>
//Some generally useful utility functions
//Maybe we should move this functions to some other places that makes more sense?
using std::string;
using std::map;

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
bool beginsWith(const string&s,const string&p){
  return s.compare(0,p.size(),p)==0;
}
bool beginsWith(const char*s,const char*p){
  while(*p){
    if(*p!=*s)return false;
    s++;
    p++;
  }
  return true;
}
bool beginsWith(const string&s,const char*p){
  return beginsWith(s.c_str(),p);
}
void stripString(string&s){
  const char*p1,*p2;
  p1=s.c_str();
  p2=p1+s.size()-1;
  while(*p1==' ')++p1;
  while(*p2==' ')--p2;
  if(p1==s.c_str()&&p2==p1+s.size()-1)return;
  s=s.substr(size_t(p1-s.c_str()),size_t(p2-p1+1));
}

extern "C"{
  extern struct {
    char outdirname[256]; // outout dir name
    char lhapdf6outdir[256]; // DEPRECATED
  } coutdirname_;
}

string xfitter::getOutDirName(){
  return stringFromFortran(coutdirname_.outdirname,256);
}
