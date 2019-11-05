#include <iostream>
#include <fstream>
#include <cstring>
#include "xfitter_cpp_base.h"
#include <string>
#include <map>
//Some generally useful utility functions
//Maybe we should move this functions to some other places that makes more sense?
using std::string;
using std::map;
using std::cerr;
using std::endl;

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
void stringToFortran(char* dest, size_t maxlen, const string& s){
  size_t len=s.size();
  if(len>maxlen){
    cerr<<"[ERROR] Failed to convert c++ string \""<<s<<"\" to Fortran character array: string length ("<<len<<") larger than size of destination array ("<<maxlen<<");"<<endl;
    hf_errlog(19081610, "F: C++ string too long to convert to Fortran, see stderr");
  }
  memcpy(dest, s.c_str(), len);
  char* end=dest+maxlen;
  for(char* p=dest+len;p<end;++p) *p=' ';
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

bool fileExists(const string&name){
  return std::ifstream(name).good();
}

extern "C"{
  extern struct {
    char outdirname[256]; // output dir name
  } coutdirname_;
}

string xfitter::getOutDirName(){
  return stringFromFortran(coutdirname_.outdirname,256);
}
