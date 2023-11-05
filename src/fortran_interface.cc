#include"xfitter_steer.h"
#include"BaseEvolution.h"
#include<cmath>

#include "zlib.h"
#include "stdio.h"
#include <string>

//Functions to access various things from fortran
//PDFs and alpha_s are taken from default evolution
using namespace std;
//PDFs
extern "C" void hf_get_pdfsq_(double const&x,double const&Q,double*pdfs){
  xfitter::defaultEvolutionInstance()->xfxQarray(x,Q,pdfs);//returns through *pdfs
}
extern "C" void hf_get_pdfs_(double const&x,double const&Q2,double*pdfs){
  xfitter::defaultEvolutionInstance()->xfxQarray(x,sqrt(Q2),pdfs);//returns through *pdfs
}
//alpha_s
extern "C" double hf_get_alphasq_(double const&Q){
  return xfitter::defaultEvolutionInstance()->getAlphaS(Q);
}
extern "C" double hf_get_alphas_(double const&Q2){
  return xfitter::defaultEvolutionInstance()->getAlphaS(sqrt(Q2));
}

extern "C" uLong checksum_(Bytef const *buff, size_t len) {
  uLong adler = adler32(0L, Z_NULL, 0);
  adler = adler32(adler, buff, len);
  return adler;
}

extern "C" uLong filechecksum_(const char text[], size_t size) {
  FILE *fileptr;
  Bytef *buffer;
  size_t filelen;

  string bla = string(text,size);
  bla.erase(bla.find_last_not_of(" ")+1, string::npos);
  
  fileptr = fopen(bla.c_str(), "rb");  // Open the file in binary mode

  fseek(fileptr, 0, SEEK_END);          // Jump to the end of the file
  filelen = ftell(fileptr);             // Get the current byte offset in the file
  rewind(fileptr);                      // Jump back to the beginning of the file
  
  buffer = new Bytef[filelen]; // Enough memory for the file
  int i = fread(buffer, filelen, 1, fileptr); // Read in the entire file
  fclose(fileptr); // Close the file

  uLong cs = checksum_(buffer, filelen);
  
  delete[] buffer;
  return cs;
}
 
