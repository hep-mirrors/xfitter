/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstdarg>
//#include <io.h>
//#include <dir.h>
//#include <conio.h>
#include "emsg.h"
using namespace std;

FILE *log_file=stderr;
string log_name;

//*************************************************************
void OpenLog(char* lonam) {
  log_name = lonam;
  log_file = fopen(lonam, "w");
}

//*************************************************************
void CloseLog() {
  fclose(log_file);
  log_name.clear();
  log_file = stderr;
}

//*************************************************************
void vemsg(int flag, const char *fmt, va_list apt) {
  fprintf(log_file,"\n");
  vfprintf(log_file,fmt,apt);
  fprintf(log_file,"\n");
  fflush(log_file);
  if(flag && log_file!=stderr) {
    fprintf(stderr,"\n");
    vfprintf(stderr,fmt,apt);
    fprintf(stderr,"\n");
    fflush(stderr);
  }
  if(flag>0) {
    fprintf(stderr,"\nFATAL error - program terminated with exit code = %i\n",flag);
    if(log_file!=stderr)
      fprintf(log_file,"\nFATAL error - program terminated with exit code = %i\n",flag);
    //exit(flag);
  }
}

//*************************************************************
void emsg(int flag, const char *fmt, ...) {
  /*
    Write message to stderr and log_file.
    flag values:
    0 log_file only,
    >0 both files, FATAL error, exit(flag),
    <0 both files.
*/
  va_list apt;
  va_start(apt,fmt);
  vemsg(flag, fmt, apt);
  va_end(apt);
  if(flag > 0) exit(flag);
}

//*************************************************************
void require(bool cond, const char *fmt, ...) {
  if(cond) return;
  va_list apt;
  va_start(apt,fmt);
  vemsg(99, fmt, apt);
  va_end(apt);
  exit(99);
}
