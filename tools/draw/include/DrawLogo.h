#ifndef DrawLog_h
#define DrawLog_h

#include <TPad.h>
#include <TPaveText.h>

#include <string>
using namespace std;

TPad * DrawLogo(string pos = "ul");
void DrawLabels(string pos = "ul");

#endif
