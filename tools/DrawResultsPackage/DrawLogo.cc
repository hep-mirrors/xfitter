#include "DrawLogo.h"

#include "logo.xpm"

#include <TImage.h>
#include <TEnv.h>

#include <string>
#include <iostream>

using namespace std;
TPad * DrawLogo(string pos)
{
  string ver = VERSION;

  TImage *logo = TImage::Create();
  logo->SetImageBuffer(logo_xpm, TImage::kXpm);

  float x, y;
  x = 0.74;
  y = 0.75;
  if (pos == "dc")
    {
      x = 0.65;
      y = 0.12;
    }

  TPad * logopad = new TPad("logopad", "", x, y, x + 0.15, y + 0.14);
  if (!logo) 
    cout << "Error, Could not find logo" << endl;
  else
    {
      TString fp = gEnv->GetValue("Root.TTFontPath", "");
      TString bc = fp + "/BlackChancery.ttf";
      TString ar = fp + "/arial.ttf";

      logo->SetConstRatio(1);
      logo->DrawText(500, 600, ver.c_str(), 200, 0, 
		     bc, TImage::kShadeBelow);
      logo->SetImageQuality(TAttImage::kImgBest);
      logopad->Draw();
      logopad->cd();
      logo->Draw("");
     }
  return logopad;
}
