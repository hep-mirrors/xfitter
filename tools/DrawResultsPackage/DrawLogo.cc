#include "DrawLogo.h"

#include "logo.xpm"

#include <TImage.h>
#include <TASImage.h>
#include "TEnv.h"
//#include "TString.h"


#include <string>
#include <iostream>

using namespace std;
TPad * DrawLogo(void)
{
  string ver = VERSION;

  TImage *logo = TImage::Create();
  logo->SetImageBuffer(logo_xpm, TImage::kXpm);


  TPad * logopad = new TPad("logopad", "", 0.74, 0.75, 0.89, 0.89);
  //TPad * logopad = new TPad("logopad", "", 0.64, 0.75, 0.79, 0.89);
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
