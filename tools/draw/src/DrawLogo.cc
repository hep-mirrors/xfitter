#include "DrawLogo.h"

#include "logo.xpm"

#include <TImage.h>
#include <TLatex.h>
#include <TEnv.h>
#include <TCanvas.h>

#include <string>
#include <iostream>
#include <sys/stat.h>

#include <CommandParser.h>

//label position
float labx;
float laby;
float xratio;

using namespace std;

// Function to find standard font to be used
string FontToUse() {
  FILE *process;
  char buff[1024];

  process = popen("fc-match --format=%{file} LiberationSans-Regular.ttf", "r");
  if (process != NULL) {
    if (fgets(buff, sizeof(buff), process)) {
    }
  }
  pclose(process);
  // Reprot if found other font
  static bool reported = false;
  if ( ! reported ) {
    if ( string(buff).find("LiberationSans-Regular.ttf") == string::npos ) {
      cout << "Could not find LiberationSans-Regular.ttf; using font " << string(buff) << endl;
      reported = true;
    }
  }
  return string(buff);
}


TPad * DrawLogo(string pos)
{
  string ver = VERSION;

  TImage *logo = TImage::Create();
  logo->SetImageBuffer(const_cast<char**>(logo_xpm), TImage::kXpm);
  if (!logo) 
    {
      cout << "Error, Could not find logo" << endl;
      return 0;
    }

  logo->SetConstRatio(1);
  logo->SetImageQuality(TAttImage::kImgBest);

  //Draw xFitter version on logo
  if (opts.version)
    {

#ifdef BUGFIXEDROOT
      TString font = FontToUse();
      struct stat st;
      //      cout << font << endl;

      if (stat(font,&st) != 0)
	{
	  cout << "Warning, cannot find font: " << font << "; xFitter version cannot be drawn on logo "<<endl;
	  opts.version = false;
	}
      else
	{
	  // first see if there is an escale character [ to indicate color. If yes, that is a sign of the release name. Draw it separately
	  if ( ver.find("[") !=std::string::npos ) {
	    string release = ver.substr(ver.find("m ")+2,ver.find("0m")-ver.find("m ")-4);
	  ver = ver.substr(0,ver.find("[")-1);
	  // std::cout << release << " " << ver << std::endl;
	  logo->DrawText(310, 80, release.c_str(), 80, 0, 
			 font, TImage::kEmbossed);	  
	  };
	}
      logo->DrawText(170, 510, ver.c_str(), 80, 0, 
		     font, TImage::kShadeBelow);
#else
      auto text = new TText(0,0,ver.c_str());
      if (opts.resolution>1000) {
	text->SetTextSize(0.04);
      }
      else {
	text->SetTextSize(0.12);
      }
      text->SetTextFont(42);
      logo->DrawText( text, 170,590);

#endif
    }

  float dx = 0.0768 * 1.5;
  float dy = 0.0597 * 1.5;

  float x, y;
  x = 1-rmarg-0.01;
  y = 1-tmarg-0.01;
  if (pos == "bc")
    {
      x = 0.74;
      y = 0.12 + dy;
    }

  TPad * logopad = new TPad("logopad", "", x-dx, y-dy, x, y);
  logopad->SetBorderSize(0);
  logopad->Draw();
  logopad->cd();
  logo->Draw("");

  return logopad;
}

void CMS()
{
  TLatex l; //l.SetTextAlign(12);
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(xratio*0.04);
  l.DrawLatex(xratio*(lmarg+0.05), 1-tmarg+0.01, "CMS");
}

void CMSpreliminary()
{
  TLatex l; //l.SetTextAlign(12);
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(xratio*0.04);
  l.DrawLatex(xratio*(lmarg+0.05), 1-tmarg+0.01, "CMS Preliminary");
}

void ATLAS()
{
  TLatex l; //l.SetTextAlign(12);
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextSize(xratio*0.04);
  l.DrawLatex(xratio*(labx+0.05), laby-0.05, "ATLAS");
}

void ATLASpreliminary()
{
  TLatex p; 
  p.SetNDC();
  p.SetTextFont(62);
  p.SetTextSize(xratio*0.04);
  p.DrawLatex(xratio*(labx+0.19), laby-0.05, "Preliminary");
}

void ATLASinternal()
{
  TLatex p; 
  p.SetNDC();
  p.SetTextFont(42);
  p.SetTextSize(xratio*0.04);
  p.DrawLatex(xratio*(labx+0.19), laby-0.05, "Internal");
}

//CDF Run II Preliminary
TPaveText * CDFIIpreliminary()
{
  TPaveText *cdfii = new TPaveText(xratio*0.45, 1-tmarg+0.02, xratio*0.95, 0.96, "NDC");
  cdfii->AddText("CDF Run II Preliminary");
  cdfii->SetTextAlign(12);
  cdfii->SetTextFont(62);
  cdfii->SetTextSize(xratio*0.04);
  cdfii->SetBorderSize(0);
  cdfii->SetFillColor(0);
  return cdfii;
}

void DrawLabels(string pos)
{
  labx = lmarg;
  laby = 1-tmarg;
  if (pos.find("ur") != string::npos)
    if (opts.atlas)
      {
      labx = 1. - rmarg - 0.19 - 0.25;
      }
    else
      {
	labx = 1. - rmarg - 0.19 - 0.20 - 0.25;
      }
  if (pos.find("bc") != string::npos)
    {
      labx = lmarg + 0.3;
      laby = bmarg + 0.05 + 0.03;
    }

  xratio = 1.;
  if (pos.find("half") != string::npos)
    xratio = 0.5;
    
  if (opts.cms)
    CMS();
  if (opts.cmspreliminary)
    CMSpreliminary();
  if (opts.atlasinternal)
    {
      ATLAS();
      ATLASinternal();
    }
  else if (opts.atlaspreliminary)
    {
      ATLAS();
      ATLASpreliminary();
    }
  else if (opts.atlas)
    ATLAS();
  if (opts.cdfiipreliminary)
    CDFIIpreliminary()->Draw();
}
