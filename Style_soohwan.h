#ifndef STYLE_SOOHWAN_H
#define STYLE_SOOHWAN_H

#include <TROOT.h>
#include <string>
#include <iostream>
#include <TLatex.h>
#include <TString.h>
#include <TH1D.h>

namespace custom{

void SetgStyle(){
  gStyle->SetOptStat(0);
  gStyle->SetHistTopMargin(0.03);
  gStyle->SetOptTitle(0);
};

//class LUM
//{
//  private:
//  	TString _use, _intL, _coll;
//	double _E = 5.02;
//  public:
//  	LUM(TString use, TString intL, TString coll, double E);
//	void DrawLumi();
//};
//LUM::LUM(){
//}
//void LUM::DrawLumi(){
//  TLatex* thetext = new TLatex();
//  thetext->SetNDC();
//  thetext->SetTextSize(0.04);
//  thetext->DrawLatex(0.1, 0.93, TString::Format("CMS #it{%s}", use.Data()));
//  thetext->DrawLatex(0.60, 0.93, TString::Format("%s #sqrt{s_{NN}} = %.2f %s ", coll.Data(),E, intL.Data()));
//}
//
//
}

#endif
