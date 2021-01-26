#include <TFile.h>
#include "Style_soohwan.h"

void plotMBtrigEff(){
  TFile* fT13 = new TFile("MBPD_Ev6p7M_TrigEff_L2_v3.root","read");
  TFile* fT14 = new TFile("MBPD_Ev6p7M_TrigEff_L3_v3.root","read");
  
  TH1D* ht13p1 = (TH1D*) fT13->Get("hEff1");
  TH1D* ht13p2 = (TH1D*) fT13->Get("hEff2");
  TH1D* ht13p3 = (TH1D*) fT13->Get("hEff3");
  TH1D* ht13p4 = (TH1D*) fT13->Get("hEff4");
  TH1D* ht14p1 = (TH1D*) fT14->Get("hEff1");
  TH1D* ht14p2 = (TH1D*) fT14->Get("hEff2");
  TH1D* ht14p3 = (TH1D*) fT14->Get("hEff3");
  TH1D* ht14p4 = (TH1D*) fT14->Get("hEff4");
  Double_t epr[5] = {0.0, 1.2, 1.8, 2.1, 2.4};
  Double_t ptc[4] = {3.5, 2.1, 1.5, 1.5};
  ht13p1->GetYaxis()->SetRangeUser(0,1);
  ht13p2->GetYaxis()->SetRangeUser(0,1);
  ht13p3->GetYaxis()->SetRangeUser(0,1);
  ht13p4->GetYaxis()->SetRangeUser(0,1);
  ht14p1->GetYaxis()->SetRangeUser(0,1);
  ht14p2->GetYaxis()->SetRangeUser(0,1);
  ht14p3->GetYaxis()->SetRangeUser(0,1);
  ht14p4->GetYaxis()->SetRangeUser(0,1);

  ht13p1->GetXaxis()->SetRangeUser(0,30);
  ht13p2->GetXaxis()->SetRangeUser(0,30);
  ht13p3->GetXaxis()->SetRangeUser(0,30);
  ht13p4->GetXaxis()->SetRangeUser(0,30);
  ht14p1->GetXaxis()->SetRangeUser(0,30);
  ht14p2->GetXaxis()->SetRangeUser(0,30);
  ht14p3->GetXaxis()->SetRangeUser(0,30);
  ht14p4->GetXaxis()->SetRangeUser(0,30);

  ht13p1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht13p2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht13p3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht13p4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht14p1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht14p2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht14p3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ht14p4->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  ht13p1->GetYaxis()->SetTitle("Single #mu Eff.");
  ht13p2->GetYaxis()->SetTitle("Single #mu Eff.");
  ht13p3->GetYaxis()->SetTitle("Single #mu Eff.");
  ht13p4->GetYaxis()->SetTitle("Single #mu Eff.");
  ht14p1->GetYaxis()->SetTitle("Single #mu Eff.");
  ht14p2->GetYaxis()->SetTitle("Single #mu Eff.");
  ht14p3->GetYaxis()->SetTitle("Single #mu Eff.");
  ht14p4->GetYaxis()->SetTitle("Single #mu Eff.");

  custom::SetgStyle();

  TCanvas *c1, *c2;
  c1 = new TCanvas("c1","Jpsi L2 Filter Eff",1600,1200);
  c2 = new TCanvas("c2","Jpsi L3 Filter Eff",1600,1200);

  TLatex* lum = new TLatex();
  lum->SetNDC();
  TLatex* ds = new TLatex();
  ds->SetTextSize(0.03);
  ds->SetNDC();

  c1->Divide(2,2);
  c2->Divide(2,2);

  c1->cd(1);
  ht13p1->SetMarkerSize(1);
  ht13p1->SetMarkerStyle(kCircle);
  ht13p1->SetMarkerColor(kRed);
  ht13p1->Draw("pe");
  lum->DrawLatex(0.60,0.55,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.50,"Data L2 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.45,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[0],epr[0],epr[1]));
  c1->cd(2);
  ht13p2->SetMarkerSize(1);
  ht13p2->SetMarkerStyle(kCircle);
  ht13p2->SetMarkerColor(kRed);
  ht13p2->Draw("pe");
  lum->DrawLatex(0.60,0.55,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.50,"Data L2 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.45,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[1],epr[1],epr[2]));
  c1->cd(3);
  ht13p3->SetMarkerSize(1);
  ht13p3->SetMarkerStyle(kCircle);
  ht13p3->SetMarkerColor(kRed);
  ht13p3->Draw("pe");
  lum->DrawLatex(0.60,0.55,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.50,"Data L2 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.45,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[2],epr[2],epr[3]));
  c1->cd(4);
  ht13p4->SetMarkerSize(1);
  ht13p4->SetMarkerStyle(kCircle);
  ht13p4->SetMarkerColor(kRed);
  ht13p4->Draw("pe");
  lum->DrawLatex(0.60,0.55,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.50,"Data L2 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.45,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[3],epr[3],epr[4]));

  c2->cd(1);
  ht14p1->SetMarkerSize(1);
  ht14p1->SetMarkerStyle(kCircle);
  ht14p1->SetMarkerColor(kRed);
  ht14p1->Draw("pe");
  lum->DrawLatex(0.60,0.45,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.40,"Data L3 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.35,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[0],epr[0],epr[1]));
  c2->cd(2);
  ht14p2->SetMarkerSize(1);
  ht14p2->SetMarkerStyle(kCircle);
  ht14p2->SetMarkerColor(kRed);
  ht14p2->Draw("pe");
  lum->DrawLatex(0.60,0.45,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.40,"Data L3 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.35,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[1],epr[1],epr[2]));
  c2->cd(3);
  ht14p3->SetMarkerSize(1);
  ht14p3->SetMarkerStyle(kCircle);
  ht14p3->SetMarkerColor(kRed);
  ht14p3->Draw("pe");
  lum->DrawLatex(0.60,0.45,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.40,"Data L3 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.35,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[2],epr[2],epr[3]));
  c2->cd(4);
  ht14p4->SetMarkerSize(1);
  ht14p4->SetMarkerStyle(kCircle);
  ht14p4->SetMarkerColor(kRed);
  ht14p4->Draw("pe");
  lum->DrawLatex(0.60,0.45,"CMS #it{Internal}");lum->DrawLatex(0.55,0.93,"PbPb #sqrt{s_{NN}} = 5.02 TeV");	
  ds->DrawLatex(0.57, 0.40,"Data L3 Jpsi Trigger Efficiency"); ds->DrawLatex(0.60, 0.35,Form("(p_{T}^{#mu}>%.1f, |#eta| #in [%.1f, %.1f])",ptc[3],epr[3],epr[4]));
  c1->Draw();
  c2->Draw();

  c1->SaveAs("plotEff_MBPD_JpsiL2_v3.pdf");
  c2->SaveAs("plotEff_MBPD_JpsiL3_v3.pdf");
}
