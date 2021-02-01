#include <iostream>
#include "Style_jaebeom.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "cutsAndBinUpsilonV2.h"
#include "tdrstyle.C"
#include "CMS_lumi_massPull.C"


using namespace std;
using namespace RooFit;

void calcSigEff(int trigFilter = 3, int bkgpol=3)
{
  const int nEta = 4;

  const char* bkgstr[4] = {"cPol1Bkg","cPol2Bkg","cPol3Bkg","cPol4Bkg"};

  TFile* fEff = new TFile(Form("MBPD_Ev6p7M_TrigEff_L%d_v3.root",trigFilter),"read");
  TH1D* hEff_bkgh[nEta];
  TH1D* hEff_bkgl[nEta];
  TH1D* hEff_sig[nEta];

  TGraphAsymmErrors* gEff_sig[nEta];
  TGraphAsymmErrors* gEff_bkg[nEta];
  TGraphAsymmErrors* gEff_sig_bkgsub[nEta];

  int nPtBins[nEta] = {10,10,11,9};

  TFile* ffit[nEta];
  RooWorkspace* ws[nEta];
  RooAbsPdf *pdf_sig[nEta];
  RooAbsPdf *pdf_bkg[nEta];
  RooAbsPdf *pdf_tot[nEta];
  RooRealVar* massVar[nEta];
  double nSig[nEta];
  double nBkg[nEta];

  double xtot_min = 2.7;
  double xtot_max = 3.4;
  double xsig_min = 3.0;
  double xsig_max = 3.2;
  double xbkgl_min = 2.7;
  double xbkgl_max = 3.0;
  double xbkgh_min = 3.2;
  double xbkgh_max = 3.4;

  RooAbsReal* integral_tot[nEta];
  RooAbsReal* integral_tot_sig[nEta];
  RooAbsReal* integral_bkg_sig[nEta];
  RooAbsReal* integral_bkgl[nEta];
  RooAbsReal* integral_bkgh[nEta];
  RooAbsReal* integral_sig[nEta];

  TH1D* heff_bkg[nEta];
  TH1D* heff_sig_bkgsub[nEta];

  for(int i=0;i<nEta;i++){
    //histogram Eff 
    hEff_bkgh[i] = (TH1D*) fEff->Get(Form("hEff_bkgh_eta%d",i));
    hEff_bkgl[i] = (TH1D*) fEff->Get(Form("hEff_bkgl_eta%d",i));
    hEff_sig[i] = (TH1D*) fEff->Get(Form("hEff_sig_eta%d",i));

    gEff_sig[i] = new TGraphAsymmErrors();
    gEff_bkg[i] = new TGraphAsymmErrors();
    gEff_sig_bkgsub[i] = new TGraphAsymmErrors();
    gEff_sig[i] -> SetName(Form("gEff_sig_eta%d",i));
    gEff_bkg[i] -> SetName(Form("gEff_bkg_eta%d",i));
    gEff_sig_bkgsub[i] -> SetName(Form("gEff_sig_bkgsub_eta%d",i));

    heff_bkg[i] = (TH1D*) hEff_bkgl[i] -> Clone(Form("hEff_bkg_avg_eta%d",i));
    heff_bkg[i] -> Reset();
    heff_sig_bkgsub[i] = (TH1D*) hEff_sig[i] -> Clone(Form("hEff_sig_bkgsub_eta%d",i));
    heff_sig_bkgsub[i]->Reset();

    //FitResult
    if(nPtBins[i] != hEff_bkgh[i]->GetNbinsX()){cout << "Inconsistent # of bins :: " << endl; return;}

    for(int ibin=1; ibin<=nPtBins[i]; ibin++){
      
      TFile* fileFit = new TFile(Form("fitRes/fitresults_MBData_DoubleCB_bkgpol%d_ptBin%d_etaRange%d.root",bkgpol,ibin,i+1),"read");
      RooWorkspace *ws = (RooWorkspace*) fileFit->Get("workspace");
      RooAbsPdf *pdf_sig =  (RooAbsPdf*) ws->pdf("cb1s");  
      RooAbsPdf *pdf_bkg =  (RooAbsPdf*) ws->pdf(Form("%s",bkgstr[bkgpol-1])); 
      RooAbsPdf *pdf_tot =  (RooAbsPdf*) ws->pdf("model"); 
      RooRealVar *massVar = ws->var("mass");
      massVar->setRange("all",xtot_min,xtot_max);
      massVar->setRange("sig",xsig_min,xsig_max);
      massVar->setRange("bkgl",xbkgl_min,xbkgl_max);
      massVar->setRange("bkgh",xbkgh_min,xbkgh_max);
      double nSig = ws->var("nSig")->getVal();
      double nBkg = ws->var("nBkg")->getVal();
      RooAbsReal* integral_tot_sig = pdf_tot->createIntegral(*massVar,RooFit::NormSet(*massVar),RooFit::Range("sig")); 
      RooAbsReal* integral_bkg_sig= pdf_bkg->createIntegral(*massVar,RooFit::NormSet(*massVar),RooFit::Range("sig"));
      RooAbsReal* integral_bkgl= pdf_bkg->createIntegral(*massVar,RooFit::NormSet(*massVar),RooFit::Range("bkgl"));
      RooAbsReal* integral_bkgh= pdf_bkg->createIntegral(*massVar,RooFit::NormSet(*massVar),RooFit::Range("bkgh"));
      RooAbsReal* integral_sig= pdf_sig->createIntegral(*massVar,RooFit::NormSet(*massVar),RooFit::Range("sig"));

      double nTot = nSig+nBkg;

      //Calculate Eff bkg 
      double bkg_sum =  nBkg * ((double) integral_bkgl->getVal() + (double) integral_bkgh->getVal());
      double scale_bkgl = nBkg * (double) integral_bkgl->getVal() / bkg_sum;
      double scale_bkgh = nBkg * (double) integral_bkgh->getVal() / bkg_sum;

      double bkgContent = scale_bkgl * hEff_bkgl[i]->GetBinContent(ibin) + scale_bkgh * hEff_bkgh[i]->GetBinContent(ibin);
      double bkgErrorl   = scale_bkgl * hEff_bkgl[i]->GetBinError(ibin);
      double bkgErrorh   = scale_bkgh * hEff_bkgh[i]->GetBinError(ibin);
      double bkgError = sqrt(bkgErrorl*bkgErrorl + bkgErrorh*bkgErrorh);
      heff_bkg[i]->SetBinContent(ibin,bkgContent);
      heff_bkg[i]->SetBinError(ibin,bkgError);
      
      double sig_sum = nSig * (double) integral_sig->getVal();
      double tot_sigsum = nTot * (double) integral_tot_sig->getVal();
      double bkg_sigsum = nBkg * (double) integral_bkg_sig->getVal();

      double sigContent = (tot_sigsum/sig_sum) * hEff_sig[i]->GetBinContent(ibin) - (bkg_sigsum/sig_sum) * heff_bkg[i]->GetBinContent(ibin);
      double sigError_s = (tot_sigsum/sig_sum) * hEff_sig[i]->GetBinError(ibin);
      double sigError_b =  (bkg_sigsum/sig_sum) * heff_bkg[i]->GetBinError(ibin);
      double sigError = sqrt(sigError_s*sigError_s + sigError_b*sigError_b);
      heff_sig_bkgsub[i]->SetBinContent(ibin,sigContent);
      heff_sig_bkgsub[i]->SetBinError(ibin,sigError);

      double binCent = heff_sig_bkgsub[i]->GetXaxis()->GetBinCenter(ibin);

      gEff_sig[i]->SetPoint(ibin-1, binCent,hEff_sig[i]->GetBinContent(ibin));
      gEff_sig[i]->SetPointError(ibin-1, 0,0,hEff_sig[i]->GetBinError(ibin), hEff_sig[i]->GetBinError(ibin));
      gEff_bkg[i]->SetPoint(ibin-1, binCent,heff_bkg[i]->GetBinContent(ibin));
      gEff_bkg[i]->SetPointError(ibin-1, 0, 0, heff_bkg[i]->GetBinError(ibin), heff_bkg[i]->GetBinError(ibin));
      gEff_sig_bkgsub[i]->SetPoint(ibin-1, binCent,heff_sig_bkgsub[i]->GetBinContent(ibin));
      gEff_sig_bkgsub[i]->SetPointError(ibin-1, 0, 0, heff_sig_bkgsub[i]->GetBinError(ibin), heff_sig_bkgsub[i]->GetBinError(ibin));
    }
  }

  TFile* fEffwr = new TFile(Form("Output_MBPD_Ev6p7M_TrigEff_L%d_v3.root",trigFilter),"recreate");
  fEffwr->cd();
  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->Divide(2,2,0,0);

  for(int ieta=0;ieta<nEta;ieta++){
    gEff_sig[ieta]->SetMarkerStyle(kFullTriangleUp);
    gEff_sig[ieta]->SetMarkerSize(1);
    gEff_sig[ieta]->SetMarkerColor(kBlack);
    gEff_sig[ieta]->SetLineColor(kBlack);
    gEff_sig[ieta]->GetXaxis()->SetLimits(0,30);
    gEff_sig[ieta]->GetXaxis()->SetRangeUser(0,30);
    gEff_sig[ieta]->GetYaxis()->SetLimits(0,1);
    gEff_sig[ieta]->GetYaxis()->SetRangeUser(0,1);
    c1->cd(ieta+1);
    gEff_sig[ieta]->Draw("AP");
    c1->Modified();
    c1->Update();
    gEff_bkg[ieta]->SetMarkerStyle(kFullTriangleDown);
    gEff_bkg[ieta]->SetMarkerSize(1);
    gEff_bkg[ieta]->SetMarkerColor(kBlue+1);
    gEff_bkg[ieta]->SetLineColor(kBlue+1);
    gEff_bkg[ieta]->Draw("pe same");
    gEff_bkg[ieta]->GetYaxis()->SetRangeUser(0,1);
    gEff_bkg[ieta]->GetXaxis()->SetLimits(0,30);
    gEff_bkg[ieta]->GetXaxis()->SetRangeUser(0,30);
    gEff_sig_bkgsub[ieta]->GetYaxis()->SetRangeUser(0,1);
    gEff_sig_bkgsub[ieta]->GetXaxis()->SetLimits(0,30);
    gEff_sig_bkgsub[ieta]->GetXaxis()->SetRangeUser(0,30);
    gEff_sig_bkgsub[ieta]->SetMarkerStyle(kFullCircle);
    gEff_sig_bkgsub[ieta]->SetMarkerSize(1);
    gEff_sig_bkgsub[ieta]->SetMarkerColor(kRed+1);
    gEff_sig_bkgsub[ieta]->SetLineColor(kRed+1);
    gEff_sig_bkgsub[ieta]->Draw("pe same");
    c1->Modified();
    c1->Update();
    c1->Modified();
    c1->Update();
    gEff_sig_bkgsub[ieta]->Write();
    gEff_sig[ieta]->Write();
    gEff_bkg[ieta]->Write();
  }
  TLegend* leg = new TLegend(0.43,0.12,0.76,0.40);
  SetLegendStyle(leg);
  leg->SetTextSize(0.042);
  leg->AddEntry(gEff_sig[0],"Signal region","pe");
  leg->AddEntry(gEff_bkg[0],"Sideband region","pe");
  leg->AddEntry(gEff_sig_bkgsub[0],"Sig. background subtracted","pe");
  c1->cd(1);
  leg->Draw("same");
  c1->Update();
  c1->SaveAs("MB_Eff.pdf");
  c1->Write();

}
