#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "cutsAndBinUpsilonV2.h"
#include "tdrstyle.C"
#include "CMS_lumi_massPull.C"

using namespace RooFit;

using namespace std;

void doFitMBMass(int trigFilter = 3, int etaR = 1, int ptBin=3, int polf= 3)
{
  gStyle->SetEndErrorSize(0);
  TFile *f1 = new TFile(Form("MBPD_Ev6p7M_TrigEff_L%d_v3.root",trigFilter),"read");

  float massLow = 2.7;
  float massHigh = 3.4;
  int nMassBin = (massHigh-massLow)*100;

  int kCheb = 1;
  int kErf = 2;
  int selBkg = kCheb;
  
  vector<vector<double>> ptBins = {
    {3.5, 4, 4.5, 5, 5.5, 6.5, 8., 10.5, 14, 18, 30.},
    {2.07, 3.0, 3.5, 4, 4.5, 5., 6., 7.5, 10, 15, 30},
    {1.5, 2.5, 3, 3.5, 4, 4.5, 5.5, 6.5, 8, 9.5, 13, 20},
    {1.5, 2.2, 2.7, 3.2, 3.7, 4.7, 6.5, 8.5, 11, 20},
  };

  const int nSize = ptBins[etaR-1].size();
  double ptBins_[nSize];
  for(int ib=0;ib<nSize;ib++){ptBins_[ib] = ptBins[etaR-1][ib];}
  
  double ptLow = ptBins_[ptBin-1];
  double ptHigh = ptBins_[ptBin];
  cout << "pt range : " << ptLow << ", " << ptHigh << endl;

  TString kineCut = Form("(etaRange == %d) && ( (pt1>%.2f && pt1<%.2f) || (pt2>%.2f && pt2<%.2f) )", etaR, ptLow, ptHigh, ptLow, ptHigh);

  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);

  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->data("reducedDS")->Print();
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.16, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  c1->SetLeftMargin(2.6);

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  RooRealVar mean("m_{J/#psi}","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi -0.1, pdgMass.JPsi + 0.1 ) ; 

  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass DF",0.01, 0.01, 0.15);

  RooRealVar *x1s = new RooRealVar("x1s","sigma ratio ", 0.5, 0, 1); 

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, *x1s) );

  RooRealVar alpha1s_1("alpha1s_1","tail shift", 3.522 , 1.426, 5.931);
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", 1.905 , 1.36332, 4.345);
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );

  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", 0.5, 0., 1);


  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean, sigma1s_2, alpha1s_2, n1s_2);

  RooAddPdf*  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );

  RooRealVar *nSig= new RooRealVar("nSig"," signals",0,100000);

  double init_mu = 2.6;
  double init_sigma = 0.6;
  double init_lambda = 3;
  RooRealVar err_mu("#mu","err_mu",init_mu,  0, 30) ;
  RooRealVar err_sigma("#sigma","err_sigma", init_sigma, 0,30);
  RooRealVar m_lambda("#lambda","m_lambda",  init_lambda, 0,30);

  RooGenericPdf *bkg;
  RooGenericPdf *bkgErfExp = new RooGenericPdf("bkgErfExp","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );

  RooRealVar ch4_k1("ch4_k1","ch4_k1",0.1,-1.,1.) ;
  RooRealVar ch4_k2("ch4_k2","ch4_k2",-0.1,-1.,1.) ;
  RooRealVar ch4_k3("ch4_k3","ch4_k3",0.1,-1.,1.) ;
  RooRealVar ch4_k4("ch4_k4","ch4_k4",-0.1,-1.,1.) ;
  RooChebychev *bkgCheb1 = new RooChebychev("cPol1Bkg","Background",*(ws->var("mass")),RooArgSet(ch4_k1));
  RooChebychev *bkgCheb2 = new RooChebychev("cPol2Bkg","Background",*(ws->var("mass")),RooArgSet(ch4_k1,ch4_k2));
  RooChebychev *bkgCheb3 = new RooChebychev("cPol3Bkg","Background",*(ws->var("mass")),RooArgSet(ch4_k1,ch4_k2,ch4_k3));
  RooChebychev *bkgCheb4 = new RooChebychev("cPol4Bkg","Background",*(ws->var("mass")),RooArgSet(ch4_k1,ch4_k2,ch4_k3,ch4_k4));

  if(selBkg == kCheb){
    if(polf==1) bkg = (RooGenericPdf*) bkgCheb1;
    else if(polf==2) bkg = (RooGenericPdf*) bkgCheb2;
    else if(polf==3) bkg = (RooGenericPdf*) bkgCheb3;
    else if(polf==4) bkg = (RooGenericPdf*) bkgCheb4;
  }
  else if(selBkg == kErf) bkg = bkgErfExp ;

  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",0,1000000);
  RooAddPdf* model = new RooAddPdf();
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *bkg),RooArgList(*nSig,*nBkg));

  ws->import(*model);


  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Name("sigPDF"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLow, massHigh,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.4);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.048);
  myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();
  fitRes2->Print("v");
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

  float pos_text_x = 0.41;
  float pos_text_y = 0.78;
  float pos_y_diff = 0.056;
  float text_size = 15;
  int text_color = 1;
  /*
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 2.5 && ptHigh==5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S)
  {
    drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
    drawText("|#eta^{#mu}| < 2.4 GeV/c", pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
    drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);
  }
  else {
    drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
    drawText("|#eta^{#mu}| < 2.4 GeV/c", pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  }

  drawText(Form("Signal Function : %s CB", SignalCB.Data() ), 0.55,0.54,1,14);
*/
  TLegend* fitleg = new TLegend(0.68,0.42,0.88,0.7); fitleg->SetTextSize(15);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
  fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  fitleg->AddEntry(myPlot2->findObject("sigPDF"),"Signal","l");
  fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
  fitleg->Draw("same");

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.23);
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.63);
  pad1->SetLeftMargin(0.15);
  pad2->SetLeftMargin(0.15);
  pad2->SetTicks(1,1);
  pad2->cd();

  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.31) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.17) ;
  pullFrame->GetYaxis()->SetLabelSize(0.13) ;
  pullFrame->GetYaxis()->SetRangeUser(-4.5,4.5) ;
  //pullFrame->GetYaxis()->SetLimits(-6,6) ;
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.20) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetLabelSize(0.183) ;
  pullFrame->GetXaxis()->SetTitleSize(0.25) ;
  pullFrame->GetXaxis()->CenterTitle();
  // pullFrame->GetXaxis()->SetTitleFont(43);
  // pullFrame->GetYaxis()->SetTitleFont(43);

  pullFrame->GetYaxis()->SetTickSize(0.02);
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw() ;

  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin;
  double *ypull = hpull->GetY();
  for(int i=0;i<nBins;i++)
  {
    if(ypull[i] == 0) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }

  int numFitPar = fitRes2->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;

  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  pad1->Update();


  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"JPsi");

  float temp1 = ws->var("nSig")->getVal();
  float temp1err = ws->var("nSig")->getError();

  outh->SetBinContent(1,  temp1 ) ;
  outh->SetBinError  (1,  temp1err ) ;

  cout << "signal    =  " << outh->GetBinContent(1) << " +/- " << outh->GetBinError(1) << endl;

  setTDRStyle();
  writeExtraText = true;
  extraText = "Internal";

  TString label;
  label="";
  CMS_lumi_massPull(pad1, 2 ,33);


  pad1->Update();
  pad2->Update();

  c1->cd();
  pad1->Draw();
  pad2->Draw();

  pad1->Update();
  pad2->Update();

  TFile* outf = new TFile(Form("fitRes/fitresults_MBData_DoubleCB_bkgpol%d_ptBin%d_etaRange%d.root",polf,ptBin,etaR),"recreate");
  outh->Write();
  c1->SaveAs(Form("fitRes/fitresults_MBData_DoubleCB_bkgpol%d_ptBin%d_etaRange%d.pdf",polf,ptBin,etaR));
  c1->Write();
  ws->Write();
}
