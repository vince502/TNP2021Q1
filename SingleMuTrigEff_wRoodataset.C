#include <ctime>

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TEfficiency.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <dirent.h>
#include <memory>
#include "cutsAndBinUpsilonV2.h"

using namespace std;
static const long MAXTREESIZE = 1000000000000;
double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);

bool accmatch(double eta, double pt){
  int acase=-1;
  bool tmp = false;
  if( eta > -2.4 && eta < 2.4){
    if (TMath::Abs(eta) < 0.3) acase=1;
    else if (TMath::Abs(eta) >=0.3 && TMath::Abs(eta) <1.1) acase =2;
    else if (TMath::Abs(eta) >=1.1 && TMath::Abs(eta) < 1.4) acase =3;
    else if (TMath::Abs(eta) >=1.4 && TMath::Abs(eta) <1.55) acase =4;
    else if (TMath::Abs(eta) >=1.55 && TMath::Abs(eta) <2.2) acase =5;
    else if (TMath::Abs(eta) >=2.2) acase =6;
  }
  switch(acase){
    case -1: tmp = false;
    case 1: tmp = (bool) (pt > 3.4);
    case 2: tmp = (bool) (pt > 3.3);
    case 3: tmp = (bool) (pt > 7.7-4.0*TMath::Abs(eta));
    case 4: tmp = (bool) (pt > 2.1);
    case 5: tmp = (bool) (pt > 4.25-1.39*TMath::Abs(eta));
    case 6: tmp = (bool) (pt > 1.2);
  }
  return tmp;
};

void SingleMuTrigEff_wRoodataset(int nevt=-1, bool isMC = false, int kTrigSel = 14, int hiHFBinEdge = 0, int PDtype = 1) 
{

  ROOT::EnableImplicitMT();
  // Example of using event plane namespace 

  TString fname1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/MinBias/HIMinimumBias_Run2018_Upsilon_PromptReco_v1.root";
  TString fnameData1 ="../Oniatree_skim_AODSkim_v2_1_part_20201214.root"; //"/u/user/vince402/root_tmpONIATREE_HIMB_20201228.root"; 
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
//
  TString fnameDataReReco = "/home/samba.old/CMS_Files/UpsilonAnalysis/Ups3S_PbPb2018/OniaTree/Data/SkimTree_MB_forUpsilon3SStudy_Data_20210119.root";
  //"/u/user/vince402/root_tmp/ONIATREE_HIMB_20201228.root";// "./Oniatree_part_20201217/Oniatree_skim_AODSkim_v2_1*.root";

  TString fnameDataReRecoPeri = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuonPsiPeri/ReReco_Oniatree_addvn_part*.root";
  TString fnameMC = "/eos/cms/store/group/phys_heavyions/dileptons/MC2018/PbPb502TeV/TTrees/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC*_v20190801.root";

  TString fPD;
  if(PDtype==1) fPD = "DB";
  TChain *mytree = new TChain("hionia/myTree");
  if(!isMC){
  if(PDtype==1) mytree->Add(fnameDataReReco.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Float_t         SumET_HF;
  Short_t           Reco_QQ_size;
  Short_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_SumET_HF;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Short_t           Reco_QQ_mupl_idx[maxBranchSize];
  Short_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);

  Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!

//  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Short_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_ptErr_global;   //!
//  mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Float_t         Reco_QQ_ctau3D[maxBranchSize];
  Float_t         Reco_QQ_ctauErr3D[maxBranchSize];
  TBranch        *b_Reco_QQ_ctau3D;
  TBranch        *b_Reco_QQ_ctauErr3D;
  mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  mytree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
  
  Int_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  if(isMC){
    mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
    mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);
  }

  
  
  const int nEP = 29;  // number of event planes in the tree
  int trigIndx=0;
  if(kTrigSel == kTrigJpsi) trigIndx=0;
  else if(kTrigSel == kTrigUps) trigIndx=1;
  else if(kTrigSel == kTrigL1DBOS40100) trigIndx=2;
  else if(kTrigSel == kTrigL1DB50100) trigIndx=3;
  else if(kTrigSel == 0) trigIndx=4;
  
  double tnp_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;

  int kL2filter = 38;
  int kL3filter = 39;

  int count =0;
  int counttnp =0;

  TString fCentSelHF = "HFNom";
  if(hiHFBinEdge==1) fCentSelHF = "HFUp";
  else if(hiHFBinEdge==-1) fCentSelHF = "HFDown";
  TFile* newfile;
  
  newfile = new TFile(Form("MBPD_Ev6p7M_TrigEff_L%d_v3.root",(int) (kTrigSel-11)),"recreate");

  int cBin;

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << endl;
  float massLow = 2.7;
  float massHigh = 3.4;

  double massLow_sig  = 3.0;
  double massHigh_sig = 3.2;

  double massLow_bkgh  = 3.2;
  double massHigh_bkgh = 3.4;
  double massLow_bkgl  = 2.7;
  double massHigh_bkgl = 3.0;

  const int nMassBin = (massHigh - massLow)/0.01;
  cout << nMassBin << endl;

  const int nEta=4;
  const int nPt[nEta] = {10, 10, 11, 9};
  double etaRange[nEta+1] = {0,1.2,1.8,2.1,2.4};

  vector<vector<double>> ptBins = {
    {3.5, 4, 4.5, 5, 5.5, 6.5, 8., 10.5, 14, 18, 30.},
    {2.07, 3.0, 3.5, 4, 4.5, 5., 6., 7.5, 10, 15, 30},
    {1.5, 2.5, 3, 3.5, 4, 4.5, 5.5, 6.5, 8, 9.5, 13, 20},
    {1.5, 2.2, 2.7, 3.2, 3.7, 4.7, 6.5, 8.5, 11, 20},
  };
  
  map<TString,TH1D*> hDen;
  map<TString,TH1D*> hNum;
  map<TString,TH1D*> hMass;

  map<TString,TH1D*> hDen_bkgl;
  map<TString,TH1D*> hDen_bkgh;
  map<TString,TH1D*> hNum_bkgl;
  map<TString,TH1D*> hNum_bkgh;
  map<TString,TH1D*> hDen_sig;
  map<TString,TH1D*> hNum_sig;
  map<TString,TH1D*> hMass_bkg;

  for(int ieta=0; ieta<nEta; ieta++){
    const int nSize = ptBins[ieta].size();
    double ptBins_[nSize];
    for(int ib=0;ib<nSize;ib++){ptBins_[ib] = ptBins[ieta][ib];}
    hDen[ieta] = new TH1D(Form("hDen_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hNum[ieta] = new TH1D(Form("hNum_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hDen[ieta] -> Sumw2();
    hNum[ieta] -> Sumw2();
    
    hDen_bkgl[ieta] = new TH1D(Form("hDen_bkgl_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hNum_bkgl[ieta] = new TH1D(Form("hNum_bkgl_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hDen_bkgh[ieta] = new TH1D(Form("hDen_bkgh_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hNum_bkgh[ieta] = new TH1D(Form("hNum_bkgh_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hDen_sig[ieta] = new TH1D(Form("hDen_sig_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hNum_sig[ieta] = new TH1D(Form("hNum_sig_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hDen_bkgl[ieta] -> Sumw2();
    hDen_bkgh[ieta] -> Sumw2();
    hDen_sig[ieta] -> Sumw2();
    hNum_bkgl[ieta] -> Sumw2();
    hNum_bkgh[ieta] -> Sumw2();
    hNum_sig[ieta] -> Sumw2();

    hMass[ieta] = new TH1D(Form("hMass_eta%d",ieta),";m_{#mu^{+}#mu^{-}};",nMassBin,massLow,massHigh);
    hMass_bkg[ieta] = new TH1D(Form("hMass_bkg_eta%d",ieta),";m_{#mu^{+}#mu^{-}};",nMassBin,massLow,massHigh);
  }

  Double_t epr[nEta+1] = {0.0, 1.2, 1.8, 2.1, 2.4};
  Double_t ptc[4] = {3.5, 2.1, 1.5, 1.5};


  //For Roodataset 
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* pt1Var    = new RooRealVar("pt1","pt1 variable", 0,100,"GeV/c");
  RooRealVar* pt2Var    = new RooRealVar("pt2","pt2 variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* etaRangeVar     = new RooRealVar("etaRange","dimuon mass, single muon eta region", -10,10,"");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var,*cBinVar, *etaRangeVar);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);


  for(int iev=0; iev<nevt ; ++iev)
  {

    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    bool jpsiEvent = false;
    cBin = getHiBinFromhiHF(SumET_HF);
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq){
      if (!(Reco_QQ_sign[irqq]==0)) continue; 
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      double mupt1 = mupl_Reco->Pt();
      double mupt2 = mumi_Reco->Pt();
      double mueta1 = mupl_Reco->Eta();
      double mueta2 = mumi_Reco->Eta();
      double mass_dimu = JP_Reco->M();

      bool isTrigMatchMuPl = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2,kTrigSel)));
      bool isTrigMatchMuMi = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2,kTrigSel)));
      bool isMassSig = (mass_dimu > massLow_sig && mass_dimu < massHigh_sig);
      bool isMassBkgl = (mass_dimu > massLow_bkgl && mass_dimu < massHigh_bkgl);
      bool isMassBkgh = (mass_dimu > massLow_bkgh && mass_dimu < massHigh_bkgh);

      if( Reco_QQ_VtxProb[irqq] <= 0.01) continue;
      if (!(mass_dimu < massHigh && mass_dimu > massLow)) continue;//v2 : 2.7~3.3, v3: 3~3.19
      if(!(cBin <180)) continue;
      if(! (accmatch(mueta1, mupt1) && accmatch(mueta2, mupt2)) ) continue;
      etaRangeVar->setVal( (int) -1);
      for(int ieta=0; ieta<nEta; ieta++){
        if( abs(mueta1) >= epr[ieta] && abs(mueta1) < epr[ieta+1] ){
          hDen[ieta]->Fill(mupt1);
          if(isTrigMatchMuPl) hNum[ieta]->Fill(mupt1);

          if(isMassSig){
            hDen_sig[ieta]->Fill(mupt1);
            if(isTrigMatchMuPl) hNum_sig[ieta]->Fill(mupt1);
          }
          if(isMassBkgl){
            hDen_bkgl[ieta]->Fill(mupt1);
            if(isTrigMatchMuPl) hNum_bkgl[ieta]->Fill(mupt1);
          }
          if(isMassBkgh){
            hDen_bkgh[ieta]->Fill(mupt1);
            if(isTrigMatchMuPl) hNum_bkgh[ieta]->Fill(mupt1);
          }
        }
        if( abs(mueta2) >= epr[ieta] && abs(mueta2) < epr[ieta+1] ){
          hDen[ieta]->Fill(mupt2);
          if(isTrigMatchMuMi) hNum[ieta]->Fill(mupt2);

          if(isMassSig){
            hDen_sig[ieta]->Fill(mupt2);
            if(isTrigMatchMuMi) hNum_sig[ieta]->Fill(mupt2);
          }
          if(isMassBkgl){
            hDen_bkgl[ieta]->Fill(mupt2);
            if(isTrigMatchMuMi) hNum_bkgl[ieta]->Fill(mupt2);
          }
          if(isMassBkgh){
            hDen_bkgh[ieta]->Fill(mupt2);
            if(isTrigMatchMuMi) hNum_bkgh[ieta]->Fill(mupt2);
          }
        }
        
        if( (abs(mueta1) >= epr[ieta] && abs(mueta1) < epr[ieta+1]) || (abs(mueta2) >= epr[ieta] && abs(mueta2) < epr[ieta+1]) ){
          hMass[ieta]->Fill(mass_dimu);
          if(mass_dimu > massHigh_sig || mass_dimu < massLow_sig) hMass_bkg[ieta]->Fill(mass_dimu);
          etaRangeVar->setVal( (int) ieta+1);
        }
      }
      
      //Fill roodataset
      if(etaRangeVar->getVal() == -1){cout << "no matched eta region!!! " << endl; return;}
      massVar -> setVal( (double) mass_dimu );
      ptVar   -> setVal( (double) JP_Reco->Pt());
      pt1Var   -> setVal( (double) mupt1);
      pt2Var   -> setVal( (double) mupt2);
      yVar    -> setVal( (double) JP_Reco->Rapidity());
      cBinVar -> setVal( (int) cBin);
      dataSet->add( *argSet);
    }
  } //end of event loop
  //  mmtree->Write();  // Don't need to call Write() for trees 
  
  map<TString,TH1D*> heff;
  map<TString,TH1D*> heff_bkgl;
  map<TString,TH1D*> heff_bkgh;
  map<TString,TH1D*> heff_sig;
  newfile->cd();
  for(int ieta=0; ieta<nEta; ieta++){
    heff[ieta] = (TH1D*) hNum[ieta] -> Clone(Form("hEff_eta%d",ieta));
    heff_bkgl[ieta] = (TH1D*) hNum_bkgl[ieta] -> Clone(Form("hEff_bkgl_eta%d",ieta));
    heff_bkgh[ieta] = (TH1D*) hNum_bkgh[ieta] -> Clone(Form("hEff_bkgh_eta%d",ieta));
    heff_sig[ieta] = (TH1D*) hNum_sig[ieta] -> Clone(Form("hEff_sig_eta%d",ieta));
    heff[ieta]->Divide(hDen[ieta]);
    heff_bkgl[ieta]->Divide(hDen_bkgl[ieta]);
    heff_bkgh[ieta]->Divide(hDen_bkgh[ieta]);
    heff_sig[ieta]->Divide(hDen_sig[ieta]);
    
    hNum[ieta]->Write();
    hDen[ieta]->Write();
    heff[ieta]->Write();

    hDen_bkgl[ieta]->Write();
    heff_bkgl[ieta]->Write();
    hNum_bkgh[ieta]->Write();
    hDen_bkgh[ieta]->Write();
    heff_bkgh[ieta]->Write();
    hNum_sig[ieta]->Write();
    hDen_sig[ieta]->Write();
    heff_sig[ieta]->Write();
    hMass[ieta]->Write();
    hMass_bkg[ieta]->Write();
  }
  dataSet->Write();
  newfile->Close();

} 
