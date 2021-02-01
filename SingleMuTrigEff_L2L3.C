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

void SingleMuTrigEff_L2L3(int nevt=-1, bool isMC = false, int hiHFBinEdge = 0, int PDtype = 1) 
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

  TString fCentSelHF = "HFNom";
  if(hiHFBinEdge==1) fCentSelHF = "HFUp";
  else if(hiHFBinEdge==-1) fCentSelHF = "HFDown";
  TFile* newfile;
  
  newfile = new TFile("MB_L2L3_Ev6p7M_TrigEff_v3.root","recreate");

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
  float massLow = 3.0;
  float massHigh = 3.2;

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
  
  int kTrigSelL3 = 14;
  int kTrigSelL2 = 13;

  for(int ieta=0; ieta<nEta; ieta++){
    const int nSize = ptBins[ieta].size();
    double ptBins_[nSize];
    for(int ib=0;ib<nSize;ib++){ptBins_[ib] = ptBins[ieta][ib];}
    hDen[ieta] = new TH1D(Form("hDen_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hNum[ieta] = new TH1D(Form("hNum_eta%d",ieta),";p_{T};",nPt[ieta],ptBins_);
    hDen[ieta] -> Sumw2();
    hNum[ieta] -> Sumw2();
  }

  Double_t epr[nEta+1] = {0.0, 1.2, 1.8, 2.1, 2.4};
  Double_t ptc[4] = {3.5, 2.1, 1.5, 1.5};

  for(int iev=0; iev<nevt ; ++iev)
  {

    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
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

      bool isTrigMatchMuPlL2 = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSelL2))) == ((ULong64_t)pow(2,kTrigSelL2)));
      bool isTrigMatchMuMiL2 = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2,kTrigSelL2))) == ((ULong64_t)pow(2,kTrigSelL2)));
      bool isTrigMatchMuPlL3 = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSelL3))) == ((ULong64_t)pow(2,kTrigSelL3)));
      bool isTrigMatchMuMiL3 = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2,kTrigSelL3))) == ((ULong64_t)pow(2,kTrigSelL3)));

      if( Reco_QQ_VtxProb[irqq] <= 0.01) continue;
      if (!(mass_dimu < massHigh && mass_dimu > massLow)) continue;//v2 : 2.7~3.3, v3: 3~3.19
      if(!(cBin <180)) continue;
      if(! (accmatch(mueta1, mupt1) && accmatch(mueta2, mupt2)) ) continue;

      for(int ieta=0; ieta<nEta; ieta++){
        if( abs(mueta1) >= epr[ieta] && abs(mueta1) < epr[ieta+1] ){
          if(isTrigMatchMuPlL2){
            hDen[ieta]->Fill(mupt1);
            if(isTrigMatchMuPlL3) hNum[ieta]->Fill(mupt1);
          }
        }
        if( abs(mueta2) >= epr[ieta] && abs(mueta2) < epr[ieta+1] ){
          if(isTrigMatchMuMiL2){
            hDen[ieta]->Fill(mupt2);
            if(isTrigMatchMuMiL3) hNum[ieta]->Fill(mupt2);
          }
        }
      }
    }
  } //end of event loop
  //  mmtree->Write();  // Don't need to call Write() for trees 
  
  map<TString,TH1D*> heff;
  newfile->cd();
  for(int ieta=0; ieta<nEta; ieta++){
    heff[ieta] = (TH1D*) hNum[ieta] -> Clone(Form("hEff_eta%d",ieta));
    heff[ieta]->Divide(hDen[ieta]);
    hNum[ieta]->Write();
    hDen[ieta]->Write();
    heff[ieta]->Write();
  }
  newfile->Close();
} 
