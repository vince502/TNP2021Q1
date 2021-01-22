#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"
#include "tnp_weight_lowptPbPb.h"

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

void SingleMuTrigEff(int nevt=-1, bool isMC = false, int kTrigSel = 14, int hiHFBinEdge = 0, int PDtype = 1) 
{

  using namespace std;
  using namespace hi;
  ROOT::EnableImplicitMT();
  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  TString fname1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/MinBias/HIMinimumBias_Run2018_Upsilon_PromptReco_v1.root";
  TString fnameData1 ="../Oniatree_skim_AODSkim_v2_1_part_20201214.root"; //"/u/user/vince402/root_tmpONIATREE_HIMB_20201228.root"; 
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
//
  TString fnameDataReReco = "/u/user/vince402/root_tmp/SkimTree_MB_forUpsilon3SStudy_Data_20210119.root"; //* Storage: T2_KNU_KR
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
  
  newfile = new TFile(Form("MBPD_Ev6p7M_TrigEff_L%d_v2.root",(int) (kTrigSel-11)),"recreate");

 // const static int nMaxDimu = 1000;
 // int evt;
 // int runN;
 // int lumi;
  int cBin;
 // int nDimu;
 // float vz;
 // float mass[nMaxDimu];
 // float pt[nMaxDimu];
 // float y[nMaxDimu];
 // float phi[nMaxDimu];
 // float eta[nMaxDimu];
 // float eta1[nMaxDimu];
 // float eta2[nMaxDimu];
 // float phi1[nMaxDimu];
 // float phi2[nMaxDimu];
 // float pt1[nMaxDimu];
 // float pt2[nMaxDimu];
 // float weight0[nMaxDimu]; 
 // float weight1[nMaxDimu]; 
 // float qxa[nMaxDimu];
 // float qya[nMaxDimu];
 // float qxb[nMaxDimu];
 // float qyb[nMaxDimu];
 // float qxc[nMaxDimu];
 // float qyc[nMaxDimu];
 // float qxdimu[nMaxDimu];
 // float qydimu[nMaxDimu];
 // float qxmupl[nMaxDimu];
 // float qxmumi[nMaxDimu];
 // float qymupl[nMaxDimu];
 // float qymumi[nMaxDimu];
 // int recoQQsign[nMaxDimu];
 // float ctau3D[nMaxDimu];
 // float ctau3DErr[nMaxDimu];
 // double TnPweight[nMaxDimu] = {1.};
 // double weight = 1;

 // TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
 // mmevttree->SetMaxTreeSize(MAXTREESIZE);
 // mmevttree->Branch("event",&evt,"event/I");
 // mmevttree->Branch("runN",&runN,"runN/I");
 // mmevttree->Branch("lumi",&lumi,"lumi/I");
 // mmevttree->Branch("cBin",&cBin,"cBin/I");
 // mmevttree->Branch("vz",&vz,"vz/F");
 // mmevttree->Branch("nDimu",&nDimu,"nDimu/I");
 // mmevttree->Branch("mass",mass,"mass[nDimu]/F");
 // mmevttree->Branch("y",y,"y[nDimu]/F");
 // mmevttree->Branch("pt",pt,"pt[nDimu]/F");
 // mmevttree->Branch("pt1",pt1,"pt1[nDimu]/F");
 // mmevttree->Branch("pt2",pt2,"pt2[nDimu]/F");
 // mmevttree->Branch("eta",eta,"eta[nDimu]/F");
 // mmevttree->Branch("eta1",eta1,"eta1[nDimu]/F");
 // mmevttree->Branch("eta2",eta2,"eta2[nDimu]/F");
 // mmevttree->Branch("qxa",qxa,"qxa[nDimu]/F");
 // mmevttree->Branch("qxb",qxb,"qxb[nDimu]/F");
 // mmevttree->Branch("qxc",qxc,"qxc[nDimu]/F");
 // mmevttree->Branch("qya",qya,"qya[nDimu]/F");
 // mmevttree->Branch("qyb",qyb,"qyb[nDimu]/F");
 // mmevttree->Branch("qyc",qyc,"qyc[nDimu]/F");
 // mmevttree->Branch("qxdimu",qxdimu,"qxdimu[nDimu]/F");
 // mmevttree->Branch("qydimu",qydimu,"qydimu[nDimu]/F");
 // mmevttree->Branch("qxmupl",qxmupl,"qxmupl[nDimu]/F");
 // mmevttree->Branch("qxmumi",qxmumi,"qxmumi[nDimu]/F");
 // mmevttree->Branch("qymupl",qymupl,"qymupl[nDimu]/F");
 // mmevttree->Branch("qymumi",qymumi,"qymumi[nDimu]/F");
 // mmevttree->Branch("recoQQsign",recoQQsign,"recoQQsign[nDimu]/I");
 // mmevttree->Branch("ctau3D",ctau3D,"ctau3D[nDimu]/F");
 // mmevttree->Branch("ctau3DErr",ctau3DErr,"ctau3DErr[nDimu]/F");
 // mmevttree->Branch("weight",&weight,"weight/D");
 // mmevttree->Branch("TnPweight",TnPweight,"TnPweight[nDimu]/D");
      


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << endl;
  TH1D* onebin = new TH1D("one", "onebin", 4,1,4);
  TH1D* md = new TH1D("md", "massdist",60, 2.7,3.3);
  int den1=0, den2=0,  den3=0, den4=0;
  int num1=0, num2=0, num3=0, num4=0;
  
  Double_t BINMAP1[11] = {3.5, 4, 4.5, 5, 5.5, 6.5, 8., 10.5, 14, 18, 30};
  Double_t BINMAP2[11] = {2.07, 3.0, 3.5, 4, 4.5, 5., 6., 7.5, 10, 15, 30};
  Double_t BINMAP3[12] = {1.5, 2.5, 3, 3.5, 4, 4.5, 5.5, 6.5, 8, 9.5, 13, 20};
  Double_t BINMAP4[10] = {1.5, 2.2, 2.7, 3.2, 3.7, 4.7, 6.5, 8.5, 11, 20};

  TH1D* h11 = new TH1D("h11","Denom;counts",10, BINMAP1);
  h11->Sumw2();
  TH1D* h21 = new TH1D("h21","Num;counts", 10, BINMAP1);
  h21->Sumw2();
  TH1D* h12 = new TH1D("h12","Denom;counts",10, BINMAP2);
  h12->Sumw2();
  TH1D* h22 = new TH1D("h22","Num;counts", 10, BINMAP2);
  h22->Sumw2();
  TH1D* h13 = new TH1D("h13","Denom;counts",11, BINMAP3);
  h13->Sumw2();
  TH1D* h23 = new TH1D("h23","Num;counts", 11, BINMAP3);
  h23->Sumw2();
  TH1D* h14 = new TH1D("h14","Denom;counts",9, BINMAP4);
  h14->Sumw2();
  TH1D* h24 = new TH1D("h24","Num;counts", 9, BINMAP4);
  h24->Sumw2();

  Double_t epr[5] = {0.0, 1.2, 1.8, 2.1, 2.4};
  Double_t ptc[4] = {3.5, 2.1, 1.5, 1.5};

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
      float mupt1 = mupl_Reco->Pt();
      float mupt2 = mumi_Reco->Pt();
      float mueta1 = mupl_Reco->Eta();
      float mueta2 = mumi_Reco->Eta();
      float mass_dimu = JP_Reco->M();

      if( Reco_QQ_VtxProb[irqq] <= 0.01) continue;
      if (!(mass_dimu < 3.3 && mass_dimu > 2.7)) continue;
      if(!(cBin <180)) continue;
      if(accmatch(mueta1, mupt1)&& accmatch(mueta2, mupt2)){
	if(( (mupt1 > ptc[0] )&&(TMath::Abs(mueta1)<=epr[1] && TMath::Abs(mueta1)>=epr[0]) )){
	h11->Fill(mupt1); den1++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h21->Fill(mupt1); num1++;}
	}
	if(( (mupt1 > ptc[1] )&&(TMath::Abs(mueta1)<=epr[2] && TMath::Abs(mueta1)>= epr[1]) )){
	h12->Fill(mupt1); den2++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h22->Fill(mupt1); num2++;}
	}
	if(( (mupt1 > ptc[2] )&&(TMath::Abs(mueta1)<= epr[3] && TMath::Abs(mueta1)>=epr[2]) )){
	h13->Fill(mupt1); den3++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h23->Fill(mupt1); num3++;}
	}
	if(( (mupt1 > ptc[3] )&&(TMath::Abs(mueta1)<= epr[4] && TMath::Abs(mueta1)>= epr[3] ) )){
	h14->Fill(mupt1); den4++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h24->Fill(mupt1); num4++;}
	}
	if(( (mupt2 > ptc[0] )&&(TMath::Abs(mueta2)<=epr[1] && TMath::Abs(mueta2)>=epr[0]) )){
	h11->Fill(mupt2); den1++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h21->Fill(mupt2); num1++;}
	}
	if(( (mupt2 > ptc[1] )&&(TMath::Abs(mueta2)<=epr[2] && TMath::Abs(mueta2)>= epr[1]) )){
	h12->Fill(mupt2); den2++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h22->Fill(mupt2); num2++; }
	}
	if(( (mupt2 > ptc[2] )&&(TMath::Abs(mueta2)<= epr[3] && TMath::Abs(mueta2)>=epr[2]) )){
	h13->Fill(mupt2); den3++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h23->Fill(mupt2); num3++; }
	}
	if(( (mupt2 > ptc[3] )&&(TMath::Abs(mueta2)<= epr[4] && TMath::Abs(mueta2)>= epr[3] ) )){
	h14->Fill(mupt2); den4++;
	if(((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2,kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)))){ h24->Fill(mupt2); num4++;}
	}
        md->Fill(mass_dimu);
      }
    }
  } //end of event loop
//  mmtree->Write();  // Don't need to call Write() for trees 
  double div1, div2, div3, div4;
  div1= ((double)num1/(double)den1);
  div2= ((double)num2/(double)den2);
  div3= ((double)num3/(double)den3);
  div4= ((double)num4/(double)den4);
  TH1D* hc1 = (TH1D*) h21->Clone("hc1");
  TH1D* hc2 = (TH1D*) h22->Clone("hc2");
  TH1D* hc3 = (TH1D*) h23->Clone("hc3");
  TH1D* hc4 = (TH1D*) h24->Clone("hc4");
  onebin->SetBinContent(1,div1);
  onebin->SetBinContent(2,div2);
  onebin->SetBinContent(3,div3);
  onebin->SetBinContent(4,div4);
  hc1->Divide(h11); 
  hc2->Divide(h12); 
  hc3->Divide(h13); 
  hc4->Divide(h14); 
  newfile->cd();
  hc1->SetName("hEff1");
  hc2->SetName("hEff2");
  hc3->SetName("hEff3");
  hc4->SetName("hEff4");
  hc1->Write();/* h11->Write();h21->Write();*/
  hc2->Write();/* h12->Write();h22->Write();*/
  hc3->Write();/* h13->Write();h23->Write();*/
  hc4->Write();/* h14->Write();h24->Write();*/
  md->Write();
  onebin->Write();
  newfile->Close();
    
} 
