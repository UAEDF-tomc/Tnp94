#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TMVA/Reader.h"
#include <algorithm>

void addMVAmuon(TString fileName, TString year){
  TFile *f = TFile::Open(fileName + ".root");
  TTree *tIn = (TTree *) f->Get("tpTree/fitter_tree");
  //// Variables for MVA 
  int medium;
  float pt, eta, sip3d, logdxy, logdz, segmComp, trackMult, miniIsoCharged, miniIsoNeutral, relIso;
  float JetPtRel, JetPtRatio, JetDeepBTagCSV, tkSigmaPtOverPt;
  //// Other variables to save 
  tIn->SetBranchAddress("Medium"                           , &medium                           );
  tIn->SetBranchAddress("pt"                               , &pt                               );
  tIn->SetBranchAddress("eta"                              , &eta                              );
  tIn->SetBranchAddress("SIP"                              , &sip3d                            );
  tIn->SetBranchAddress("logdxy"                           , &logdxy                           );
  tIn->SetBranchAddress("logdz"                            , &logdz                            );
  tIn->SetBranchAddress("segmentCompatibility"             , &segmComp                         );
  tIn->SetBranchAddress("JetNDauCharged"                   , &trackMult                        );
  tIn->SetBranchAddress("miniIsoCharged"                   , &miniIsoCharged                   );
  tIn->SetBranchAddress("miniIsoNeutral"                   , &miniIsoNeutral                   );
  tIn->SetBranchAddress("JetPtRel"                         , &JetPtRel                         );
  tIn->SetBranchAddress("JetPtRatio"                       , &JetPtRatio                       );
  tIn->SetBranchAddress("JetDeepBTagCSV"                   , &JetDeepBTagCSV                   );
  tIn->SetBranchAddress("relIso"                           , &relIso                           );
  tIn->SetBranchAddress("tkSigmaPtOverPt"                  , &tkSigmaPtOverPt                  );

  TFile *fOut = new TFile(fileName + "_withIsoAndMva.root", "RECREATE");
  fOut->mkdir("tpTree")->cd();
  TTree *tOut = tIn->CloneTree(0);

  int TTVLoose, TTVLeptonMvaTTZ3l, TTVLeptonMvaTTZ4l, TTVLeptonMvaTTW, TTVLeptonMvatZq, tkSigmaPtOverPtCut;
  tOut->Branch("TTVLoose",           &TTVLoose,           "TTVLoose/I");
  tOut->Branch("TTVLeptonMvaTTZ3l",  &TTVLeptonMvaTTZ3l,  "TTVLeptonMvaTTZ3l/I");
  tOut->Branch("TTVLeptonMvaTTZ4l",  &TTVLeptonMvaTTZ4l,  "TTVLeptonMvaTTZ4l/I");
  tOut->Branch("TTVLeptonMvaTTW",    &TTVLeptonMvaTTW,    "TTVLeptonMvaTTW/I");
  tOut->Branch("TTVLeptonMvatZq",    &TTVLeptonMvatZq,    "TTVLeptonMvatZq/I");
  tOut->Branch("tkSigmaPtOverPtCut", &tkSigmaPtOverPtCut, "tkSigmaPtOverPtCut/I");

  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  reader->AddVariable( "pt",                   &pt);
  reader->AddVariable( "eta",                  &eta);
  reader->AddVariable( "trackMultClosestJet",  &trackMult);
  reader->AddVariable( "miniIsoCharged",       &miniIsoCharged);
  reader->AddVariable( "miniIsoNeutral",       &miniIsoNeutral);
  reader->AddVariable( "pTRel",                &JetPtRel);
  reader->AddVariable( "ptRatio",              &JetPtRatio);
  reader->AddVariable( "relIso",               &relIso);
  reader->AddVariable( "deepCsvClosestJet",    &JetDeepBTagCSV);
  reader->AddVariable( "sip3d",                &sip3d);
  reader->AddVariable( "dxy",                  &logdxy);
  reader->AddVariable( "dz",                   &logdz);
  reader->AddVariable( "segmentCompatibility", &segmComp);

  reader->BookMVA("BDTG method", "../data/mu_BDTG_TTV_" + year + ".weights.xml");

  int step = tIn->GetEntries()/1000;
  double evDenom = 100.0/double(tIn->GetEntries());
  TStopwatch timer;
  timer.Start();
  for (int i=0, n=tIn->GetEntries(); i<n; ++i) {
    tIn->GetEntry(i);

    JetPtRatio         = std::min((double)JetPtRatio, 1.5);
    JetDeepBTagCSV     = std::max((std::isnan(JetDeepBTagCSV) ? 0. : JetDeepBTagCSV), 0.);

    float mva          = reader->EvaluateMVA("BDTG method");
    TTVLoose           = std::abs(std::exp(logdxy)) < 0.05 and std::abs(std::exp(logdz)) < 0.1 and std::abs(sip3d) < 8 and (miniIsoCharged+miniIsoNeutral) < 0.4;
    TTVLeptonMvaTTZ4l  = mva > -0.4 and JetDeepBTagCSV < (year=="2016" ? 0.8958 : 0.8001);
    TTVLeptonMvaTTZ3l  = mva >  0.4 and JetDeepBTagCSV < (year=="2016" ? 0.8958 : 0.8001);
    TTVLeptonMvaTTW    = mva >  0.6 and JetDeepBTagCSV < (year=="2016" ? 0.8958 : 0.8001);
    TTVLeptonMvatZq    = mva >  0.8 and JetDeepBTagCSV < (year=="2016" ? 0.8958 : 0.8001);
    tkSigmaPtOverPtCut = tkSigmaPtOverPt < 0.2;
    // Fill the output tree
    tOut->Fill();
    //if (i > 10000) break;
    if ((i+1) % step == 0) { 
      double totalTime = timer.RealTime()/60.; timer.Continue();
      double fraction = double(i+1)/double(n+1), remaining = totalTime*(1-fraction)/fraction;
      printf("Done %9d/%9d   %5.1f%%   (elapsed %5.1f min, remaining %5.1f min)\n", i, n, i*evDenom, totalTime, remaining); 
      fflush(stdout); 
    }
  }

  tOut->AutoSave(); // according to root tutorial this is the right thing to do
  fOut->Close();
}
