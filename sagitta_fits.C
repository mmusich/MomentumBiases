// Standalone code to fit for sagitta

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

using namespace ROOT;
using namespace ROOT::VecOps;

// MuonGenMatchedIndex
// works as it is because there is only one muon pair selected  per event
float posMuonGenMatchedIndex(RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecD GenPart_pt){
  int index = 999;
  float posGenPt = 0.0;
  for(int i=0;i<GenPart_status.size();i++){
    if (GenPart_status[i] == 1 && GenPart_pdgId[i] == 13 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 23){
      index = i;
    }
  }
  posGenPt = GenPart_pt[index];
  return posGenPt;
}

float negMuonGenMatchedIndex(RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecD GenPart_pt){
  int index = 999;
  float negGenPt = 0.0;
  for(int i=0;i<GenPart_status.size();i++){
    if (GenPart_status[i] == 1 && GenPart_pdgId[i] == -13 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 23){
      index = i;
    }
  }
  negGenPt = GenPart_pt[index];
  return negGenPt;
}

//dxy_significance

RVecD dxy_significance(RVecD Muon_dxy, RVecD Muon_dxyErr){
  return abs(Muon_dxy)/Muon_dxyErr;
}

// MuonisGood

RVecB MuonisGood(RVecD Muon_pt, RVecD Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecD Muon_pfRelIso04_all, RVecD dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
       if (Muon_pt[i] > 10 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
            muonisgood.push_back(1);
       }
  }
  return muonisgood;
}

//Pairs

RVec<std::tuple<int,int,double>> pairs(RVecD Muon_pt, RVecD Muon_charge, RVecD Muon_eta, RVecD Muon_phi, RVecB MuonisGood, RVecD Muon_dxy, RVecD Muon_dz, double rest_mass){

  RVec<std::tuple<int,int,double>> pairs; // <pos_muon_index, neg_muon_index, pair_inv_mass>
  std::tuple<int,int,double> temp;

  for(int i=1;i<Muon_pt.size();i++){
  for(int j=0;j<i;j++){
            if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
                 // Define opening angle
                 double gamma_angle = 0;
                 gamma_angle = acos(sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*cos(Muon_phi[i])*cos(Muon_phi[j]) + sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*sin(Muon_phi[i])*sin(Muon_phi[j]) + cos(2*atan(exp((-1)*Muon_eta[i])))*cos(2*atan(exp((-1)*Muon_eta[j]))));
                 if(gamma_angle > 3.141592653/4){
                      TLorentzVector firstTrack, secondTrack, mother;
                      firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);  
                      secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
                      mother = firstTrack + secondTrack;
                      double mll = mother.M();
                      if(75<mll && mll<105){
                           if(Muon_charge[i]==1){
                                temp=make_tuple(i,j,mll);
                           } else {
                                temp=make_tuple(j,i,mll);
                           }
                           pairs.push_back(temp);
                      }
                 }
            }
       }
  }
  if(pairs.size()>1){
       double diff=100.0;
       int best=0;
       for(int i=0;i<pairs.size();i++){
            if(abs(get<2>(pairs.at(i))-91)<diff){
                 diff=(abs(get<2>(pairs.at(i))-91));
                 best=i;
            }
       }
       RVec<std::tuple<int,int,double>> pairss;
       pairss.push_back(pairs.at(best));
       return pairss;
  }
  return pairs;
}

int frame(){

  auto outFileName = "snapshot_output.root";
  
  TChain chain("Events");
  //chain.Add("/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_5.root");
  //chain.Add("/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_6.root");
  chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_773.root");
  chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_772.root");
  
  RDataFrame df(chain);

  auto d1 = df.Filter("HLT_IsoMu24 == 1")
                .Filter("nMuon >= 2")
                .Filter("PV_npvsGood >= 1");

  auto d2 = d1.Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)")
                .Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)")
                .Define("pairs", "pairs(Muon_pt, Muon_charge, Muon_eta, Muon_phi, MuonisGood, Muon_dxy, Muon_dz, 0.105658)"); // muMass = 0.105658 GeV
  
  auto d3 = d2.Filter("pairs.size()==1");

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  auto d4 = d3.Define("posTrackPt","float posTrackPt; posTrackPt=Muon_pt[get<0>(pairs.at(0))]; return posTrackPt;")
              .Define("negTrackPt","float negTrackPt; negTrackPt=Muon_pt[get<1>(pairs.at(0))]; return negTrackPt;")
              .Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs.at(0))]; return posTrackEta;")
	      .Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs.at(0))]; return negTrackEta;")
              .Define("mll","std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){temporary.push_back(get<2>(pairs.at(i)));} return temporary;");

  auto d5 = d4.Define("posGenPt", "posMuonGenMatchedIndex(GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)")
              .Define("negGenPt", "negMuonGenMatchedIndex(GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)");

  auto d6 = d5.Define("posPtDiff","float posPtDiff  = posTrackPt - posGenPt; return posPtDiff;")
              .Define("negPtDiff","float negPtDiff  = negTrackPt - negGenPt; return negPtDiff;");
  
  d6.Snapshot("Events", outFileName, {"Muon_pt", "Muon_charge", "posTrackPt", "posGenPt", "negTrackPt", "negGenPt", "posPtDiff"});
  
  return 0; 

}
