// Standalone code to make multiD histograms for the sagitta fit

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

using namespace ROOT;
using namespace ROOT::VecOps;

//dxy_significance
RVecF dxy_significance(RVecF Muon_dxy, RVecF Muon_dxyErr){
  return abs(Muon_dxy)/Muon_dxyErr;
}

// MuonisGood
RVecB MuonisGood(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
    if (Muon_pt[i] > 10 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
      muonisgood.push_back(1);
    }
  }
  return muonisgood;
}

//Pairs
std::tuple<int,int,float, float> pairs(RVecF Muon_pt, RVecF Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, float rest_mass, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi){

  RVec<std::tuple<int,int,float,float>> pairs; // <pos_muon_index, neg_muon_index, pair_inv_mass>
  std::tuple<int,int,float,float> temp, pair_to_return;

  for(int i=1;i<Muon_pt.size();i++){
    for(int j=0;j<i;j++){
      if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	// Define opening angle
	float gamma_angle = 0;
	gamma_angle = acos(sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*cos(Muon_phi[i])*cos(Muon_phi[j]) + sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*sin(Muon_phi[i])*sin(Muon_phi[j]) + cos(2*atan(exp((-1)*Muon_eta[i])))*cos(2*atan(exp((-1)*Muon_eta[j]))));
	if(gamma_angle > 3.141592653/4){
	  TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen;
	  firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);  
	  secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	  mother = firstTrack + secondTrack;
	  float mll = mother.M();
	  if(75<mll && mll<105){ //Cut in mll
	    //Gen match
	    bool firstGenMatched = false, secondGenMatched = false;
	    for (int k=0;k<GenPart_eta.size();k++){
	      if (GenPart_status[k]==1 && abs(GenPart_pdgId[k])==13 && GenPart_pdgId[GenPart_genPartIdxMother[k]]==23){ // mu(-) has PDGID 13
		if(abs(pow(GenPart_eta[k],2) + pow(GenPart_phi[k],2) - pow(Muon_eta[i],2) - pow(Muon_phi[i],2))<0.15){ 
		  //CHECK IF deltaR value is large enough, maybe more than a gen muon satisfies the condition
		  firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		  firstGenMatched = true;
		} else if(abs(pow(GenPart_eta[k],2) + pow(GenPart_phi[k],2) - pow(Muon_eta[j],2) - pow(Muon_phi[j],2))<0.15){
		  secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		  secondGenMatched = true;
		}
	      }
	    }	  
	    if(firstGenMatched == false || secondGenMatched == false){
	      continue; 
	    }
	    motherGen = firstGenTrack + secondGenTrack;
	    float mllGen = motherGen.M();
	    float mll_diff = mll - mllGen;
	    
	    if(Muon_charge[i]==1){
	      temp=make_tuple(i,j,mll,mll_diff);
	    } else {
	      temp=make_tuple(j,i,mll,mll_diff);
	    }
	    pairs.push_back(temp);
	  }
	}
      }
    }
  }

  if(pairs.size()==1){
    pair_to_return=pairs.at(0);
  } else if(pairs.size()>1){
    float diff=100.0;
    int best=0;
    for(int i=0;i<pairs.size();i++){
      if(abs(get<2>(pairs.at(i))-91)<diff){
	diff=(abs(get<2>(pairs.at(i))-91));
	best=i;
      }
    }
    pair_to_return=pairs.at(best);
  } else {
    pair_to_return=make_tuple(0,0,0.0,0.0);
 }
  return pair_to_return;
}

int frame(){

  ROOT::EnableImplicitMT(128);

  TChain chain("Events");
  //chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_773.root");
  //chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_772.root");
  
  string line;

  ifstream file("MCfilenames.txt");
  if (file.is_open()) {
    while (getline(file, line)) {
      chain.Add(line.c_str());
    }
    file.close();
  }  

  RDataFrame df(chain);

  auto d1 = df.Filter("HLT_IsoMu24 == 1")
    .Filter("nMuon >= 2")
    .Filter("PV_npvsGood >= 1");

  auto d2 = d1.Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)")
    .Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)")
    .Define("pairs", "pairs(Muon_pt, Muon_charge, Muon_eta, Muon_phi, MuonisGood, Muon_dxy, Muon_dz, 0.105658, GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi)") 
    // muMass = 0.105658 GeV
    .Define("mll_diff","return get<3>(pairs);")
    .Define("mll","return get<2>(pairs);");  
  auto d3 = d2.Filter("mll>70"); // this means only events with one mu pair are kept 

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  auto d4 = d3.Define("posTrackPt","float posTrackPt; posTrackPt=Muon_pt[get<0>(pairs)]; return posTrackPt;")
    .Define("negTrackPt","float negTrackPt; negTrackPt=Muon_pt[get<1>(pairs)]; return negTrackPt;")
    .Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs)]; return posTrackEta;")
    .Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs)]; return negTrackEta;");

  TFile *f1 = new TFile("snapshot_output.root","RECREATE");
  d4.Snapshot("Events", "snapshot_output.root", {"Muon_pt", "GenPart_pdgId", "GenPart_pt", "GenPart_status", "mll", "mll_diff"});
  
  /*
  //number of eta pt bins should match or be larger than the one chosen later
  int nbinseta=24, nbinspt=50;
  float ptlow=25.0, pthigh=55.0;
  //larger bins while debugging
  //int nbinseta=6, nbinspt=2;
    
  auto multi_hist = d4.HistoND<float, float, float, float, RVecF, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, 24}, {-2.4,ptlow,-2.4,ptlow,-6}, {2.4,pthigh,2.4,pthigh,6}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff","genWeight"});

  TFile f2("multiD_histo.root","recreate");
  multi_hist->Write(); 
  f2.Close();
  */
  return 0;

}
  
