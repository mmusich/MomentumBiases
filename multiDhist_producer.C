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
RVecB MuonisGood(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVec<UChar_t> Muon_genPartFlav, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
    if (Muon_pt[i] > 10 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && Muon_genPartFlav[i]==1 && dxy_significance[i] < 4){
      muonisgood.push_back(1);
    }
  }
  return muonisgood;
}

//Rand_pt
RVecF Rand_pt(RVecF GenPart_pt, RVecF GenPart_eta){
  RVecF rand_pt;
  float x, width;
  TRandom *r = new TRandom();
  for(int i=0;i<GenPart_pt.size();i++){
    width = (0.0083*abs(GenPart_eta[i])+0.01)*GenPart_pt[i];
    x = r->Gaus(GenPart_pt[i], width);
    rand_pt.push_back(x);
  }
  return rand_pt;
}


//Pairs
std::tuple<int,int,float,float,float> pairs(RVecF Muon_pt, RVecF Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, float rest_mass, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF Rand_pt, int option){

  RVec<std::tuple<int,int,float,float,float>> pairs; // <pos_muon_index, neg_muon_index, mll_reco, mll_gen, mll_diff>
  std::tuple<int,int,float,float,float> temp, pair_to_return;

  for(int i=1;i<Muon_pt.size();i++){ //TODO change the loop to check MuonisGood[i] only once
    for(int j=0;j<i;j++){
      if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen;
	firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);  
	secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	mother = firstTrack + secondTrack;
	float mll_reco = mother.M();
	if(75.0<mll_reco && mll_reco<105.0){ //Cut in mll
	  //Gen match
	  bool firstGenMatched = false, secondGenMatched = false;
	  for (int k=0;k<GenPart_eta.size();k++){
	    if (GenPart_status[k]==1 && abs(GenPart_pdgId[k])==13 && GenPart_pdgId[GenPart_genPartIdxMother[k]]==23){ // mu(-) has PDGID 13
	      if(pow(pow(GenPart_eta[k]-Muon_eta[i],2) + pow(GenPart_phi[k]-Muon_phi[i],2),0.5)<0.3){ 
		firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		if(option==1){ firstTrack.SetPtEtaPhiM( Rand_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass); }
		firstGenMatched = true;
		if(secondGenMatched == true){break;}
	      } else if(pow(pow(GenPart_eta[k]-Muon_eta[j],2) + pow(GenPart_phi[k]-Muon_phi[j],2),0.5)<0.3){
		secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		if(option==1){ secondTrack.SetPtEtaPhiM( Rand_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass); }
		secondGenMatched = true;
		if(firstGenMatched == true){break;}
	      }
	    }
	  }	  
	  if(firstGenMatched == false || secondGenMatched == false){
	    continue; 
	  }
	  if(option==1){
	    mother = firstTrack + secondTrack;
	    mll_reco = mother.M();
	  }
	  motherGen = firstGenTrack + secondGenTrack;
	  float mll_gen = motherGen.M();
	  float mll_diff = mll_reco - mll_gen;
	  
	  if(Muon_charge[i]==1){
	    temp=make_tuple(i,j,mll_reco,mll_gen,mll_diff);
	  } else {
	    temp=make_tuple(j,i,mll_reco,mll_gen,mll_diff);
	  }
	  pairs.push_back(temp);
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
    pair_to_return=make_tuple(0,0,0.0,0.0,0.0);
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
    .Define("Rand_pt","Rand_pt(GenPart_pt, GenPart_eta)")
    .Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, Muon_genPartFlav, dxy_significance)")
    .Define("pairs", "pairs(Muon_pt, Muon_charge, Muon_eta, Muon_phi, MuonisGood, Muon_dxy, Muon_dz, 0.105658, GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, Rand_pt, 1)") //ATTENTION to last argument, the option for rand
    // muMass = 0.105658 GeV
    .Define("mll_reco","return get<2>(pairs);"); 
 
  auto d3 = d2.Filter("mll_reco>70"); // this means only events with one mu pair are kept 

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  auto d4 = d3.Define("mll_gen","return get<3>(pairs);")
    .Define("mll_diff","return get<4>(pairs);")
    .Define("weight", "std::copysign(1.0, genWeight)")
    .Define("posTrackPt","float posTrackPt; posTrackPt=Muon_pt[get<0>(pairs)]; return posTrackPt;")
    .Define("negTrackPt","float negTrackPt; negTrackPt=Muon_pt[get<1>(pairs)]; return negTrackPt;")
    .Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs)]; return posTrackEta;")
    .Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs)]; return negTrackEta;");

  //Save tree for debugging
  TFile *f1 = new TFile("snapshot_output.root","RECREATE");
  d4.Snapshot("Events", "snapshot_output.root", {"GenPart_status", "GenPart_pdgId", "Rand_pt", "GenPart_pt", "Muon_pt" ,"mll_reco", "mll_gen", "mll_diff"});

  //Tree to pass to fitting script
  //TFile *f2 = new TFile("tree_output.root","RECREATE");
  //d4.Snapshot("Events", "tree_output.root",{"posTrackPt", "posTrackEta", "negTrackPt", "negTrackEta", "mll_reco", "mll_gen", "mll_diff", "weight"});
  
  /*
  //Control histograms
  auto mll_reco_hist = d4.Histo1D({"mll_reco", "mll inclusive all bins", 45, 75.0, 105.0},"mll_reco");
  auto mll_diff_hist = d4.Histo1D({"mll_diff", "mll diff inclusive all bins", 24, -6.0, 6.0},"mll_diff");
  TFile f3("control_histo.root","recreate");
  mll_reco_hist->Write();
  mll_diff_hist->Write();
  f3.Close();
  */
  
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=18, nbinsmll=18, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, ptbinranges, mll_diffbinranges, mllbinranges;

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  double myetaboundaries[etabinranges.size()];
  for (int i=0; i<etabinranges.size(); i++){
    myetaboundaries[i] = etabinranges[i];
  }

  for (int i=0; i<=nbinsmll_diff; i++) mll_diffbinranges.push_back(-6.0 + i * 12.0/nbinsmll_diff);
  for (int i=0; i<=nbinsmll; i++) mllbinranges.push_back(75.0 + i * (105.0-75.0)/nbinsmll);

  std::unique_ptr<TFile> f4( TFile::Open("control_bin_histo.root", "RECREATE") );  

  //pT bin optimisation starts
  auto pt_pos_uni = d4.Histo1D({"pt_pos_uni", "pt mu+", nbinspt*3, ptlow, pthigh},"posTrackPt");
  f4->WriteObject(pt_pos_uni.GetPtr(), "pt_pos_uni");
  // Get quartiles
  double xq[nbinspt+1], myptboundaries[nbinspt+1];
  for (int i=0;i<=nbinspt;i++) xq[i] = float(i)/nbinspt;
  pt_pos_uni->GetQuantiles(nbinspt+1,myptboundaries, xq);

  ptbinranges.clear();
  std::cout<<"ptbinranges = [";
  for (int i=0; i<=nbinspt; i++){
    ptbinranges.push_back(myptboundaries[i]);
    std::cout<<ptbinranges[i]<<", ";
  }
  std::cout<<"] \n";

  //TH1 in pT+ with variable bin size -> should be uniform
  auto pt_pos = d4.Histo1D({"pt_pos", "pt mu+", nbinspt, myptboundaries},"posTrackPt");
  f4->WriteObject(pt_pos.GetPtr(), "pt_pos");

  auto pt_eta_pos = d4.Histo2D({"pt_eta_pos", "pt eta mu+", nbinseta, myetaboundaries, nbinspt, myptboundaries},"posTrackEta", "posTrackPt");
  f4->WriteObject(pt_eta_pos.GetPtr(), "pt_eta_pos");

  std::unique_ptr<TFile> f5( TFile::Open("multiD_histo.root", "RECREATE") );

  //  auto mDh = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_reco", "multi_data_histo_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","weight"});
  // f5->WriteObject(mDh.GetPtr(), "multi_data_histo_reco");
  //  delete gROOT->FindObject("multi_data_histo_reco");

  auto mDh = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_rand", "multi_data_histo_rand", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","weight"});
  f5->WriteObject(mDh.GetPtr(), "multi_data_histo_rand");
  delete gROOT->FindObject("multi_data_histo_rand");

  mDh = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_gen", "multi_data_histo_gen", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges,ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_gen","weight"});
  f5->WriteObject(mDh.GetPtr(), "multi_data_histo_gen");
  delete gROOT->FindObject("multi_data_histo_gen");

  mDh = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff", "multi_data_histo_diff", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","weight"});
  f5->WriteObject(mDh.GetPtr(), "multi_data_histo_diff"); 
  
  return 0;

}
  
