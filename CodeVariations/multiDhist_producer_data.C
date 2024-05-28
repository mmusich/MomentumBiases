// Standalone code to make multiD histograms from data for the sagitta fit

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

int frame(){

  ROOT::EnableImplicitMT(128);

  TChain chain("Events");
  //chain.Add("/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_5.root");
  //chain.Add("/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_6.root");
  
  string line;

  ifstream file("data_filenames.txt");
  if (file.is_open()) {
    while (getline(file, line)) {
      chain.Add(line.c_str());
    }
    file.close();
  }  

  RDataFrame df(chain);

  auto d0 = df.Filter("HLT_IsoMu24 == 1")
    .Filter("nMuon >= 2")
    .Filter("PV_npvsGood >= 1");

  auto d1 = d0.Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)")
    .Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)");

  //Pairs
  auto pairs = [&](RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz)->std::tuple<int,int,float,float,float>{

    RVec<std::tuple<int,int,float,float,float>> pairs; // <pos_muon_index, neg_muon_index, mll_data, pos_muon_pt, neg_muon_pt>
    std::tuple<int,int,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float width, firstPt_data, secondPt_data, mll_data=0.0;
    
    for(int i=1;i<Muon_pt.size();i++){
      if(MuonisGood[i]){
	for(int j=0;j<i;j++){
	  if(MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	    TLorentzVector firstTrack, secondTrack, mother;
	    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	    mother = firstTrack + secondTrack;
	    mll_data = mother.M();
	    firstPt_data = Muon_pt[i];
	    secondPt_data = Muon_pt[j];

	    if(75.0<mll_data && mll_data<105.0){ //Cut in mll
	      if(Muon_charge[i]==1){
		temp=make_tuple(i,j,mll_data,firstPt_data,secondPt_data);
	      } else {
		temp=make_tuple(j,i,mll_data,secondPt_data,firstPt_data);
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
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0);
    }
    return pair_to_return;
  };

  auto d2 = d1.Define("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz"});
    
  auto dz = d2.Define("mll_data","return get<2>(pairs);"); 
  auto d3 = dz.Filter("mll_data>70"); // this means only events with one mu pair are kept 

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  auto d4 = d3.Define("posTrackPt","float posTrackPt = get<3>(pairs); return posTrackPt;")
    .Define("negTrackPt","float negTrackPt = get<4>(pairs); return negTrackPt;")
    .Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs)]; return posTrackEta;")
    .Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs)]; return negTrackEta;");

  //Save tree for debugging
  //TFile *f1 = new TFile("snapshot_output_data.root","RECREATE");
  //d4.Snapshot("Events", "snapshot_output.root", {"Muon_pt" ,"mll_data"});
  
  /*
  //Control histograms
  auto mll_data_hist = d4.Histo1D({"mll_data", "mll inclusive all bins", 45, 75.0, 105.0},"mll_data");
  TFile f3("control_histo_data.root","recreate");
  mll_data_hist->Write();
  f3.Close();
  */
  
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll=5, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, mllbinranges, ptbinranges{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};
  
  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  for (int i=0; i<=nbinsmll; i++) mllbinranges.push_back(75.0 + i * (105.0-75.0)/nbinsmll);

  std::unique_ptr<TFile> f5( TFile::Open("multiD_histo_data.root", "RECREATE") );
  
  auto mDh_data = d4.HistoND<float, float, float, float, float>({"multi_data_histo_data", "multi_data_histo_data", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_data"});
  f5->WriteObject(mDh_data.GetPtr(), "multi_data_histo_data");
  
  return 0;

}
  
