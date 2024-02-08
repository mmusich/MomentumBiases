// Standalone code to make multiD histograms for the sagitta fit

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

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

  auto d0 = df.Filter("HLT_IsoMu24 == 1")
    .Filter("nMuon >= 2")
    .Filter("PV_npvsGood >= 1");

  auto d1 = d0.Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)")
    .Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, Muon_genPartFlav, dxy_significance)");

  unsigned int nslots = d1.GetNSlots();
  std::vector<TRandom3*> rans = {}, rans_2 = {}, rans_3 = {}, rans_4 = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(4357 + i*10) );
    rans_2.emplace_back( new TRandom3(3951 + i*10) );
    rans_3.emplace_back( new TRandom3(5193 + i*10) );
    rans_4.emplace_back( new TRandom3(9361 + i*10) );
  }

  //Pairs
  auto pairs = [&](unsigned int nslot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi)->std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float>{

    RVec<std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float>> pairs; // <pos_muon_index, neg_muon_index, mll_reco, posPt_reco, negPt_reco, mll_diff_reco, mll_smear, posPt_smear, negPt_smear, mll_diff_smear, mll_gen, posPt_gen, negPt_gen, mll_diff_squared_smear, smear_beta_weight, posPt_smear_beta_val, negPt_smear_beta_val, mll_diff_squared_reco>
    std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco, firstPt_smear, secondPt_smear, mll_smear, firstPt_gen, secondPt_gen, mll_gen, firstPt_smear_beta_val, secondPt_smear_beta_val, smear_beta_weight;
    float smear_pt, mean, width, beta=0.001, smear_beta_weight_first_term, smear_beta_weight_second_term;
    
    for(int i=1;i<Muon_pt.size();i++){
      if(MuonisGood[i]){
	for(int j=0;j<i;j++){
	  if(MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	    TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen, firstSmearTrack, secondSmearTrack, motherSmear;
	    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	    mother = firstTrack + secondTrack;
	    mll_reco = mother.M();
	    firstPt_reco = Muon_pt[i];
	    secondPt_reco = Muon_pt[j];
	    
	    if(75.0<mll_reco && mll_reco<105.0){ //Cut in mll
	      //Gen match
	      bool firstGenMatched = false, secondGenMatched = false;
	      for (int k=0;k<GenPart_eta.size();k++){
		if (GenPart_status[k]==1 && abs(GenPart_pdgId[k])==13 && GenPart_pdgId[GenPart_genPartIdxMother[k]]==23){ // mu(-) has PDGID 13
		  if(pow(pow(GenPart_eta[k]-Muon_eta[i],2) + pow(GenPart_phi[k]-Muon_phi[i],2),0.5)<0.3){
		    firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    firstPt_gen = GenPart_pt[k];
		    
		    //smear 1st muon
		    mean = GenPart_pt[k]; //beta = 1 
		    width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		    firstPt_smear = rans[nslot]->Gaus(mean, width);
		    firstSmearTrack.SetPtEtaPhiM(firstPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_val, weight for beta != 1
		    smear_beta_weight_first_term = TMath::Gaus(firstPt_smear, mean*beta, width) / TMath::Gaus(firstPt_smear, mean, width);
		    firstPt_smear_beta_val = rans_3[nslot]->Gaus(mean*beta, width);
		    		    
		    firstGenMatched = true;
		    if(secondGenMatched == true){break;}
		  } else if(pow(pow(GenPart_eta[k]-Muon_eta[j],2) + pow(GenPart_phi[k]-Muon_phi[j],2),0.5)<0.3){
		    secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    secondPt_gen = GenPart_pt[k];
		    
		    //smear 2nd muon
		    mean = GenPart_pt[k]; //beta = 1
		    width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		    secondPt_smear = rans_2[nslot]->Gaus(mean, width);
		    secondSmearTrack.SetPtEtaPhiM(secondPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_val, weight for beta != 1
		    smear_beta_weight_second_term = TMath::Gaus(secondPt_smear, mean*beta, width) / TMath::Gaus(secondPt_smear, mean, width); 
		    secondPt_smear_beta_val = rans_4[nslot]->Gaus(mean*beta, width);
		    
		    secondGenMatched = true;
		    if(firstGenMatched == true){break;}
		  }
		}
	      }
	      if(firstGenMatched == false || secondGenMatched == false){
		continue;
	      }
	      
	      motherSmear = firstSmearTrack + secondSmearTrack;
	      mll_smear = motherSmear.M();
	      smear_beta_weight = smear_beta_weight_first_term * smear_beta_weight_second_term;
	      
	      motherGen = firstGenTrack + secondGenTrack;
	      float mll_gen = motherGen.M();
	      float mll_diff_reco = mll_reco - mll_gen;
	      float mll_diff_smear = mll_smear - mll_gen;
	      // save for jacobians
	      float mll_diff_squared_smear = (mll_smear - mll_gen)*(mll_smear - mll_gen);
	      float mll_diff_squared_reco = (mll_reco - mll_gen)*(mll_reco - mll_gen);
              
	      if(Muon_charge[i]==1){
		temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco,mll_diff_reco,mll_smear,firstPt_smear,secondPt_smear,mll_diff_smear,mll_gen,firstPt_gen,secondPt_gen,mll_diff_squared_smear,smear_beta_weight,firstPt_smear_beta_val,secondPt_smear_beta_val, mll_diff_squared_reco);
	      } else {
		temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco,mll_diff_reco,mll_smear,secondPt_smear,firstPt_smear,mll_diff_smear,mll_gen,secondPt_gen,firstPt_gen,mll_diff_squared_smear,smear_beta_weight,secondPt_smear_beta_val,firstPt_smear_beta_val, mll_diff_squared_reco);
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
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    return pair_to_return;
  };

  auto d2 = d1.DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"});
    
  auto dz = d2.Define("mll_reco","return get<2>(pairs);"); 
  auto d3 = dz.Filter("mll_reco>1.0"); // this means only events with one mu pair are kept 

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco, 5 mll_diff_reco, 6 mll_smear, 7 posPt_smear, 8 negPt_smear, 9 mll_diff_smear, 10 mll_gen, 11 posPt_gen, 12 negPt_gen, 13 mll_diff_squared_smear, 14 smear_beta_weight, 15 posPt_smear_beta_val, 16 negPt_smear_beta_val, 17 mll_diff_squared_reco>

  auto d4 = d3.Define("posTrackPt","return get<3>(pairs);")
    .Define("negTrackPt","return get<4>(pairs);")
    .Define("mll_diff_reco","return get<5>(pairs);")
    .Define("jacobian_weight_mll_diff_reco","return get<5>(pairs)*std::copysign(1.0, genWeight);")
    .Define("mll_smear","return get<6>(pairs);")
    .Define("posPtSmear","return get<7>(pairs);")
    .Define("negPtSmear","return get<8>(pairs);")
    .Define("mll_diff_smear","return get<9>(pairs);")
    .Define("jacobian_weight_mll_diff_smear", "return get<9>(pairs)*std::copysign(1.0, genWeight);")
    .Define("mll_gen","return get<10>(pairs);")
    .Define("posPtGen","return get<11>(pairs);")
    .Define("negPtGen","return get<12>(pairs);")
    .Define("jacobian_weight_mll_diff_squared_smear","return get<13>(pairs)*std::copysign(1.0, genWeight);")
    .Define("smear_beta_weight","return get<14>(pairs)*std::copysign(1.0, genWeight);")
    .Define("posPtSmearBetaVal","return get<15>(pairs);")
    .Define("negPtSmearBetaVal","return get<16>(pairs);")
    .Define("jacobian_weight_mll_diff_squared_reco","return get<17>(pairs)*std::copysign(1.0, genWeight);")
    .Define("weight", "std::copysign(1.0, genWeight)")
    .Define("posTrackEta","return Muon_eta[get<0>(pairs)];")
    .Define("negTrackEta","return Muon_eta[get<1>(pairs)];");

  //Save tree for debugging
  //TFile *f1 = new TFile("snapshot_output.root","RECREATE");
  //d4.Snapshot("Events", "snapshot_output.root", {"GenPart_status"});

  //Control histograms
  auto mll_smear = d4.Histo1D({"mll_smear", "mll inclusive all bins", 20, 75.0, 105.0},"mll_smear","weight"); 
  auto mll_smear_beta_val = d4.Histo1D({"mll_smear_beta_val", "mll inclusive all bins", 20, 75.0, 105.0},"mll_smear","smear_beta_weight");
  auto pt_smear = d4.Histo1D({"pt_smear", "pt smear beta = 1", 15, 25.0, 55.0},"posPtSmear","weight");
  auto pt_smear_beta_val = d4.Histo1D({"pt_smear_beta_val", "pt smear beta = 0.001", 15, 25.0, 55.0},"posPtSmear","smear_beta_weight");
  TFile f3("control_histo.root","recreate");
  mll_smear->Write();
  mll_smear_beta_val->Write();
  pt_smear->Write();
  pt_smear_beta_val->Write();
  f3.Close();
  
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=8, nbinsmll=8, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, ptbinranges, mllbinranges;
  vector<double> mll_diffbinranges;

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  double myetaboundaries[etabinranges.size()];
  for (int i=0; i<etabinranges.size(); i++){
    myetaboundaries[i] = etabinranges[i];
  }

  for (int i=0; i<=nbinsmll_diff; i++) mll_diffbinranges.push_back(-5.0 + i * 10.0/nbinsmll_diff);
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
  
  /*
  std::unique_ptr<TFile> f7( TFile::Open("multiD_histo_smear_beta_val_easy.root", "RECREATE") );
  // mll_diff_smear_beta_val_easy 
  auto mDh_diff_smear_beta_val_easy = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_beta_val_easy", "multi_data_histo_diff_smear_beta_val_easy", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","weight"});
  f7->WriteObject(mDh_diff_smear_beta_val_easy.GetPtr(), "multi_data_histo_diff_smear_beta_val_easy");
  
  */
  int option = 1; //ATTENTION, 1 saves SMEAR plots (eventually can save both reco and smear plots and give up this option thing)  
  if (option==1){
    std::unique_ptr<TFile> f6( TFile::Open("multiD_histo_smear.root", "RECREATE") ); 
    
    // mll_smear
    auto mDh_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_smear", "multi_data_histo_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","weight"});
    f6->WriteObject(mDh_smear.GetPtr(), "multi_data_histo_smear");

    // mll_diff_smear
    auto mDh_diff_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear", "multi_data_histo_diff_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","weight"});
    f6->WriteObject(mDh_diff_smear.GetPtr(), "multi_data_histo_diff_smear");

    // mll_smear weighted by jacobian_weight_mll_diff_squared_smear
    auto mDh_jac_diff_squared_smear_mll = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_squared_smear_mll", "multi_data_histo_jac_diff_squared_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_diff_squared_smear"});
    f6->WriteObject(mDh_jac_diff_squared_smear_mll.GetPtr(), "multi_data_histo_jac_diff_squared_smear_mll");

    // mll_diff_smear weighted by jacobian_weight_mll_diff_squared_smear
    auto mDh_jac_diff_squared_smear_mll_diff = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_squared_smear_mll_diff", "multi_data_histo_jac_diff_squared_smear_mll_diff", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","jacobian_weight_mll_diff_squared_smear"});
    f6->WriteObject(mDh_jac_diff_squared_smear_mll_diff.GetPtr(), "multi_data_histo_jac_diff_squared_smear_mll_diff");
    
    // mll_diff_smear weighted by jacobian_weight_mll_diff_smear
    auto mDh_jac_diff_smear_mll_diff = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_smear_mll_diff", "multi_data_histo_jac_diff_smear_mll_diff", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","jacobian_weight_mll_diff_smear"});
    f6->WriteObject(mDh_jac_diff_smear_mll_diff.GetPtr(), "multi_data_histo_jac_diff_smear_mll_diff");
    //TODO continue changing names
    // mll_smear weighted by jacobian_weight_mll_diff_smear
    auto mDh_jac_beta_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_beta_smear", "multi_data_histo_jac_beta_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jac_beta_weight_smear"});
    f6->WriteObject(mDh_jac_beta_smear.GetPtr(), "multi_data_histo_jac_beta_smear");

    // mll_diff_smear weighted by jac_beta_weight_smear
    auto mDh_jac_beta_smear_control = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_beta_smear_control", "multi_data_histo_jac_beta_smear_control", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","jac_beta_weight_smear"});
    f6->WriteObject(mDh_jac_beta_smear_control.GetPtr(), "multi_data_histo_jac_beta_smear_control");

    // mll_smear weighted by smear_beta_weight
    auto mDh_smear_beta_val = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_smear_beta_val", "multi_data_histo_smear_beta_val", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","smear_beta_weight"});
    f6->WriteObject(mDh_smear_beta_val.GetPtr(), "multi_data_histo_smear_beta_val");

    // mll_diff_smear weighted by smear_beta_weight
    auto mDh_diff_smear_beta_val = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_beta_val", "multi_data_histo_diff_smear_beta_val", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","smear_beta_weight"});
    f6->WriteObject(mDh_diff_smear_beta_val.GetPtr(), "multi_data_histo_diff_smear_beta_val");
    
    // mll_diff_smear weighted by jacobian_beta_weight_smear_mll_diff
    auto mDh_jac_beta_smear_mll_diff = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_beta_smear_mll_diff", "multi_data_histo_jac_beta_smear_mll_diff", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","jacobian_beta_weight_smear_mll_diff"});
    f6->WriteObject(mDh_jac_beta_smear_mll_diff.GetPtr(), "multi_data_histo_jac_beta_smear_mll_diff");
    
  } else { // option==0 saves RECO plots
    std::unique_ptr<TFile> f5( TFile::Open("multiD_histo_reco.root", "RECREATE") );

    auto mDh_reco = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_reco", "multi_data_histo_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","weight"});
    f5->WriteObject(mDh_reco.GetPtr(), "multi_data_histo_reco");
    
    auto mDh_gen = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_gen", "multi_data_histo_gen", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges,ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_gen","weight"});
    f5->WriteObject(mDh_gen.GetPtr(), "multi_data_histo_gen");
    
    auto mDh_diff_reco = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_reco", "multi_data_histo_diff_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff_reco","weight"});
    f5->WriteObject(mDh_diff_reco.GetPtr(), "multi_data_histo_diff_reco");
  }
  
  return 0;

}
  
