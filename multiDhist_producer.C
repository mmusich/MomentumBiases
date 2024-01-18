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
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(4357 + i*10) );
  }

  //Pairs
  int option = 1; //ATTENTION, 1 gives pT smear

  auto pairs = [&](unsigned int nslot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi)->std::tuple<int,int,float,float,float,float,float,float,float,float,float,float>{

    RVec<std::tuple<int,int,float,float,float,float,float,float,float,float,float,float>> pairs; // <pos_muon_index, neg_muon_index, mll_reco (or gen smeared if option==1), mll_gen, mll_diff, posPt_reco(or smeared), negPt_reco(or smeared), posPt_gen, negPt_gen, mll_jac_alpha_weight, mll_jac_beta_weight, smear_beta_weight>
    std::tuple<int,int,float,float,float,float,float,float,float,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float smear_pt, mean, width, beta=0.95, firstPt_reco, secondPt_reco, mll_reco=0.0, firstPt_gen, secondPt_gen, smear_beta_weight=0.0, smear_beta_weight_first_term, smear_beta_weight_second_term;
    
    for(int i=1;i<Muon_pt.size();i++){
      if(MuonisGood[i]){
	for(int j=0;j<i;j++){
	  if(MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	    TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen;
	    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	    mother = firstTrack + secondTrack;
	    mll_reco = mother.M();
	    firstPt_reco = Muon_pt[i];
	    secondPt_reco = Muon_pt[j];
	    
	    //float mll_reco = r->Gaus(mother.M(), 0.03*mother.M());
	    if(75.0<mll_reco && mll_reco<105.0){ //Cut in mll
	      //Gen match
	      bool firstGenMatched = false, secondGenMatched = false;
	      for (int k=0;k<GenPart_eta.size();k++){
		if (GenPart_status[k]==1 && abs(GenPart_pdgId[k])==13 && GenPart_pdgId[GenPart_genPartIdxMother[k]]==23){ // mu(-) has PDGID 13
		  if(pow(pow(GenPart_eta[k]-Muon_eta[i],2) + pow(GenPart_phi[k]-Muon_phi[i],2),0.5)<0.3){
		    firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    firstPt_gen = GenPart_pt[k];
		    if(option==1){
		      mean = GenPart_pt[k]; //beta = 1 
		      width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		      smear_pt = rans[nslot]->Gaus(mean, width);
		      //Define pt+ (not +, but it doesn't matter, I mean the 1st side) side of the smear_beta_weight 
		      smear_beta_weight_first_term = (rans[nslot]->Gaus(mean*beta, width)) / (rans[nslot]->Gaus(mean, width)) ;
		      //Overwriting reco track
		      firstTrack.SetPtEtaPhiM(smear_pt, GenPart_eta[k], GenPart_phi[k], rest_mass);
		      firstPt_reco = smear_pt; //overwriting reco pt as well, for bining purposes
		    }
		    firstGenMatched = true;
		    if(secondGenMatched == true){break;}
		  } else if(pow(pow(GenPart_eta[k]-Muon_eta[j],2) + pow(GenPart_phi[k]-Muon_phi[j],2),0.5)<0.3){
		    secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    secondPt_gen = GenPart_pt[k];
		    if(option==1){
		      //width = (0.0084*abs(GenPart_eta[k])+0.01)*GenPart_pt[k];
                      mean = GenPart_pt[k]; //beta = 1
		      width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		      smear_pt = rans[nslot]->Gaus(mean, width); //???????????????????? is it a problem that this is the same nslot?
		      //Define pt- (not -, but it doesn't matter, I mean the 2nd side) side of the smear_beta_weight
		      smear_beta_weight_second_term = (rans[nslot]->Gaus(mean*beta, width)) / (rans[nslot]->Gaus(mean, width)) ;
		      //Overwriting reco track
		      secondTrack.SetPtEtaPhiM(smear_pt, GenPart_eta[k], GenPart_phi[k], rest_mass);
		      secondPt_reco = smear_pt; //Overwriting reco pt as well, for bining purposes
		    }
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
		smear_beta_weight = smear_beta_weight_first_term * smear_beta_weight_second_term;
	      }

	      motherGen = firstGenTrack + secondGenTrack;
	      float mll_gen = motherGen.M();
	      //attention
	      //mll_reco = rans[nslot]->Gaus(motherGen.M(), 0.03*motherGen.M());
	      float mll_diff = mll_reco - mll_gen;
	      float mll_jac_alpha_weight = (mll_reco - mll_gen)*(mll_reco - mll_gen);
	      float mll_jac_beta_weight = (mll_reco - mll_gen)*mll_gen;

	      if(Muon_charge[i]==1){
		temp=make_tuple(i,j,mll_reco,mll_gen,mll_diff,firstPt_reco,secondPt_reco,firstPt_gen,secondPt_gen,mll_jac_alpha_weight,mll_jac_beta_weight,smear_beta_weight);
	      } else {
		temp=make_tuple(j,i,mll_reco,mll_gen,mll_diff,secondPt_reco,firstPt_reco,secondPt_gen,firstPt_gen,mll_jac_alpha_weight,mll_jac_beta_weight,smear_beta_weight);
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
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    return pair_to_return;
  };

  auto d2 = d1.DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"});
    
  auto dz = d2.Define("mll_reco","return get<2>(pairs);"); 
  auto d3 = dz.Filter("mll_reco>70"); // this means only events with one mu pair are kept 

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  auto d4 = d3.Define("mll_gen","return get<3>(pairs);")
    .Define("mll_diff","return get<4>(pairs);")
    .Define("jac_alpha_weight","return get<9>(pairs)*std::copysign(1.0, genWeight);")
    .Define("jac_beta_weight","return get<10>(pairs)*std::copysign(1.0, genWeight);")
    .Define("smear_beta_weight","return get<11>(pairs)*std::copysign(1.0, genWeight);")
    .Define("weight", "std::copysign(1.0, genWeight)")
    .Define("posTrackPt","float posTrackPt = get<5>(pairs); return posTrackPt;")
    .Define("negTrackPt","float negTrackPt = get<6>(pairs); return negTrackPt;")
    .Define("posPtGen","float posPtGen = get<7>(pairs); return posPtGen;")
    .Define("negPtGen","float negPtGen = get<8>(pairs); return negPtGen;")
    .Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs)]; return posTrackEta;")
    .Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs)]; return negTrackEta;");

  //Save tree for debugging
  //TFile *f1 = new TFile("snapshot_output.root","RECREATE");
  //d4.Snapshot("Events", "snapshot_output.root", {"GenPart_status", "GenPart_pdgId", "GenPart_pt", "Muon_pt" ,"mll_reco", "mll_gen", "mll_diff","jac_weight"});

  //Tree to pass to fitting script
  //TFile *f2 = new TFile("tree_output.root","RECREATE");
  //d4.Snapshot("Events", "tree_output.root",{"posTrackPt", "posTrackEta", "negTrackPt", "negTrackEta", "mll_reco", "mll_gen", "mll_diff", "weight"});
  
  
  //Control histograms
  auto mll_reco = d4.Histo1D({"mll_smear", "mll inclusive all bins", 20, 75.0, 105.0},"mll_reco","weight"); //ATTENTION change title reco/smear according to option
  auto mll_reco_beta_val = d4.Histo1D({"mll_smear_beta_val", "mll inclusive all bins", 20, 75.0, 105.0},"mll_reco","smear_beta_weight");
  //auto mll_diff_hist = d4.Histo1D({"mll_diff", "mll diff inclusive all bins", 24, -5.0, 5.0},"mll_diff");
  auto pt_smear = d4.Histo1D({"pt_smear", "pt smear beta = 1", 15, 25.0, 55.0},"posTrackPt","weight");
  auto pt_smear_beta_val = d4.Histo1D({"pt_smear_beta_val", "pt smear beta = 0.95", 15, 25.0, 55.0},"posTrackPt","smear_beta_weight");
  TFile f3("control_histo.root","recreate");
  mll_reco->Write();
  mll_reco_beta_val->Write();
  //mll_diff_hist->Write();
  pt_smear->Write();
  pt_smear_beta_val->Write();
  f3.Close();
  
  
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=10, nbinsmll=10, nbinseta=24, nbinspt=5;
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
    
  if (option==1){
    std::unique_ptr<TFile> f6( TFile::Open("multiD_histo_smear.root", "RECREATE") );
    //ATTENTION do not change variable name mll_reco, it's used for both reco and smear depending on option
    auto mDh_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_smear", "multi_data_histo_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","weight"});
    f6->WriteObject(mDh_smear.GetPtr(), "multi_data_histo_smear");

    //one more 5D weighted histogram where the weight is the jac weight
    //ATTENTION do not change variable name jac_alpha_weight, it's used for both reco and smear depending on option
    auto mDh_diff_squared_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_squared_smear", "multi_data_histo_diff_squared_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","jac_alpha_weight"});
    f6->WriteObject(mDh_diff_squared_smear.GetPtr(), "multi_data_histo_diff_squared_smear");

    //ATTENTION do not change variable name jac_alpha_weight, it's used for both reco and smear depending on option
    auto mDh_diff_squared_smear_control = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_squared_smear_control", "multi_data_histo_diff_squared_smear_control", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","jac_alpha_weight"});
    f6->WriteObject(mDh_diff_squared_smear_control.GetPtr(), "multi_data_histo_diff_squared_smear_control");

    //ATTENTION do not change variable name jac_beta_weight, it's used for both reco and smear depending on option
    auto mDh_jac_beta_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_beta_smear", "multi_data_histo_jac_beta_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","jac_beta_weight"});
    f6->WriteObject(mDh_jac_beta_smear.GetPtr(), "multi_data_histo_jac_beta_smear");

    //ATTENTION do not change variable name jac_beta_weight, it's used for both reco and smear depending on option
    auto mDh_smear_beta_val = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_smear_beta_val", "multi_data_histo_smear_beta_val", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","smear_beta_weight"});
    f6->WriteObject(mDh_smear_beta_val.GetPtr(), "multi_data_histo_smear_beta_val");

    //ATTENTION do not change variable name mll_diff, it's used for both reco and smear depending on option
    auto mDh_diff_smear = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear", "multi_data_histo_diff_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","weight"});
    f6->WriteObject(mDh_diff_smear.GetPtr(), "multi_data_histo_diff_smear"); 
  } else { // option==0 case
    std::unique_ptr<TFile> f5( TFile::Open("multiD_histo_reco.root", "RECREATE") );

    auto mDh_reco = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_reco", "multi_data_histo_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","weight"});
    f5->WriteObject(mDh_reco.GetPtr(), "multi_data_histo_reco");
    
    auto mDh_gen = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_gen", "multi_data_histo_gen", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges,ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_gen","weight"});
    f5->WriteObject(mDh_gen.GetPtr(), "multi_data_histo_gen");

    auto mDh_diff_reco = d4.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_reco", "multi_data_histo_diff_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","weight"});
    f5->WriteObject(mDh_diff_reco.GetPtr(), "multi_data_histo_diff_reco");
  }
  
  return 0;

}
  
