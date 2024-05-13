// Standalone code to make multiD histograms for the A,e,M fit

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

#include <sys/time.h>

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
    } else {
      muonisgood.push_back(0);
    }
  }
  return muonisgood;
}

int multiDhist_producer(){

  double t1(0.);
  // Get start time
  struct timeval tv_start, tv_stop;
  gettimeofday(&tv_start, 0);

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
  auto pairs = [&](unsigned int nslot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi)->std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float>{

    // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
    RVec<std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float>> pairs; 
    std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco, firstPt_smear, secondPt_smear, mll_smear, firstPt_gen, secondPt_gen, mll_gen, firstPt_smear_beta_val, secondPt_smear_beta_val, smear_beta_weight;
    float smear_pt, mean, width, smeared_mean, smear_beta_weight_first_term, smear_beta_weight_second_term;
    float A = 0.0004, e = 0.002, M = 0.00001; //for now constant per eta bin
        
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
		  if(pow(pow(GenPart_eta[k]-Muon_eta[i],2) + pow(TVector2::Phi_mpi_pi(GenPart_phi[k]-Muon_phi[i]),2),0.5)<0.3){
		    firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    firstPt_gen = GenPart_pt[k];
		            
		    //smear 1st muon
		    mean = GenPart_pt[k]; //beta = 1 
		    width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		    firstPt_smear = rans[nslot]->Gaus(mean, width);
		    firstSmearTrack.SetPtEtaPhiM(firstPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_val, weight for beta != 1 use Eqn 22 in note AN2021_131_v5
		    if (GenPart_pdgId[k]>0){ //mu(-)
		      smeared_mean = 2.0*e/((1.0+A)-pow((1.0+A)*(1.0+A)+4.0*e*(-M-1.0/mean),0.5));
		    } else { //mu(+)
		      smeared_mean = 2.0*e/((1.0+A)-pow((1.0+A)*(1.0+A)+4.0*e*(M-1.0/mean),0.5));
		    }
		    smear_beta_weight_first_term = TMath::Gaus(firstPt_smear, smeared_mean, width) / TMath::Gaus(firstPt_smear, mean, width);
		    firstPt_smear_beta_val = rans[nslot]->Gaus(smeared_mean, width);
		            
		    firstGenMatched = true;
		    if(secondGenMatched == true){break;}
		  } else if(pow(pow(GenPart_eta[k]-Muon_eta[j],2) + pow(TVector2::Phi_mpi_pi(GenPart_phi[k]-Muon_phi[j]),2),0.5)<0.3){
		    secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    secondPt_gen = GenPart_pt[k];
		            
		    //smear 2nd muon
		    mean = GenPart_pt[k]; //beta = 1
		    width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		    secondPt_smear = rans[nslot]->Gaus(mean, width);
		    secondSmearTrack.SetPtEtaPhiM(secondPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_val, weight for beta != 1
		    if (GenPart_pdgId[k]>0){ //mu(-)
                      smeared_mean = 2.0*e/((1.0+A)-pow((1.0+A)*(1.0+A)+4.0*e*(-M-1.0/mean),0.5));
                    } else { //mu(+)
                      smeared_mean = 2.0*e/((1.0+A)-pow((1.0+A)*(1.0+A)+4.0*e*(M-1.0/mean),0.5));
                    }
		    smear_beta_weight_second_term = TMath::Gaus(secondPt_smear, smeared_mean, width) / TMath::Gaus(secondPt_smear, mean, width); 
		    secondPt_smear_beta_val = rans[nslot]->Gaus(smeared_mean, width);
		            
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
	      //--------------------------------------------------------------------------------
              
	      // ATTENTION OVERWRITING smear
              // smear mass directly
              //mll_smear = rans[nslot]->Gaus(mll_gen, 0.02*mll_gen);
              //smear_beta_weight = TMath::Gaus(mll_smear, mll_gen*(1-0.001), 0.02*mll_gen) / TMath::Gaus(mll_smear, mll_gen, 0.02*mll_gen);
              
	      //--------------------------------------------------------------------------------
                    

	      if(Muon_charge[i]==1){
		// <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val >
		temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco,mll_smear,firstPt_smear,secondPt_smear,mll_gen,firstPt_gen,secondPt_gen,smear_beta_weight,firstPt_smear_beta_val,secondPt_smear_beta_val);
	      } else {
		temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco,mll_smear,secondPt_smear,firstPt_smear,mll_gen,secondPt_gen,firstPt_gen,smear_beta_weight,secondPt_smear_beta_val,firstPt_smear_beta_val);
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
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    return pair_to_return;
  };

  auto d2 = d1.DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"});
    
  auto d3 = d2.Define("mll_reco","return get<2>(pairs);"); 
  auto d4 = d3.Filter("mll_reco>1.0"); // this means only events with one mu pair are kept 

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  //  <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>

  auto d5 = d4.Define("posTrackPt","return get<3>(pairs);")
    .Define("negTrackPt","return get<4>(pairs);")
    .Define("posTrackEta","return Muon_eta[get<0>(pairs)];")
    .Define("negTrackEta","return Muon_eta[get<1>(pairs)];")
    .Define("mll_smear","return get<5>(pairs);")
    .Define("posPtSmear","return get<6>(pairs);")
    .Define("negPtSmear","return get<7>(pairs);")
    .Define("mll_gen","return get<8>(pairs);")
    .Define("posPtGen","return get<9>(pairs);")
    .Define("negPtGen","return get<10>(pairs);")
    .Define("smear_beta_weight","return get<11>(pairs)*std::copysign(1.0, genWeight);")
    .Define("posPtSmearBetaVal","return get<12>(pairs);")
    .Define("negPtSmearBetaVal","return get<13>(pairs);")
    .Define("weight", "std::copysign(1.0, genWeight)");

auto d6 = d5.Define("mll_diff_reco_over_gen","return (mll_reco - mll_gen) / mll_gen;")
    .Define("mll_diff_reco","return mll_reco - mll_gen;")
    .Define("mll_diff_smear_over_gen","return (mll_smear - mll_gen) / mll_gen;")
    .Define("mll_diff_smear","return mll_smear - mll_gen;")
    .Define("mll_diff_smear_plus_offset","float offset = 0.1; return (mll_smear - mll_gen) + offset;") // offset goes here
    .Define("mll_diff_smear_minus_offset","float offset = 0.1; return (mll_smear - mll_gen) - offset;"); // offset goes here

auto d7 = d6.Define("jacobian_weight_mll_diff_smear", "return mll_diff_smear*weight;")   
    .Define("jacobian_weight_mll_diff_squared_smear","return mll_diff_smear*mll_diff_smear*weight;")
    .Define("jacobian_weight_mll_over_gen_all_squared_smear", "return (mll_smear/mll_gen)*(mll_smear/mll_gen)*weight;")
    .Define("jacobian_weight_mll_over_gen_smear", "return mll_smear/mll_gen*weight;")
    .Define("jacobian_weight_mll_over_gen_all_squared_reco", "return (mll_reco/mll_gen)*(mll_reco/mll_gen)*weight;")
    .Define("jacobian_weight_mll_over_gen_reco", "return mll_reco/mll_gen*weight;")
    // for mixed approach pdf
    .Define("jacobian_weight_mll_squared_smear","return mll_smear*mll_smear*weight;")
    .Define("jacobian_weight_mll_squared_gen","return mll_gen*mll_gen*weight;")
    .Define("jacobian_weight_mll_times_gen_smear","return mll_smear*mll_gen*weight;")
    // jacobians for unique pdf across k bins
    .Define("jacobian_weight_mll_diff_times_gen_smear", "return mll_diff_smear*mll_gen*weight;");
    
  //Save tree for debugging
  std::unique_ptr<TFile> f1( TFile::Open("snapshot_output.root", "RECREATE") );
  d7.Snapshot("Events", "snapshot_output.root", {"GenPart_status", "GenPart_pt", "posPtSmearBetaVal", "negPtSmearBetaVal", "Muon_charge", "GenPart_pdgId", "GenPart_genPartIdxMother"});
  f1->Close();

  //Control histograms
  auto mll_smear = d7.Histo1D({"mll_smear", "mll inclusive all bins", 20, 75.0, 105.0},"mll_smear","weight");
  auto mll_diff_smear = d7.Histo1D({"mll_diff_smear", "mll_diff inclusive all bins", 20, -5.0, 5.0},"mll_diff_smear","weight");
  auto mll_smear_beta_val = d7.Histo1D({"mll_smear_beta_val", "mll inclusive all bins", 20, 75.0, 105.0},"mll_smear","smear_beta_weight");
  auto pt_smear = d7.Histo1D({"pt_smear", "pt smear beta = 1", 15, 25.0, 55.0},"posPtSmear","weight");
  auto pt_smear_beta_val = d7.Histo1D({"pt_smear_beta_val", "pt smear beta = 0.001", 15, 25.0, 55.0},"posPtSmear","smear_beta_weight");

  std::unique_ptr<TFile> f2( TFile::Open("control_histo.root", "RECREATE") ); 
  mll_smear->Write();
  mll_diff_smear->Write();
  mll_smear_beta_val->Write();
  pt_smear->Write();
  pt_smear_beta_val->Write();
  f2->Close();
  
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff_over_gen=32, nbinsmll_diff=16, nbinsmll=32, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, ptbinranges, mllbinranges;
  vector<double> mll_diffbinranges, mll_diff_over_genbinranges;

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  double myetaboundaries[etabinranges.size()];
  for (int i=0; i<etabinranges.size(); i++){
    myetaboundaries[i] = etabinranges[i];
  }

  for (int i=0; i<=nbinsmll_diff; i++) mll_diffbinranges.push_back(-5.0 + i*10.0/nbinsmll_diff);
  for (int i=0; i<=nbinsmll; i++) mllbinranges.push_back(75.0 + i*(105.0-75.0)/nbinsmll);
  for (int i=0; i<=nbinsmll_diff_over_gen; i++) mll_diff_over_genbinranges.push_back(-0.1 + i*2.0*0.1/nbinsmll_diff_over_gen);

  std::unique_ptr<TFile> f3( TFile::Open("control_bin_histo.root", "RECREATE") );  

  //pT bin optimisation starts
  auto pt_pos_uni = d7.Histo1D({"pt_pos_uni", "pt mu+", nbinspt*3, ptlow, pthigh},"posTrackPt");
  f3->WriteObject(pt_pos_uni.GetPtr(), "pt_pos_uni");
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
  auto pt_pos = d7.Histo1D({"pt_pos", "pt mu+", nbinspt, myptboundaries},"posTrackPt");
  f3->WriteObject(pt_pos.GetPtr(), "pt_pos");

  auto pt_eta_pos = d7.Histo2D({"pt_eta_pos", "pt eta mu+", nbinseta, myetaboundaries, nbinspt, myptboundaries},"posTrackEta", "posTrackPt");
  f3->WriteObject(pt_eta_pos.GetPtr(), "pt_eta_pos");
  
  f3->Close();

  /*
  std::unique_ptr<TFile> f4( TFile::Open("multiD_histo_smear_beta_val_easy.root", "RECREATE") );
  // mll_diff_smear_beta_val_easy 
  auto mDh_diff_smear_beta_val_easy = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_beta_val_easy", "multi_data_histo_diff_smear_beta_val_easy", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","weight"});
  f4->WriteObject(mDh_diff_smear_beta_val_easy.GetPtr(), "multi_data_histo_diff_smear_beta_val_easy");
  f4->Close();
  */
  
  std::unique_ptr<TFile> f5( TFile::Open("multiD_histo_smear.root", "RECREATE") ); 
  
  //--------------------------------------------------------------------------------------
  // Smear histograms
  //--------------------------------------------------------------------------------------

  // Mass and mll_diff distributions
  
  // mll_smear
  auto mDh_smear = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_mll_smear", "multi_data_histo_mll_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","weight"});
  f5->WriteObject(mDh_smear.GetPtr(), "multi_data_histo_mll_smear");

  // mll_diff_smear_over_gen 
  auto mDh_diff_smear_over_gen = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_over_gen", "multi_data_histo_diff_smear_over_gen", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff_over_gen}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diff_over_genbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear_over_gen","weight"});
  f5->WriteObject(mDh_diff_smear_over_gen.GetPtr(), "multi_data_histo_diff_smear_over_gen");

  // mll_diff_smear
  auto mDh_diff_smear = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear", "multi_data_histo_diff_smear", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","weight"});
  f5->WriteObject(mDh_diff_smear.GetPtr(), "multi_data_histo_diff_smear");

  // mll_diff_smear_plus_offset
  auto mDh_diff_smear_plus_offset = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_plus_offset", "multi_data_histo_diff_smear_plus_offset", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear_plus_offset","weight"});
  f5->WriteObject(mDh_diff_smear_plus_offset.GetPtr(), "multi_data_histo_diff_smear_plus_offset");

  // mll_diff_smear_minus_offset
  auto mDh_diff_smear_minus_offset = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_minus_offset", "multi_data_histo_diff_smear_minus_offset", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear_minus_offset","weight"});
  f5->WriteObject(mDh_diff_smear_minus_offset.GetPtr(), "multi_data_histo_diff_smear_minus_offset");
 
  //--------------------------------------------------------------------------------------
  // Jacobian terms

  // mll_smear weighted by jacobian_weight_mll_over_gen_all_squared_smear
  auto mDh_jac_mll_over_gen_all_squared_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_over_gen_all_squared_smear_mll", "multi_data_histo_jac_mll_over_gen_all_squared_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_over_gen_all_squared_smear"});
  f5->WriteObject(mDh_jac_mll_over_gen_all_squared_smear_mll.GetPtr(), "multi_data_histo_jac_mll_over_gen_all_squared_smear_mll");  

  // mll_smear weighted by jacobian_weight_mll_over_gen_smear
  auto mDh_jac_mll_over_gen_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_over_gen_smear_mll", "multi_data_histo_jac_mll_over_gen_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_over_gen_smear"});
  f5->WriteObject(mDh_jac_mll_over_gen_smear_mll.GetPtr(), "multi_data_histo_jac_mll_over_gen_smear_mll");  

  // mll_smear weighted by jacobian_weight_mll_squared_smear
  auto mDh_jac_mll_squared_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_squared_smear_mll", "multi_data_histo_jac_mll_squared_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_squared_smear"});
  f5->WriteObject(mDh_jac_mll_squared_smear_mll.GetPtr(), "multi_data_histo_jac_mll_squared_smear_mll");

  // mll_smear weighted by jacobian_weight_mll_squared_gen
  auto mDh_jac_mll_gen_squared_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_gen_squared_smear_mll", "multi_data_histo_jac_mll_gen_squared_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_squared_gen"});
  f5->WriteObject(mDh_jac_mll_gen_squared_smear_mll.GetPtr(), "multi_data_histo_jac_mll_gen_squared_smear_mll");

  // mll_smear weighted by jacobian_weight_mll_times_gen_smear
  auto mDh_jac_mll_times_gen_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_times_gen_smear_mll", "multi_data_histo_jac_mll_times_gen_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_times_gen_smear"});
  f5->WriteObject(mDh_jac_mll_times_gen_smear_mll.GetPtr(), "multi_data_histo_jac_mll_times_gen_smear_mll");

  // mll_smear weighted by jacobian_weight_mll_diff_times_gen_smear
  auto mDh_jac_diff_times_gen_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_times_gen_smear_mll", "multi_data_histo_jac_diff_times_gen_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_diff_times_gen_smear"});
  f5->WriteObject(mDh_jac_diff_times_gen_smear_mll.GetPtr(), "multi_data_histo_jac_diff_times_gen_smear_mll");

  // mll_smear weighted by jacobian_weight_mll_diff_squared_smear
  auto mDh_jac_diff_squared_smear_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_squared_smear_mll", "multi_data_histo_jac_diff_squared_smear_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","jacobian_weight_mll_diff_squared_smear"});
  f5->WriteObject(mDh_jac_diff_squared_smear_mll.GetPtr(), "multi_data_histo_jac_diff_squared_smear_mll");

  // mll_diff_smear weighted by jacobian_weight_mll_diff_squared_smear
  auto mDh_jac_diff_squared_smear_mll_diff = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_squared_smear_mll_diff", "multi_data_histo_jac_diff_squared_smear_mll_diff", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","jacobian_weight_mll_diff_squared_smear"});
  f5->WriteObject(mDh_jac_diff_squared_smear_mll_diff.GetPtr(), "multi_data_histo_jac_diff_squared_smear_mll_diff");
  
  // mll_diff_smear weighted by jacobian_weight_mll_diff_smear
  auto mDh_jac_diff_smear_mll_diff = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_diff_smear_mll_diff", "multi_data_histo_jac_diff_smear_mll_diff", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","jacobian_weight_mll_diff_smear"});
  f5->WriteObject(mDh_jac_diff_smear_mll_diff.GetPtr(), "multi_data_histo_jac_diff_smear_mll_diff");
  
  //--------------------------------------------------------------------------------------
  // Weights to shift mass or mll_diff
  
  // mll_smear weighted by smear_beta_weight
  auto mDh_smear_beta_val = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_mll_smear_beta_val", "multi_data_histo_mll_smear_beta_val", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_smear","smear_beta_weight"});
  f5->WriteObject(mDh_smear_beta_val.GetPtr(), "multi_data_histo_mll_smear_beta_val");
  
  // mll_diff_smear weighted by smear_beta_weight
  auto mDh_diff_smear_beta_val = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_beta_val", "multi_data_histo_diff_smear_beta_val", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","smear_beta_weight"});
  f5->WriteObject(mDh_diff_smear_beta_val.GetPtr(), "multi_data_histo_diff_smear_beta_val");
  
  f5->Close();

  //--------------------------------------------------------------------------------------
  // Reco and gen histos
  //--------------------------------------------------------------------------------------
  
  std::unique_ptr<TFile> f6( TFile::Open("multiD_histo_reco.root", "RECREATE") );
  
  // mll_reco
  auto mDh_reco = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_mll_reco", "multi_data_histo_mll_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","weight"});
  f6->WriteObject(mDh_reco.GetPtr(), "multi_data_histo_mll_reco");

  // mll_diff_reco_over_gen 
  auto mDh_diff_reco_over_gen = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_reco_over_gen", "multi_data_histo_diff_reco_over_gen", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff_over_gen}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diff_over_genbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff_reco_over_gen","weight"});
  f6->WriteObject(mDh_diff_reco_over_gen.GetPtr(), "multi_data_histo_diff_reco_over_gen");
  
  // mll_diff_reco
  auto mDh_diff_reco = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_reco", "multi_data_histo_diff_reco", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff_reco","weight"});
  f6->WriteObject(mDh_diff_reco.GetPtr(), "multi_data_histo_diff_reco");
  
  //--------------------------------------------------------------------------------------
  // Jacobian terms

  // mll_reco weighted by jacobian_weight_mll_over_gen_all_squared_reco
  auto mDh_jac_mll_over_gen_all_squared_reco_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_over_gen_all_squared_reco_mll", "multi_data_histo_jac_mll_over_gen_all_squared_reco_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","jacobian_weight_mll_over_gen_all_squared_reco"});
  f6->WriteObject(mDh_jac_mll_over_gen_all_squared_reco_mll.GetPtr(), "multi_data_histo_jac_mll_over_gen_all_squared_reco_mll");   

  // mll_reco weighted by jacobian_weight_mll_over_gen_reco
  auto mDh_jac_mll_over_gen_reco_mll = d7.HistoND<float, float, float, float, float, double>({"multi_data_histo_jac_mll_over_gen_reco_mll", "multi_data_histo_jac_mll_over_gen_reco_mll", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mllbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_reco","jacobian_weight_mll_over_gen_reco"});
  f6->WriteObject(mDh_jac_mll_over_gen_reco_mll.GetPtr(), "multi_data_histo_jac_mll_over_gen_reco_mll"); 

  f6->Close();

  gettimeofday(&tv_stop, 0);
  t1 = (tv_stop.tv_sec - tv_start.tv_sec)*1000.0 + (tv_stop.tv_usec - tv_start.tv_usec)/1000.0;
  
  cout << "\n" << "# Calculation time: " << t1/60000.0 << " min" << "\n";  

  return 0;

}
  
