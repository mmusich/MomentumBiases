// Standalone code to fit for biases in mass distribtions per eta, pT, eta, pT bins

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TVectorT.h"
#include <string.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <sys/time.h>

using namespace ROOT;
using namespace ROOT::VecOps;
using ROOT::RDF::RNode;

using Eigen::MatrixXd;
using Eigen::VectorXd;

//--------------------------------------------------------
// Functions for TTree

//dxy_significance
RVecF dxy_significance(RVecF Muon_dxy, RVecF Muon_dxyErr){
  return abs(Muon_dxy)/Muon_dxyErr;
}

// MuonisGood
RVecB MuonisGood(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVec<UChar_t> Muon_genPartFlav, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
    if (Muon_pt[i] > 10.0 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && Muon_genPartFlav[i]==1 && dxy_significance[i] < 4){
      muonisgood.push_back(1);
    } else {
      muonisgood.push_back(0);
    }
  }
  return muonisgood;
}

// MuonisGoodData
RVecB MuonisGoodData(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
    if (Muon_pt[i] > 10.0 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
      muonisgood.push_back(1);
    } else {
      muonisgood.push_back(0);
    }
  }
  return muonisgood;
}

//--------------------------------------------------------
// Functions for histograms

// Functions for label names
string stringify_name(int a, int b, int c, int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

string stringify_title(int a, int b, int c, int d, vector<double> e, vector<double> p){
  string comma = ",";
  string txt1 = "#eta+ in [", txt2="] pT+ in [", txt3="] #eta- in [", txt4="] pT- in [", txt5 = "]";
  return txt1+to_string(e[a-1]).substr(0,4)+comma+to_string(e[a]).substr(0,4)+txt2+to_string(p[b-1]).substr(0,4)+comma+to_string(p[b]).substr(0,4)+txt3+to_string(e[c-1]).substr(0,4)+comma+to_string(e[c]).substr(0,4)+txt4+to_string(p[d-1]).substr(0,4)+comma+to_string(p[d]).substr(0,4)+txt5;
}

// Function for Gaussian fit and integral 
vector<double> fitHisto(TH1* histogram, int draw_option, int color, int nsigma){

  vector<double> fitresult;

  double mean = histogram->GetMean();
  double sigma = histogram->GetRMS();
  double mean_err, sigma_err, integral, chi2;

  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", mean - 5 * sigma, mean + 5 * sigma);
  // first fit
  if(0 == histogram->Fit(gaussianFunc, "QNR")){
    mean = gaussianFunc->GetParameter(1);
    sigma = gaussianFunc->GetParameter(2);

    // second fit: few sigma of first fit around mean of first fit
    gaussianFunc->SetRange(mean - nsigma * sigma, mean + nsigma * sigma);
    if (0 == histogram->Fit(gaussianFunc, "Q0R")) { // don't draw yet
      if (histogram->GetFunction(gaussianFunc->GetName())) { 
	histogram->GetFunction(gaussianFunc->GetName())->SetLineColor(color);
        if (draw_option == 1){ histogram->GetFunction(gaussianFunc->GetName())->ResetBit(TF1::kNotDraw);} // draw fit
      }
      mean = gaussianFunc->GetParameter(1);
      sigma = gaussianFunc->GetParameter(2);
      mean_err = gaussianFunc->GetParError(1);
      sigma_err = gaussianFunc->GetParError(2);
      chi2 = gaussianFunc->GetChisquare();

      TF1 *gaussianNewFunc = new TF1("gaussianNewFunc", "gaus", mean - 5 * sigma, mean + 5 * sigma);
      gaussianNewFunc->SetParameter(0, 1.0/(sigma*pow(2*M_PI,0.5)));
      gaussianNewFunc->SetParameter(1, mean);
      gaussianNewFunc->SetParameter(2, sigma);

      // integral of fitted gaussian function in range of the histogram
      integral = gaussianNewFunc->Integral(75.0, 105.0);
      // TODO implement line below
      //integral = gaussianNewFunc->Integral(histogram->GetXaxis()->GetBinLowEdge(1), histogram->GetXaxis()->GetBinLowEdge(histogram->GetNbins()+1)); // TODO histogram->GetNbins() includes under/over flow??

    } else {
      mean = -90.0;
      sigma = -5.0;
      mean_err = 90.0;
      sigma_err = 5.0;
      integral = -100.0;
      chi2 = -100.0;
    } 
  } else {
    mean = -90.0;
    sigma = -5.0;
    mean_err = 90.0;
    sigma_err = 5.0;
    integral = -100.0;
    chi2 = -100.0;
  }
  
  fitresult.push_back(mean);
  fitresult.push_back(sigma);
  fitresult.push_back(mean_err);
  fitresult.push_back(sigma_err);
  fitresult.push_back(integral);
  fitresult.push_back(chi2);
  
  return fitresult;
}

//--------------------------------------------------------
// Main
int bias_fitter(){

  double t1(0.);
  // Get start time
  struct timeval tv_start, tv_stop;
  gettimeofday(&tv_start, 0);

  ROOT::EnableImplicitMT();

  // Choose validation/analysis mode
  //string mode_option("simulation"); //TODO pass as command line argument

  //TODO make it so that label of data sample is attached to output files

  //--------------------------------------------------------
  // Define dataframe
  
  TChain chain("Events");
  chain.Add("/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_1.root");
  chain.Add("/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0000/NanoV9MCPostVFP_1.root");
  //chain.Add("/scratch/wmass/y2017/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2017_TrackFitV722_NanoProdv3/NanoV9MC2017_1.root");
  /*
  string line;
  
  ifstream file("InOutputFiles/MCFilenames_2016.txt");
  if (file.is_open()) {
    while (getline(file, line)) {
      chain.Add(line.c_str());
    }
    file.close();
  }
  */
  RDataFrame df(chain); // TODO write for analysis mode, will need 2 data frames, but honestly should be a different script that people use

  // Random numbers for pairs
  unsigned int nslots = df.GetNSlots();
  std::cout<<"nslots"<<nslots<<"\n";
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(4357 + i*10) );
  }
  
  auto dlast = std::make_unique<RNode>(df);

  dlast = std::make_unique<RNode>(dlast->Filter("HLT_IsoMu24 == 1"));
  dlast = std::make_unique<RNode>(dlast->Filter("nMuon >= 2"));
  dlast = std::make_unique<RNode>(dlast->Filter("PV_npvsGood >= 1"));
  
  //if (mode_option.compare("simulation") == 0) {
  dlast = std::make_unique<RNode>(dlast->Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)"));
  dlast = std::make_unique<RNode>(dlast->Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, Muon_genPartFlav, dxy_significance)"));
    //} else if (mode_option.compare("data") == 0){
    //dlast = std::make_unique<RNode>(dlast->Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)"));
    //dlast = std::make_unique<RNode>(dlast->Define("MuonisGood", "MuonisGoodData(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)"));
    //}
  
  // Binning (pT done later)
  double pt_low = 25.0, pt_high = 55.0, eta_low = -2.4, eta_high = 2.4, mll_diff_low = -7.0,  mll_diff_high = 7.0, mll_low = 75.0, mll_high = 105.0;
  int nbinsmll_diff_over_gen=32, nbinsmll_diff=22, nbinsmll=32, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, ptbinranges, mllbinranges;
  vector<double> mll_diffbinranges, mll_diff_over_genbinranges;
  vector<float> A_values(nbinseta), e_values(nbinseta), M_values(nbinseta);

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){
    etabinranges.push_back(eta_low + i * (eta_high - eta_low)/nbinseta);
    std::cout<<etabinranges[i]<<", ";
  }
  std::cout<<"] \n";

  double myetaboundaries[etabinranges.size()];
  for (int i=0; i<etabinranges.size(); i++){
    myetaboundaries[i] = etabinranges[i];
  }

  for (int i=0; i<=nbinsmll_diff; i++) mll_diffbinranges.push_back(mll_diff_low + i*(mll_diff_high - mll_diff_low)/nbinsmll_diff);
  double mymll_diffboundaries[mll_diffbinranges.size()];
  for (int i=0; i<mll_diffbinranges.size(); i++){
    mymll_diffboundaries[i] = mll_diffbinranges[i];
  }
  
  for (int i=0; i<=nbinsmll; i++) mllbinranges.push_back(mll_low + i*(mll_high-mll_low)/nbinsmll);
  for (int i=0; i<=nbinsmll_diff_over_gen; i++) mll_diff_over_genbinranges.push_back(-0.1 + i*2.0*0.1/nbinsmll_diff_over_gen);

  int nbins = nbinseta*nbinseta*nbinspt*nbinspt;
  int nbins_plus_one = nbins+1;
  double mybinsboundaries[nbins_plus_one];
  for (int i=0; i<=nbins; i++){
    mybinsboundaries[i] = i;
  }
  
  for (int i=0; i<nbinseta; i++){
    //std::cout<< "\n" << "i: " << i << "\n";
    // A from -0.0002 to 0.0005 and back
    A_values[i] = (-7.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 5.0)*0.0001;
    //A_values[i] = 0.0004; 
    //std::cout<< "A: " << A_values[i] << ", ";
    // e from 0.01 to 0.001 and back
    e_values[i] = (9.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 1.0)*0.001;
    //e_values[i] = 0.002; 
    //std::cout<< "e: " << e_values[i] << ", ";
    // M from 4*10^-5 to -2*10^-5 and back
    M_values[i] = (6.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) - 2.0)*0.00001;
    //M_values[i] = 0.00001; 
    //std::cout<< "M: " << M_values[i] << ", ";
  }

  auto GetEtaBin = [&](float eta)->int{
    for (int i=0; i<nbinseta; i++){
      if( etabinranges[i] <= eta && eta < etabinranges[i+1]) { return i; }
    }
    return -9;
  };

  //Pairs
  auto pairs = [&](unsigned int nslot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi)->std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float>{

    // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
    RVec<std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float>> pairs; 
    std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco, firstPt_smear, secondPt_smear, mll_smear, firstPt_gen, secondPt_gen, mll_gen, firstPt_smear_beta_val, secondPt_smear_beta_val, smear_beta_weight;
    float smear_pt, mean, width, smeared_mean, smeared_curvature, smear_beta_weight_first_term, smear_beta_weight_second_term;
    //float A = 0.0004, e = 0.002, M = 0.00001; //for now constant per eta bin
    float A, e, M;
    int eta_bin_idx;

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
	            
	    if(mll_low<mll_reco && mll_reco<mll_high){ //Cut in mll
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
		    eta_bin_idx = GetEtaBin(Muon_eta[i]);
		    if (eta_bin_idx < 0 || eta_bin_idx > 23){
		      std::cout<<"\n"<<"WARNING eta out of bounds";
		      pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
		      return pair_to_return ;
		    }
		    A = A_values[eta_bin_idx];
		    e = e_values[eta_bin_idx];
		    M = M_values[eta_bin_idx];
		    if (GenPart_pdgId[k]>0){ //mu(-)
		      if ( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) >= 0.0 ){
			smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) ,0.5) )/(-2.0)/e;
			if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
			  std::cout<<"\n"<<"WARNING WRONG CURVATURE";
			  pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
			  return pair_to_return ;
			}
		      } else {
			std::cout<<"\n"<<"WARNING negative delta";
			pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
			return pair_to_return ;
		      }
		      //smeared_curvature = 1.0/mean + M; // for M only
		    } else { //mu(+)
		      if ( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) >= 0.0 ){
			smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) ,0.5) )/(-2.0)/e;
			if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
			  std::cout<<"\n"<<"WARNING WRONG CURVATURE";
			  pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
			  return pair_to_return ;
			}
                      } else {
                        std::cout<<"\n"<<"WARNING negative delta";
                        pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                        return pair_to_return ;
                      }
		      //smeared_curvature = 1.0/mean - M; // for M only
		    }
		    smeared_mean = 1.0/smeared_curvature;
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
		    eta_bin_idx = GetEtaBin(Muon_eta[j]);
		    if (eta_bin_idx < 0 || eta_bin_idx > 23){
                      std::cout<<"\n"<<"WARNING eta out of bounds";
                      pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                      return pair_to_return ;
                    }
                    A = A_values[eta_bin_idx];
                    e = e_values[eta_bin_idx];
                    M = M_values[eta_bin_idx];
		    if (GenPart_pdgId[k]>0){ //mu(-)
                      if ( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) >= 0.0 ){
                        smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) ,0.5) )/(-2.0)/e;
                        if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
                          std::cout<<"\n"<<"WARNING WRONG CURVATURE";
                          pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                          return pair_to_return ;
                        }
                      } else {
                        std::cout<<"\n"<<"WARNING negative delta";
                        pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                        return pair_to_return ;
                      }
		      //smeared_curvature = 1.0/mean + M; // for M only
                    } else { //mu(+)
                      if ( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) >= 0.0 ){
                        smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) ,0.5) )/(-2.0)/e;
                        if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
                          std::cout<<"\n"<<"WARNING WRONG CURVATURE";
                          pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                          return pair_to_return ;
                        }
                      } else {
                        std::cout<<"\n"<<"WARNING negative delta";
                        pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
                        return pair_to_return ;
                      }
		      //smeared_curvature = 1.0/mean - M; // for M only
                    }
		    smeared_mean = 1.0/smeared_curvature;
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
	      //if(std::isnan(smear_beta_weight)){
	      //std::cout << "slot" << nslot << ": smear_beta_weight =" << smear_beta_weight << "; first smeared_mean =" << firstPt_smear_beta_val << "; first mean = " << firstPt_smear << "; second smeared_mean =" << secondPt_smear_beta_val << "; second mean = " << secondPt_smear << "\n";
	      //}
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

  /*
  //PairsData
  auto pairsData = [&](RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz)->std::tuple<int,int,float,float,float>{
    
    // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco>
    RVec<std::tuple<int,int,float,float,float>> pairs; 
    std::tuple<int,int,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco;
    
    for(int i=1;i<Muon_pt.size();i++){
      if(MuonisGood[i]){
	for(int j=0;j<i;j++){
	  if(MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){ //TODO better deltaR matching
	    TLorentzVector firstTrack, secondTrack, mother;
	    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	    mother = firstTrack + secondTrack;
	    mll_reco = mother.M();
	    firstPt_reco = Muon_pt[i];
	    secondPt_reco = Muon_pt[j];
	            
	    if(mll_low<mll_reco && mll_reco<mll_high){ //Cut in mll
	      if(Muon_charge[i]==1){
		// <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco >
		temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco);
	      } else {
		temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco);
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
  */
  
  //if (mode_option.compare("simulation") == 0) {
  dlast = std::make_unique<RNode>(dlast->DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"}));
    //} else if (mode_option.compare("data") == 0){
    //dlast = std::make_unique<RNode>(dlast->Define("pairs", pairsData, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz"}));
    //}
  
  dlast = std::make_unique<RNode>(dlast->Define("mll_reco","return get<2>(pairs);")); 
  dlast = std::make_unique<RNode>(dlast->Filter("mll_reco>1.0")); // this means only events with one mu pair are kept
    
  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
  
  dlast = std::make_unique<RNode>(dlast->Define("posTrackPt","return get<3>(pairs);")); //TOOD make binning on smear not reco
  dlast = std::make_unique<RNode>(dlast->Define("negTrackPt","return get<4>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Define("posTrackEta","return Muon_eta[get<0>(pairs)];"));
  dlast = std::make_unique<RNode>(dlast->Define("negTrackEta","return Muon_eta[get<1>(pairs)];"));

  //if (mode_option.compare("simulation") == 0) { //NOTE ATTENTION can't yet run analysis mode with optimised pT binning
  dlast = std::make_unique<RNode>(dlast->Define("weight", "std::copysign(1.0, genWeight)"));
  //}
  
  //pT bin optimisation starts
  auto pt_pos_uni = dlast->Histo1D({"pt_pos_uni", "pt mu+", nbinspt*3, pt_low, pt_high},"posTrackPt","weight");
  
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

  //---------------------------------------------------
  // Sort events in eta pt eta pt bins

  // GetEtaPtEtaPtBin
  auto GetEtaPtEtaPtBin = [&](float pos_eta, float pos_pt, float neg_eta, float neg_pt )->int{
    int idx = 0;
    for (int i=0; i<nbinseta; i++){
      for (int j=0; j<nbinspt; j++){
	for (int k=0; k<nbinseta; k++){
	  for (int l=0; l<nbinspt; l++){
	    if(etabinranges[i] <= pos_eta && pos_eta < etabinranges[i+1] && ptbinranges[j] <= pos_pt && pos_pt < ptbinranges[j+1] && etabinranges[k] <= neg_eta && neg_eta < etabinranges[k+1] && ptbinranges[l] <= neg_pt && neg_pt < ptbinranges[l+1] ) { return idx; } 
	    idx++;
	  }
	}
      }
    }
    std::cout<<"WARNING miscounting k bins"<<"\n";
    return -9;
  };

  //dlast = std::make_unique<RNode>(dlast->Define("bin_index_reco", GetEtaPtEtaPtBin ,{"posTrackEta", "posTrackPt", "negTrackEta", "negTrackPt"}));
  
  auto d_mc = std::make_unique<RNode>(df);
  auto d_sim_data = std::make_unique<RNode>(df);  
  
  //if (mode_option.compare("simulation") == 0) {
  
  dlast = std::make_unique<RNode>(dlast->Define("mll_gen","return get<8>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Define("mll_diff_reco","return mll_reco - mll_gen;"));
  dlast = std::make_unique<RNode>(dlast->Define("mll_smear","return get<5>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Define("mll_diff_smear","return mll_smear - mll_gen;"));
  dlast = std::make_unique<RNode>(dlast->Define("posPtSmear","return get<6>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Filter("posPtSmear >= 25.0 && posPtSmear < 55.0")); // TODO write a filter function with [&] to use pt_high rather than value directly
  dlast = std::make_unique<RNode>(dlast->Define("negPtSmear","return get<7>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Filter("negPtSmear >= 25.0 && negPtSmear < 55.0"));
  dlast = std::make_unique<RNode>(dlast->Define("smear_beta_weight","return get<11>(pairs)*weight;"));
  dlast = std::make_unique<RNode>(dlast->Define("posPtSmearBetaVal","return get<12>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Filter("posPtSmearBetaVal >= 25.0 && posPtSmearBetaVal < 55.0"));
  dlast = std::make_unique<RNode>(dlast->Define("negPtSmearBetaVal","return get<13>(pairs);"));
  dlast = std::make_unique<RNode>(dlast->Filter("negPtSmearBetaVal >= 25.0 && negPtSmearBetaVal < 55.0"));
  
  dlast = std::make_unique<RNode>(dlast->Define("bin_index_smear", GetEtaPtEtaPtBin ,{"posTrackEta", "posPtSmear", "negTrackEta", "negPtSmear"}));
  dlast = std::make_unique<RNode>(dlast->Define("bin_index_smear_beta_val", GetEtaPtEtaPtBin ,{"posTrackEta", "posPtSmearBetaVal", "negTrackEta", "negPtSmearBetaVal"}));
    
  // define a frame that plays the role of data -> odd events, MC -> even events
  d_sim_data = std::make_unique<RNode>(dlast->Filter("event%2==1"));
  d_mc = std::make_unique<RNode>(dlast->Filter("event%2==0"));

  // Jacobians TODO if you don't save individual jac terms, save the entire jacobian (4 validation, 2 ana) in pairs fct
  // Smear
  d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_diff_smear", "return mll_diff_smear*weight;"));
  d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_diff_squared_smear","return mll_diff_smear*mll_diff_smear*weight;"));
  d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_minus_2gen_smear","return (mll_smear-2.0*mll_gen)*weight;"));
  d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_diff_times_gen_smear", "return mll_diff_smear*mll_gen*weight;"));
  // Reco
  //d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_diff_reco", "return mll_diff_reco*weight;"));
  //d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_diff_squared_reco","return mll_diff_reco*mll_diff_reco*weight;"));
  //d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_minus_2gen_reco","return (mll_reco-2.0*mll_gen)*weight;"));
  //d_mc = std::make_unique<RNode>(d_mc->Define("jacobian_weight_mll_diff_times_gen_reco", "return mll_diff_reco*mll_gen*weight;"));
  
  // save 2d histo of mll diff, mass, jacobian terms vs bin_index_(smear/smear_beta_val)
  //NOTE in analysis mode we only need mass jacobians, so 4 2D histos, so it might be worth it saving jac terms over running over the TTree again

  // binned in mll_diff
  // MC
  auto mll_diff_smear_diffbin_4Dbin = d_mc->Histo2D({"mll_diff_smear_diffbin_4Dbin", "mll_diff_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"mll_diff_smear", "bin_index_smear", "weight");
  // Jac
  auto jac_weight_mll_diff_smear_diffbin_4Dbin = d_mc->Histo2D({"jac_weight_mll_diff_smear_diffbin_4Dbin", "jac_weight_mll_diff_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"jacobian_weight_mll_diff_smear", "bin_index_smear");
  auto jac_weight_mll_diff_squared_smear_diffbin_4Dbin = d_mc->Histo2D({"jac_weight_mll_diff_squared_smear_diffbin_4Dbin", "jac_weight_mll_diff_squared_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"jacobian_weight_mll_diff_squared_smear", "bin_index_smear");
  // Dummy data
  auto mll_diff_smear_beta_val_diffbin_4Dbin = d_sim_data->Histo2D({"mll_diff_smear_beta_val_diffbin_4Dbin", "mll_diff_smear_beta_val_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"mll_diff_smear", "bin_index_smear_beta_val", "smear_beta_weight");
  
  // binned in mll
  // MC
  auto mll_smear_mllbin_4Dbin = d_mc->Histo2D({"mll_smear_mllbin_4Dbin", "mll_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"mll_smear", "bin_index_smear", "weight");
  // Jac
  auto jac_weight_mll_diff_smear_mllbin_4Dbin = d_mc->Histo2D({"jac_weight_mll_diff_smear_mllbin_4Dbin", "jac_weight_mll_diff_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"jacobian_weight_mll_diff_smear", "bin_index_smear");
  auto jac_weight_mll_diff_squared_smear_mllbin_4Dbin = d_mc->Histo2D({"jac_weight_mll_diff_squared_smear_mllbin_4Dbin", "jac_weight_mll_diff_squared_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"jacobian_weight_mll_diff_squared_smear", "bin_index_smear");
  auto jac_weight_mll_minus_2gen_smear_mllbin_4Dbin = d_mc->Histo2D({"jac_weight_mll_minus_2gen_smear_mllbin_4Dbin", "jac_weight_mll_minus_2gen_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"jacobian_weight_mll_minus_2gen_smear", "bin_index_smear");
  auto jac_weight_mll_diff_times_gen_smear_mllbin_4Dbin = d_mc->Histo2D({"jac_weight_mll_diff_times_gen_smear_mllbin_4Dbin", "jac_weight_mll_diff_times_gen_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"jacobian_weight_mll_diff_times_gen_smear", "bin_index_smear");
  // Dummy data
  auto mll_smear_beta_val_mllbin_4Dbin = d_sim_data->Histo2D({"mll_smear_beta_val_mllbin_4Dbin", "mll_smear_beta_val_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"mll_smear", "bin_index_smear_beta_val", "smear_beta_weight");
  
  std::unique_ptr<TFile> fz( TFile::Open("InOutputFiles/test.root", "RECREATE") );
  mll_diff_smear_diffbin_4Dbin->Write();
  mll_smear_mllbin_4Dbin->Write();
  fz->Close();
    
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Mass fitting part
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Variables declaration 

  //TODO comment explain where the variables are needed
  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0; // to check number of passing events
  double value_alpha, error_alpha, value_beta, error_beta, value_epsilon, error_epsilon, diff_squared, diff, evts_in_bin, error_evts_in_bin, error_diff_squared, error_diff, mll_minus_2gen, error_mll_minus_2gen, diff_times_gen, error_diff_times_gen, jac_b_weight, error_jac_b_weight; // for jac calculation
  double mean_mc = 0.0, error_mean_mc = 0.0, sigma_mc=0.0, error_sigma_mc = 0.0, integral_mc, chi2_mc;
  double mean_diff_data = 1.0, error_mean_diff_data = 1.0, sigma_diff_data = 1.0, error_sigma_diff_data = 1.0, integral_diff_data = 1.0; //TODO just initialise with 0, are they called in analysis mode at all?
  double fit_beta_error, value_corrected, error_corrected, mean_corrected_diff, error_mean_corrected_diff, sigma_corrected_diff, error_sigma_corrected_diff;
  double max_hist_mll_diff, max_hist_mll;
  int filled_bins_mll, filled_bins_mll_diff, position_to_fill;
  string name, leg_entry, title;
  vector<double> fitresult;
  
  auto mDh_proj_diff_mc = mll_diff_smear_diffbin_4Dbin->ProjectionX("mDh_proj_diff_mc", 0, 1, "e");
  auto mDh_proj_mll_mc = mll_smear_mllbin_4Dbin->ProjectionX("mDh_proj_mll_mc", 0, 1, "e");
  auto mDh_proj_jac_diff_mll = jac_weight_mll_diff_smear_mllbin_4Dbin->ProjectionX("mDh_proj_jac_diff_mll", 0, 1, "e");
  auto mDh_proj_jac_mll_minus_2gen_mll = jac_weight_mll_minus_2gen_smear_mllbin_4Dbin->ProjectionX("mDh_proj_jac_mll_minus_2gen_mll", 0, 1, "e");
  auto mDh_proj_jac_diff_times_gen_mll = jac_weight_mll_diff_times_gen_smear_mllbin_4Dbin->ProjectionX("mDh_proj_jac_diff_times_gen_mll", 0, 1, "e");
  auto mDh_proj_jac_diff_squared_mll = jac_weight_mll_diff_squared_smear_mllbin_4Dbin->ProjectionX("mDh_proj_jac_diff_squared_mll", 0, 1, "e");
  auto mDh_proj_mll_data = mll_smear_beta_val_mllbin_4Dbin->ProjectionX("mDh_proj_mll_data", 0, 1, "e");
    
  auto mDh_proj_jac_diff_squared_mll_diff = jac_weight_mll_diff_squared_smear_diffbin_4Dbin->ProjectionX("mDh_proj_jac_diff_squared_mll_diff", 0, 1, "e");
  auto mDh_proj_jac_diff_mll_diff = jac_weight_mll_diff_smear_diffbin_4Dbin->ProjectionX("mDh_proj_jac_diff_mll_diff", 0, 1, "e");
  auto mDh_proj_diff_data = mll_diff_smear_beta_val_diffbin_4Dbin->ProjectionX("mDh_proj_diff_data", 0, 1, "e");
    
  // for each projection for bin_idx check stats, fit, save mean rms
  for (int i=1; i<10; i++){ // TODOD i<nbins

    all_histos_count++;

    // Require enough stats in diff_mc for sigma_MC fit
    delete gROOT->FindObject("mDh_proj_diff_mc");
    mDh_proj_diff_mc = mll_diff_smear_diffbin_4Dbin->ProjectionX("mDh_proj_diff_mc", i, i+1, "e");
    nevents = mDh_proj_diff_mc->Integral(1,nbinsmll_diff); 
    
    if (nevents < 300.0){ // reject low stats
      empty_histos_count++;
    } else {

      // Require most of mll_mc Gaussian to be in fit window
      // TODO and require same from mll_data
      delete gROOT->FindObject( "mDh_proj_mll_mc" );
      mDh_proj_mll_mc = mll_smear_mllbin_4Dbin->ProjectionX("mDh_proj_mll_mc", i, i+1, "e");
      mDh_proj_mll_mc->GetXaxis()->SetTitle("mll [GeV]");
      mDh_proj_mll_mc->GetYaxis()->SetTitle("Events");
      fitresult = fitHisto(mDh_proj_mll_mc, 0, 8, 5);
      
      if (fitresult[4] < 0.75){ // reject small gaus integral 
	empty_histos_count++;
      } else {
	// do mlldiff fits, save mean, rms, in a histo

	// do mlldiff fits, save alpha, beta, nu in a histo
	
	// do mass fits, save alpha, beta, nu in a histo
      }
      
    }
  }
  
  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Dataframe control histograms
  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  /*
      //Save tree for debugging
      std::unique_ptr<TFile> f1( TFile::Open("InOutputFiles/snapshot_output.root", "RECREATE") );
      dlast->Snapshot("Events", "snapshot_output.root", {"GenPart_status", "GenPart_pt", "posPtSmearBetaVal", "negPtSmearBetaVal", "Muon_charge", "GenPart_pdgId", "GenPart_genPartIdxMother"});
      f1->Close();
      
      //Control histograms
      auto mll_smear = d_mc->Histo1D({"mll_smear", "mll inclusive all bins", 20, mll_low, mll_high},"mll_smear","weight");
      auto mll_diff_smear = d_mc->Histo1D({"mll_diff_smear", "mll_diff inclusive all bins", 20, mll_diff_low, mll_diff_high},"mll_diff_smear","weight");
      auto mll_smear_beta_val = d_sim_data->Histo1D({"mll_smear_beta_val", "mll inclusive all bins #beta != 1", 20, mll_low, mll_high},"mll_smear","smear_beta_weight");
      auto pt_smear = d_mc->Histo1D({"pt_smear", "pt smear beta = 1", 15, pt_low, pt_high},"posPtSmear","weight");
      auto pt_smear_beta_val = d_sim_data->Histo1D({"pt_smear_beta_val", "pt smear #beta != 1", 15, pt_low, pt_high},"posPtSmear","smear_beta_weight");
      
      std::unique_ptr<TFile> f2( TFile::Open("InOutputFiles/control_histo.root", "RECREATE") );
      mll_smear->Write();
      mll_diff_smear->Write();
      mll_smear_beta_val->Write();
      pt_smear->Write();
      pt_smear_beta_val->Write();
      f2->Close();
      
      std::unique_ptr<TFile> f3( TFile::Open("InOutputFiles/control_bin_histo.root", "RECREATE") );
      f3->WriteObject(pt_pos_uni.GetPtr(), "pt_pos_uni");
      
      //TH1 in pT+ with variable bin size -> should be uniform
      auto pt_pos = d_mc->Histo1D({"pt_pos", "pt mu+", nbinspt, myptboundaries},"posTrackPt","weight");
      f3->WriteObject(pt_pos.GetPtr(), "pt_pos");
      
      auto pt_eta_pos = d_mc->Histo2D({"pt_eta_pos", "pt eta mu+", nbinseta, myetaboundaries, nbinspt, myptboundaries},"posTrackEta", "posTrackPt", "weight");
      f3->WriteObject(pt_eta_pos.GetPtr(), "pt_eta_pos");
      
      // small uniform pt binning
      vector<double> smallptunibinning;
      for (int i=0; i<=nbinspt*10; i++){
      smallptunibinning.push_back(pt_low + i * (pt_high - pt_low)/nbinspt/10);
      }
      
      double mysmallptunibinning[smallptunibinning.size()];
      for (int i=0; i<smallptunibinning.size(); i++){
      mysmallptunibinning[i] = smallptunibinning[i];
      }
      
      auto pt_pos_bin = d_mc->Histo2D({"pt_pos_bin", "pt mu+ distr in bin", nbinspt*10, mysmallptunibinning, nbinspt, myptboundaries}, "posTrackPt", "posTrackPt", "weight");
      f3->WriteObject(pt_pos_bin.GetPtr(), "pt_pos_bin");
      
      auto pt_neg_bin = d_mc->Histo2D({"pt_neg_bin", "pt mu- distr in bin", nbinspt*10, mysmallptunibinning, nbinspt, myptboundaries}, "negTrackPt", "negTrackPt", "weight");
      f3->WriteObject(pt_neg_bin.GetPtr(), "pt_neg_bin");
      
      f3->Close();
  */
  /*
      std::unique_ptr<TFile> f4( TFile::Open("InOutputFiles/multiD_histo_smear_beta_val_easy.root", "RECREATE") );
      // mll_diff_smear_beta_val_easy 
      auto mDh_diff_smear_beta_val_easy = d_sim_data->HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_beta_val_easy", "multi_data_histo_diff_smear_beta_val_easy", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","weight"});
      f4->WriteObject(mDh_diff_smear_beta_val_easy.GetPtr(), "multi_data_histo_diff_smear_beta_val_easy");
      f4->Close();
  */
  
  //} closes if (mode_option.compare("simulation") == 0) {
  
  gettimeofday(&tv_stop, 0);
  t1 = (tv_stop.tv_sec - tv_start.tv_sec)*1000.0 + (tv_stop.tv_usec - tv_start.tv_usec)/1000.0;
  
  cout << "\n" << "# Calculation time: " << t1/60000.0 << " min" << "\n";  
  
  return 0;

}
  
