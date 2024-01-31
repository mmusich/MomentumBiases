// Standalone code to fit for sagitta bias corrections

#include "TFile.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TVectorT.h"
#include <string.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace ROOT;
using namespace ROOT::VecOps;

string stringify_name(int a, int b, int c, int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

string stringify_title(int a, int b, int c, int d, vector<double> e, vector<double> p){
  string comma = ",";
  string txt1 = "#eta+ in [", txt2="] pT+ in [", txt3="] #eta- in [", txt4="] pT- in [", txt5 = "]";
  return txt1+to_string(e[a-1]).substr(0,4)+comma+to_string(e[a]).substr(0,4)+txt2+to_string(p[b-1]).substr(0,4)+comma+to_string(p[b]).substr(0,4)+txt3+to_string(e[c-1]).substr(0,4)+comma+to_string(e[c]).substr(0,4)+txt4+to_string(p[d-1]).substr(0,4)+comma+to_string(p[d]).substr(0,4)+txt5;
}

/////////////// Function for Gaussian fit and integral ////////////////////////////////////////
vector<double> fitHisto(TH1* histogram, int option, int color){

  vector<double> fitresult;

  double mean = histogram->GetMean();
  double sigma = histogram->GetRMS();
  double mean_err, sigma_err, integral;

  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", mean - 3 * sigma, mean + 3 * sigma); 
  if(0 == histogram->Fit(gaussianFunc, "QNR")){
    mean = gaussianFunc->GetParameter(1);
    sigma = gaussianFunc->GetParameter(2);

    // second fit: two sigma of first fit around mean of first fit
    gaussianFunc->SetRange(mean - 5 * sigma, mean + 5 * sigma);
    if (0 == histogram->Fit(gaussianFunc, "Q0R")) {
      if (histogram->GetFunction(gaussianFunc->GetName())) { // Take care that it is later on drawn:
	histogram->GetFunction(gaussianFunc->GetName())->SetLineColor(color);
        if (option == 1){ histogram->GetFunction(gaussianFunc->GetName())->ResetBit(TF1::kNotDraw);}
      }
      mean = gaussianFunc->GetParameter(1);
      sigma = gaussianFunc->GetParameter(2);
      mean_err = gaussianFunc->GetParError(1);
      sigma_err = gaussianFunc->GetParError(2);

      TF1 *gaussianNewFunc = new TF1("gaussianNewFunc", "gaus", mean - 5 * sigma, mean + 5 * sigma);
      gaussianNewFunc->SetParameter(0, 1.0/(sigma*pow(2*M_PI,0.5)));
      gaussianNewFunc->SetParameter(1, mean);
      gaussianNewFunc->SetParameter(2, sigma);

      integral = gaussianNewFunc->Integral(75.0, 105.0);

    } else {
      mean = -90.0;
      sigma = -5.0;
      mean_err = 90.0;
      sigma_err = 5.0;
      integral = -100.0;
    } 
  } else {
    mean = -90.0;
    sigma = -5.0;
    mean_err = 90.0;
    sigma_err = 5.0;
    integral = -100.0;
  }

  fitresult.push_back(mean);
  fitresult.push_back(sigma);
  fitresult.push_back(mean_err);
  fitresult.push_back(sigma_err);
  fitresult.push_back(integral);

  return fitresult;
}

int frame(){

  ROOT::EnableImplicitMT(128);

  /////////////// Input files /////////////////////////////////////////////////////////

  //reco
  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo_reco.root") );
  std::unique_ptr<THnD> mDh_reco(myFile->Get<THnD>("multi_data_histo_reco"));
  std::unique_ptr<THnD> mDh_gen(myFile->Get<THnD>("multi_data_histo_gen"));
  std::unique_ptr<THnD> mDh_diff_reco(myFile->Get<THnD>("multi_data_histo_diff_reco")); // it's reco - gen

  //smear
  std::unique_ptr<TFile> myFile2( TFile::Open("multiD_histo_smear.root") );
  std::unique_ptr<THnD> mDh_smear(myFile2->Get<THnD>("multi_data_histo_smear"));
  std::unique_ptr<THnD> mDh_diff_smear(myFile2->Get<THnD>("multi_data_histo_diff_smear")); //it's smeared - gen
  std::unique_ptr<THnD> mDh_smear_beta_val(myFile2->Get<THnD>("multi_data_histo_smear_beta_val"));
  std::unique_ptr<THnD> mDh_diff_smear_beta_val(myFile2->Get<THnD>("multi_data_histo_diff_smear_beta_val")); //it's smeared - gen
  //weights
  std::unique_ptr<THnD> mDh_diff_squared_smear(myFile2->Get<THnD>("multi_data_histo_diff_squared_smear")); // it's smeared - gen  
  std::unique_ptr<THnD> mDh_diff_squared_smear_control(myFile2->Get<THnD>("multi_data_histo_diff_squared_smear_control")); // it's smeared - gen
  std::unique_ptr<THnD> mDh_jac_beta_smear(myFile2->Get<THnD>("multi_data_histo_jac_beta_smear")); 
  std::unique_ptr<THnD> mDh_jac_beta_smear_mll_diff(myFile2->Get<THnD>("multi_data_histo_jac_beta_smear_mll_diff"));

  //smear easy way
  std::unique_ptr<TFile> myFile3( TFile::Open("multiD_histo_smear_beta_val_easy.root") );
  std::unique_ptr<THnD> mDh_diff_smear_beta_val_easy(myFile3->Get<THnD>("multi_data_histo_diff_smear_beta_val_easy")); //it's smeared - gen

  //these must match how the 5D histo was produced
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=8, nbinsmll=8, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, mllbinranges, ptbinranges{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  std::cout<<"\n mllbinranges = [";
  for (int i=0; i<=nbinsmll; i++){mllbinranges.push_back(75.0 + i * (105.0 - 75.0)/nbinsmll); std::cout<<mllbinranges[i]<<", ";}
  std::cout<<"] \n";

  /////////////// Prepare variables, histograms //////////////////////////////////////////////////
  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0;
  double max_hist_mll_diff, max_hist_mll, value, error, value_beta, error_beta, diff_squared, jac_b_weight, evts_in_bin, mean_mc, error_mean_mc, sigma_mc=0.0, error_sigma_mc, error_diff_squared, error_jac_b_weight, normalisation_smear, normalisation_smear_beta_val, fit_beta_error, value_corrected, mean_diff_smear_beta_val, error_mean_diff_smear_beta_val, mean_corrected_diff_smear, error_mean_corrected_diff_smear;
  int filled_bins_mll, position_to_fill, filled_bins_mll_diff;
  string name, leg_entry;

  std::map<string, vector<double>> GenRecoFit;
  vector<double> fitresult;

  TH1D *mean_reco = new TH1D("mean", "reco-gen mll mean", 3, 0, 3);
  mean_reco->SetCanExtend(TH1::kAllAxes);
  mean_reco->GetXaxis()->SetTitle("Bin number");
  mean_reco->GetYaxis()->SetTitle("reco-gen mll mean [GeV]");

  TH1D *sigma_reco = new TH1D("sigma_reco", "reco-gen mll sigma", 3, 0, 3);
  sigma_reco->SetCanExtend(TH1::kAllAxes);
  sigma_reco->GetXaxis()->SetTitle("Bin number");
  sigma_reco->GetYaxis()->SetTitle("reco-gen mll sigma [GeV]");

  TH1D *mean_smear = new TH1D("mean_smear", "smear-gen mll mean", 3, 0, 3);
  mean_smear->SetCanExtend(TH1::kAllAxes);
  mean_smear->SetMarkerStyle(kPlus);
  mean_smear->GetXaxis()->SetTitle("Bin number");
  mean_smear->GetYaxis()->SetTitle("smear-gen mll mean [GeV]");

  TH1D *sigma_smear = new TH1D("sigma_smear", "smear-gen mll sigma", 3, 0, 3);
  sigma_smear->SetCanExtend(TH1::kAllAxes);
  sigma_smear->GetXaxis()->SetTitle("Bin number");
  sigma_smear->GetYaxis()->SetTitle("smear-gen mll sigma [GeV]");

  TH1D *gaus_integral = new TH1D("gaus_integral", "reco mll integral 75.5-105.0", 3, 0, 3);
  gaus_integral->SetCanExtend(TH1::kAllAxes);
  gaus_integral->SetMarkerStyle(kPlus);
  gaus_integral->SetMarkerColor(kBlue);
  gaus_integral->SetLineColor(kBlue);
  gaus_integral->GetXaxis()->SetTitle("Bin number");
  gaus_integral->GetYaxis()->SetTitle("integral [-]");

  TH1D *occupancy = new TH1D("bin_occupancy", "bin occupancy", 3, 0, 3);
  occupancy->SetCanExtend(TH1::kAllAxes);
  occupancy->SetMarkerStyle(kPlus);
  occupancy->SetMarkerColor(kBlue);
  occupancy->SetLineColor(kBlue);
  occupancy->GetXaxis()->SetTitle("Bin number");
  occupancy->GetYaxis()->SetTitle("Events");

  TH1D *jac_inclusive = new TH1D("jacobian_inclusive", "jacobian alpha inclusive in eta, pt", nbinsmll, 75.0, 105.0);
  jac_inclusive->GetXaxis()->SetTitle("mll smear");
  jac_inclusive->GetYaxis()->SetTitle("jacobian alpha");

  TH1D *jac_beta_inclusive = new TH1D("jacobian_beta_inclusive", "jacobian beta inclusive in eta, pt", nbinsmll, 75.0, 105.0);
  jac_beta_inclusive->GetXaxis()->SetTitle("mll smear");
  jac_beta_inclusive->GetYaxis()->SetTitle("jacobian beta");

  TH1D *alpha = new TH1D("alpha", "alpha", 3, 0, 3);
  alpha->SetCanExtend(TH1::kAllAxes);
  alpha->SetMarkerStyle(kPlus);
  alpha->SetMarkerColor(kBlue);
  alpha->GetXaxis()->SetTitle("Bin number");
  alpha->GetYaxis()->SetTitle("alpha");

  TH1D *beta = new TH1D("beta", "beta", 3, 0, 3);
  beta->SetCanExtend(TH1::kAllAxes);
  beta->SetMarkerStyle(kPlus);
  beta->SetMarkerColor(kBlue);
  beta->GetXaxis()->SetTitle("Bin number");
  beta->GetYaxis()->SetTitle("beta");

  TH1D *beta_control = new TH1D("beta_control", "mll_diff shifted - mll_diff corrected", 3, 0, 3);
  beta_control->SetCanExtend(TH1::kAllAxes);
  beta_control->SetMarkerStyle(kPlus);
  beta_control->SetMarkerColor(kBlue);
  beta_control->GetXaxis()->SetTitle("Bin number");
  beta_control->GetYaxis()->SetTitle(" ");

  TH1D *pull_beta_control = new TH1D("pull_beta_control", "Pull epsilon mll_diff", 12, -6.0, 6.0);
  pull_beta_control->GetXaxis()->SetTitle("Pull");
  pull_beta_control->GetYaxis()->SetTitle("Events");

  TH1D *nu = new TH1D("nu", "nu", 3, 0, 3);
  nu->SetCanExtend(TH1::kAllAxes);
  nu->SetMarkerStyle(kPlus);
  nu->SetMarkerColor(kBlue);
  nu->GetXaxis()->SetTitle("Bin number");
  nu->GetYaxis()->SetTitle("nu");

  TH1D *nu_control = new TH1D("nu_control", "nu_control", 3, 0, 3);
  nu_control->SetCanExtend(TH1::kAllAxes);
  nu_control->SetMarkerStyle(kPlus);
  nu_control->SetMarkerColor(kBlue);
  nu_control->GetXaxis()->SetTitle("Bin number");
  nu_control->GetYaxis()->SetTitle("nu from mll diff");

  auto multi_hist_proj_diff_reco = mDh_diff_reco->Projection(4);
  auto multi_hist_proj_reco = mDh_reco->Projection(4);
  auto multi_hist_proj_gen = mDh_gen->Projection(4);
  auto multi_hist_proj_smear = mDh_smear->Projection(4);
  auto multi_hist_proj_diff_smear = mDh_diff_smear->Projection(4);
  auto multi_hist_proj_diff_squared_smear = mDh_diff_squared_smear->Projection(4);
  auto multi_hist_proj_diff_squared_smear_control = mDh_diff_squared_smear_control->Projection(4);
  auto multi_hist_proj_jac_beta_smear = mDh_jac_beta_smear->Projection(4);
  auto multi_hist_proj_smear_beta_val = mDh_smear_beta_val->Projection(4);
  auto multi_hist_proj_diff_smear_beta_val = mDh_diff_smear_beta_val->Projection(4);
  auto multi_hist_proj_diff_smear_beta_val_easy = mDh_diff_smear_beta_val_easy->Projection(4);
  auto multi_hist_proj_jac_beta_smear_mll_diff = mDh_jac_beta_smear_mll_diff->Projection(4);

  // Files to write results
  std::unique_ptr<TFile> f1( TFile::Open("control_bin_histo.root", "RECREATE") );
  std::unique_ptr<TFile> f2( TFile::Open("reco_gen_histos.root", "RECREATE") );
  ofstream f3("passed_regions.txt");

  f1->WriteObject(multi_hist_proj_diff_smear, "multi_hist_proj_diff_smear"); 
  f1->WriteObject(multi_hist_proj_diff_squared_smear, "multi_hist_proj_diff_squared_smear");

  // Canvas
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  c1->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  auto leg2 = new TLegend(0.58, 0.68, 0.90, 0.90);

  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.015);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.015);
  leg2->SetFillColor(10);
  leg2->SetNColumns(1);
  leg2->SetHeader("");

  // Prepare counting events
  total_nevents = multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff);
  std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";
  
  all_histos_count = 0.0;
  empty_histos_count = 0.0;
  hfrac = -1.0;  
  efrac = -1.0;
  remaining_nevents = 0.0;

  /////////////// Loop over eta+,pt+,eta-,pt- //////////////////////////////////////////////////
  for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
    mDh_diff_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_gen->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_squared_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_squared_smear_control->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_beta_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_smear_beta_val->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear_beta_val->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear_beta_val_easy->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_beta_smear_mll_diff->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){   
      mDh_diff_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_gen->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_squared_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_squared_smear_control->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_beta_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_smear_beta_val->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear_beta_val->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear_beta_val_easy->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_beta_smear_mll_diff->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	mDh_diff_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_gen->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_squared_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_diff_squared_smear_control->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_beta_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_smear_beta_val->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear_beta_val->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear_beta_val_easy->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_beta_smear_mll_diff->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	  mDh_diff_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_gen->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_squared_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_diff_squared_smear_control->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_jac_beta_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_smear_beta_val->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_smear_beta_val->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_smear_beta_val_easy->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_jac_beta_smear_mll_diff->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);

	  all_histos_count++;
	  
	  delete gROOT->FindObject("multi_data_histo_diff_reco_proj_4");
	  multi_hist_proj_diff_reco = mDh_diff_reco->Projection(4);
	  multi_hist_proj_diff_reco->GetXaxis()->SetTitle("mll_diff [GeV]");
	  multi_hist_proj_diff_reco->GetYaxis()->SetTitle("Events");
	  nevents = multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff); // get stats for sigma_MC fit
	  if (nevents < 100.0){ // reject low stats
	    empty_histos_count++;
	  } else {
	    delete gROOT->FindObject("multi_data_histo_reco_proj_4");
	    multi_hist_proj_reco = mDh_reco->Projection(4);
	    multi_hist_proj_reco->GetXaxis()->SetTitle("mll [GeV]");
	    multi_hist_proj_reco->GetYaxis()->SetTitle("Events");
	    fitresult = fitHisto(multi_hist_proj_reco, 0, 4); // get integral of m_reco Gaussian
	    
	    if (fitresult[4] < 0.75){ // reject small gaus integral 
	      empty_histos_count++;
	    } else {
	      remaining_nevents += nevents;
	      
	      occupancy->SetBinError(occupancy->Fill(name.c_str(), nevents), 100);
	      
	      name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	      f3 << name << "\n";
	      //std::cout<< name <<"\n";

	      /////////////// Reco sigma_MC fit //////////////////////////////////////////////////
	      //multi_hist_proj_diff_reco->Scale( 1.0 / multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff) ); //normalise to unity
	      fitresult = fitHisto(multi_hist_proj_diff_reco, 1, 4);

	      GenRecoFit[name] = fitresult;
	      if (fitresult[0] > -90.0){ mean_reco->SetBinError(mean_reco->Fill(name.c_str(), fitresult[0]), fitresult[2]); }
	      if (fitresult[1] > -5.0){ sigma_reco->SetBinError(sigma_reco->Fill(name.c_str(), fitresult[1]), fitresult[3]); }
	      
	      multi_hist_proj_diff_reco->SetName(name.c_str());
	      multi_hist_proj_diff_reco->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      
	      c1->cd(1); //prepare to draw mll_diff
	      max_hist_mll_diff = -1.0;
	      multi_hist_proj_diff_reco->SetLineColor(kBlue);
	      
	      max_hist_mll_diff = multi_hist_proj_diff_reco->GetBinContent(multi_hist_proj_diff_reco->GetMaximumBin());
	      
	      /////////////// Smear sigma_MC fit //////////////////////////////////////////////////
	      delete gROOT->FindObject("multi_data_histo_diff_smear_proj_4");
	      multi_hist_proj_diff_smear = mDh_diff_smear->Projection(4);
	      multi_hist_proj_diff_smear->GetXaxis()->SetTitle("mll_diff [GeV]");
	      multi_hist_proj_diff_smear->GetYaxis()->SetTitle("Events");
	      //multi_hist_proj_diff_smear->Scale( 1.0 / multi_hist_proj_diff_smear->Integral(1,nbinsmll_diff) ); //normalise to unity
	      fitresult = fitHisto(multi_hist_proj_diff_smear, 1, 8);
	      if (fitresult[0] > -90.0){ mean_smear->SetBinError(mean_smear->Fill(name.c_str(), fitresult[0]), fitresult[2]); }
	      //save for alpha, beta fit
	      mean_mc = fitresult[0]; 
	      error_mean_mc = fitresult[2];
	      sigma_mc = fitresult[1]; 
	      error_sigma_mc = fitresult[3]; 

	      if (fitresult[1] > -5.0){ sigma_smear->SetBinError(sigma_smear->Fill(name.c_str(), fitresult[1]), fitresult[3]); }
	      multi_hist_proj_diff_smear->SetLineColor(kGreen);

	      delete gROOT->FindObject("multi_data_histo_diff_smear_beta_val_proj_4");
	      multi_hist_proj_diff_smear_beta_val = mDh_diff_smear_beta_val->Projection(4);
	      multi_hist_proj_diff_smear_beta_val->GetXaxis()->SetTitle("mll_diff [GeV]");
              multi_hist_proj_diff_smear_beta_val->GetYaxis()->SetTitle("Events");
	      //multi_hist_proj_diff_smear_beta_val->Scale( 1.0 / multi_hist_proj_diff_smear_beta_val->Integral(1,nbinsmll_diff) ); //normalise to unity
	      fitresult = fitHisto(multi_hist_proj_diff_smear_beta_val, 1, 5);
	      // save these below for beta_control histogram 
	      mean_diff_smear_beta_val = fitresult[0]; 
	      error_mean_diff_smear_beta_val = fitresult[2];

	      multi_hist_proj_diff_smear_beta_val->SetLineColor(kYellow);
	      if(max_hist_mll_diff < multi_hist_proj_diff_smear_beta_val->GetBinContent(multi_hist_proj_diff_smear_beta_val->GetMaximumBin())){
                max_hist_mll_diff = multi_hist_proj_diff_smear_beta_val->GetBinContent(multi_hist_proj_diff_smear_beta_val->GetMaximumBin());
              }

	      delete gROOT->FindObject("multi_data_histo_diff_smear_beta_val_easy_proj_4");
	      multi_hist_proj_diff_smear_beta_val_easy = mDh_diff_smear_beta_val_easy->Projection(4);
	      multi_hist_proj_diff_smear_beta_val_easy->GetXaxis()->SetTitle("mll_diff [GeV]");
              multi_hist_proj_diff_smear_beta_val_easy->GetYaxis()->SetTitle("Events");
	      fitresult = fitHisto(multi_hist_proj_diff_smear_beta_val_easy, 1, 51);
	      multi_hist_proj_diff_smear_beta_val_easy->SetLineColor(kViolet);
	      if(max_hist_mll_diff < multi_hist_proj_diff_smear_beta_val_easy->GetBinContent(multi_hist_proj_diff_smear_beta_val_easy->GetMaximumBin())){
		max_hist_mll_diff = multi_hist_proj_diff_smear_beta_val_easy->GetBinContent(multi_hist_proj_diff_smear_beta_val_easy->GetMaximumBin());
              }
	      
	      /////////////// Draw mll_diff ///////////////////////////////////////////////////////
	      if(max_hist_mll_diff < multi_hist_proj_diff_smear->GetBinContent(multi_hist_proj_diff_smear->GetMaximumBin())){
		max_hist_mll_diff = multi_hist_proj_diff_smear->GetBinContent(multi_hist_proj_diff_smear->GetMaximumBin());
	      } 
	      
	      max_hist_mll_diff = max_hist_mll_diff * 1.3;
              multi_hist_proj_diff_smear->SetMaximum(max_hist_mll_diff);
	      multi_hist_proj_diff_smear->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      multi_hist_proj_diff_smear->Draw();
	      multi_hist_proj_diff_smear_beta_val->Draw("SAME");
	      //multi_hist_proj_diff_smear_beta_val_easy->Draw("SAME");
	      //multi_hist_proj_diff_reco->Draw("SAME");

	      leg1->Clear();
	      leg_entry = "Region " + name + ": " + nevents + " events";
	      leg1->SetHeader(leg_entry.c_str(),"C");
              
	      //leg1->AddEntry(multi_hist_proj_diff_reco, "reco-gen", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear, "smeared-gen, #varepsilon=0", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear_beta_val, "smeared-gen, #varepsilon=-0.15,by weight", "l");
	      //leg1->AddEntry(multi_hist_proj_diff_smear_beta_val_easy, "smeared-gen, #varepsilon=-0.15,by sampling", "l");
	      leg_entry = "#varepsilon=0: #mu=" + to_string(mean_mc).substr(0, 5) + ", #sigma=" + to_string(sigma_mc).substr(0, 5);
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      //leg1->Draw("");
	      
	      c1->cd(2); //prepare to draw mll
	      max_hist_mll = -1.0;
	      	      
	      multi_hist_proj_reco->SetLineColor(kBlue);
	      //multi_hist_proj_reco->Scale(1.0 / multi_hist_proj_reco->Integral(1,nbinsmll)); //normalise to unity
	      fitresult = fitHisto(multi_hist_proj_reco, 0, 4);
	      multi_hist_proj_reco->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      max_hist_mll = multi_hist_proj_reco->GetBinContent(multi_hist_proj_reco->GetMaximumBin());
	      if (fitresult[4] > -100.0){ gaus_integral->SetBinError(gaus_integral->Fill(name.c_str(), fitresult[4]), 0.01); }
	      delete gROOT->FindObject("multi_data_histo_gen_proj_4");
	      multi_hist_proj_gen = mDh_gen->Projection(4);
	      //multi_hist_proj_gen->Scale(1.0 / multi_hist_proj_gen->Integral(1,nbinsmll)); //normalise to unity
	      multi_hist_proj_gen->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      multi_hist_proj_gen->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_gen->GetYaxis()->SetTitle("Events");
	      //fitresult = fitHisto(multi_hist_proj_gen, 0, 2);
	      multi_hist_proj_gen->SetLineColor(kRed);
	      if(max_hist_mll < multi_hist_proj_gen->GetBinContent(multi_hist_proj_gen->GetMaximumBin())){
		max_hist_mll = multi_hist_proj_gen->GetBinContent(multi_hist_proj_gen->GetMaximumBin());
	      }
	      delete gROOT->FindObject("multi_data_histo_smear_proj_4");
	      multi_hist_proj_smear = mDh_smear->Projection(4);
	      normalisation_smear = 1.0 / multi_hist_proj_smear->Integral(1,nbinsmll); 
	      //multi_hist_proj_smear->Scale(normalisation_smear); //normalise to unity
	      multi_hist_proj_smear->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_smear->GetYaxis()->SetTitle("Events");
	      if(max_hist_mll < multi_hist_proj_smear->GetBinContent(multi_hist_proj_smear->GetMaximumBin())){
                max_hist_mll = multi_hist_proj_smear->GetBinContent(multi_hist_proj_smear->GetMaximumBin());
	      }

	      delete gROOT->FindObject("multi_data_histo_smear_beta_val_proj_4");
              multi_hist_proj_smear_beta_val = mDh_smear_beta_val->Projection(4);
	      normalisation_smear_beta_val = 1.0 / multi_hist_proj_smear_beta_val->Integral(1,nbinsmll);
	      //multi_hist_proj_smear_beta_val->Scale(normalisation_smear_beta_val); //normalise to unity
              multi_hist_proj_smear_beta_val->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_smear_beta_val->GetYaxis()->SetTitle("Events");
              if(max_hist_mll < multi_hist_proj_smear_beta_val->GetBinContent(multi_hist_proj_smear_beta_val->GetMaximumBin())){
                max_hist_mll = multi_hist_proj_smear_beta_val->GetBinContent(multi_hist_proj_smear_beta_val->GetMaximumBin());
	      }

	      /////////////// Start drawing mll ///////////////////////////////////////////////////////
	      //fitresult = fitHisto(multi_hist_proj_smear, 0, 8);
	      multi_hist_proj_smear->SetLineColor(kGreen);
	      //fitresult = fitHisto(multi_hist_proj_smear_beta_val, 0, 5);
	      multi_hist_proj_smear_beta_val->SetLineColor(kYellow);
	      
	      max_hist_mll = max_hist_mll * 1.3;
	      multi_hist_proj_gen->SetMaximum(max_hist_mll);
	      //multi_hist_proj_reco->Draw();
	      multi_hist_proj_gen->Draw();
	      multi_hist_proj_smear->Draw("SAME");
	      multi_hist_proj_smear_beta_val->Draw("SAME");
	      	      
	      leg2->Clear();
	      leg_entry = "Region " + name + ": " + nevents + " events";
	      leg2->SetHeader(leg_entry.c_str(),"C"); 

	      //leg2->AddEntry(multi_hist_proj_reco, "reco", "l");
	      leg2->AddEntry(multi_hist_proj_gen, "gen", "l");
	      leg2->AddEntry(multi_hist_proj_smear, "smeared, #beta=1", "l");
	      leg2->AddEntry(multi_hist_proj_smear_beta_val, "smeared, #beta=0.995", "l"); //value of beta goes here !!!!!!
	      	      
	      /////////////// Fill vectors and variance for minimisation //////////////////////////////////////////
	      filled_bins_mll=0;
	      vector<int> good_indices_mll;
	      for(int i=1; i<=nbinsmll; i++){
		if (multi_hist_proj_smear_beta_val->GetBinContent(i) > 10 && multi_hist_proj_smear->GetBinContent(i) > 10 && multi_hist_proj_gen->GetBinContent(i) > 0){
		  good_indices_mll.push_back(i);
		  filled_bins_mll++;
		}
	      }
	      if (filled_bins_mll <= 1 ){continue;} //This needs refined

	      filled_bins_mll_diff=0;
	      vector<int> good_indices_mll_diff;
	      for(int i=1; i<=nbinsmll_diff; i++){
		if( multi_hist_proj_diff_smear->GetBinContent(i) >= 10 && multi_hist_proj_diff_smear_beta_val->GetBinContent(i) > 0 ){ //asking jacobian to be average over at least 10 events
		  good_indices_mll_diff.push_back(i);
		  filled_bins_mll_diff++;
		}
	      }
	      if (filled_bins_mll_diff <= 1 ){continue;} //This needs refined
	      
	      VectorXd h_smear_minus_smear_vector(filled_bins_mll);
	      Eigen::MatrixXd V_inv_sqrt(filled_bins_mll, filled_bins_mll), J(filled_bins_mll, 2), J_beta(filled_bins_mll, 2);

              V_inv_sqrt=MatrixXd::Zero(filled_bins_mll, filled_bins_mll);
	      position_to_fill=0;
              for(int i=1; i<=nbinsmll; i++){
		if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
		  h_smear_minus_smear_vector(position_to_fill) = multi_hist_proj_smear_beta_val->GetBinContent(i) - multi_hist_proj_smear->GetBinContent(i);
		  J(position_to_fill,0) = multi_hist_proj_smear->GetBinContent(i);
		  J_beta(position_to_fill,0) = multi_hist_proj_smear->GetBinContent(i);
		  V_inv_sqrt(position_to_fill,position_to_fill) = 1 / multi_hist_proj_smear_beta_val->GetBinErrorLow(i);
		  position_to_fill++;
		}
	      }
	      if (position_to_fill != filled_bins_mll){ std::cout<<"problem counting vector size \n"; }
              
              VectorXd h_smear_diff_minus_smear_diff_vector(filled_bins_mll_diff);
	      Eigen::MatrixXd V_inv_sqrt_control(filled_bins_mll_diff, filled_bins_mll_diff), J_beta_control(filled_bins_mll_diff, 2);

	      V_inv_sqrt_control=MatrixXd::Zero(filled_bins_mll_diff, filled_bins_mll_diff);
	      position_to_fill=0;
	      for(int i=1; i<=nbinsmll_diff; i++){
		if( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		  h_smear_diff_minus_smear_diff_vector(position_to_fill) = multi_hist_proj_diff_smear_beta_val->GetBinContent(i) - multi_hist_proj_diff_smear->GetBinContent(i);
		  J_beta_control(position_to_fill,0) = multi_hist_proj_diff_smear->GetBinContent(i);
		  V_inv_sqrt_control(position_to_fill,position_to_fill) = 1 / multi_hist_proj_diff_smear_beta_val->GetBinErrorLow(i);
		  position_to_fill++;
		}
	      }
	      if (position_to_fill != filled_bins_mll_diff){ std::cout<<"problem counting vector size \n"; }

	      /////////////// Jacobian ///////////////////////////////////////////////////////
	      delete gROOT->FindObject("jacobian");
	      delete gROOT->FindObject("jacobian_control");
	      TH1D *jac = new TH1D("jacobian", "jacobian", nbinsmll, 75.0, 105.0);
	      jac->GetXaxis()->SetTitle("mll_smear [GeV]");
	      jac->GetYaxis()->SetTitle("GeV");
              TH1D *jac_control = new TH1D("jacobian_control", "jacobian", nbinsmll_diff, -5.0, 5.0);
              jac_control->GetXaxis()->SetTitle("mll_smear(#beta=1) - mll_gen [GeV]");
              jac_control->GetYaxis()->SetTitle("GeV");

	      delete gROOT->FindObject("jacobian_beta");
	      delete gROOT->FindObject("jacobian_beta_control");
	      TH1D *jac_beta = new TH1D("jacobian_beta", "jacobian beta", nbinsmll, 75.0, 105.0);
              jac_beta->GetXaxis()->SetTitle("mll_smear [GeV]");
              jac_beta->GetYaxis()->SetTitle("GeV");
	      TH1D *jac_beta_control = new TH1D("jacobian_beta_control", "jacobian beta", nbinsmll_diff, -5.0, 5.0);
              jac_beta_control->GetXaxis()->SetTitle("mll_smear(#beta=1) - mll_gen [GeV]");
              jac_beta_control->GetYaxis()->SetTitle("GeV");
	      
	      //Fill jacobian bin by bin
	      delete gROOT->FindObject("multi_data_histo_diff_squared_smear_proj_4");
	      multi_hist_proj_diff_squared_smear = mDh_diff_squared_smear->Projection(4);
	      delete gROOT->FindObject("multi_data_histo_jac_beta_smear_proj_4");
              multi_hist_proj_jac_beta_smear = mDh_jac_beta_smear->Projection(4);

	      position_to_fill=0;
	      for(int i=1; i<=nbinsmll; i++){
		diff_squared = multi_hist_proj_diff_squared_smear->GetBinContent(i);
		error_diff_squared = multi_hist_proj_diff_squared_smear->GetBinErrorLow(i);

		jac_b_weight = multi_hist_proj_jac_beta_smear->GetBinContent(i);
                error_jac_b_weight = multi_hist_proj_jac_beta_smear->GetBinErrorLow(i);
	
		if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
		  evts_in_bin = multi_hist_proj_smear->GetBinContent(i);
		  value =  diff_squared /  (evts_in_bin * sigma_mc * sigma_mc) - 1;
		  error = 1 / (evts_in_bin * sigma_mc*sigma_mc) * pow((4 * diff_squared*diff_squared * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_diff_squared*error_diff_squared) , 0.5);
		  J(position_to_fill,1) = value;

		  jac->SetBinContent(i, value);
		  jac->SetBinError(i, error);
		  jac_inclusive->Fill(mllbinranges[i-1], value);

		  value_beta = jac_b_weight / (evts_in_bin * sigma_mc * sigma_mc);
		  error_beta = 1 / (evts_in_bin * sigma_mc*sigma_mc) * pow((4 * jac_b_weight*jac_b_weight * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_jac_b_weight*error_jac_b_weight) , 0.5);
		  J_beta(position_to_fill,1) = value_beta;

		  jac_beta->SetBinContent(i, value_beta);
		  jac_beta->SetBinError(i, error_beta);
		  jac_beta_inclusive->Fill(mllbinranges[i-1], value_beta);
		  
		  position_to_fill++;
		} 
	      }
	      //jac_beta->Scale(1.0 / jac_beta->Integral(1,nbinsmll)); //normalise jacobian to unity
	      //position_to_fill = 0;
	      //for(int i=1; i<=nbinsmll; i++){
	      //if (jac_beta->GetBinContent(i) != 0 ){
	      //  J_beta(position_to_fill) = jac_beta->GetBinContent(i);
	      //  position_to_fill++;
	      //}
	      //}
	      if (position_to_fill != filled_bins_mll){ std::cout<<"problem counting jac size \n"; }

	      ////////// Jacobian alpha vs m_diff /////////////////////////////////////////
              delete gROOT->FindObject("multi_data_histo_diff_squared_smear_control_proj_4");
              multi_hist_proj_diff_squared_smear_control = mDh_diff_squared_smear_control->Projection(4);
	      for(int i=1; i<=nbinsmll_diff; i++){
		diff_squared = multi_hist_proj_diff_squared_smear_control->GetBinContent(i);
		error_diff_squared = multi_hist_proj_diff_squared_smear_control->GetBinErrorLow(i);
		evts_in_bin = multi_hist_proj_diff_smear->GetBinContent(i);
		if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		  value =  diff_squared /  (evts_in_bin * sigma_mc * sigma_mc) - 1;
		  error = 1 / (evts_in_bin * sigma_mc*sigma_mc) * pow((4 * diff_squared*diff_squared * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_diff_squared*error_diff_squared) , 0.5);
		  
		  jac_control->SetBinContent(i, value);
		  jac_control->SetBinError(i, error);
		
		}
	      }

	      ////////// Jacobian beta vs m_diff /////////////////////////////////////////
              delete gROOT->FindObject("multi_data_histo_jac_beta_smear_mll_diff_proj_4");
              multi_hist_proj_jac_beta_smear_mll_diff = mDh_jac_beta_smear_mll_diff->Projection(4);
	      position_to_fill = 0;
              for(int i=1; i<=nbinsmll_diff; i++){
                jac_b_weight = multi_hist_proj_jac_beta_smear_mll_diff->GetBinContent(i);
                error_jac_b_weight = multi_hist_proj_jac_beta_smear_mll_diff->GetBinErrorLow(i);
                evts_in_bin = multi_hist_proj_diff_smear->GetBinContent(i);
                
                if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
                  //value_beta = ( jac_b_weight/evts_in_bin - mean_mc ) / pow(sigma_mc,2);  
		  value_beta = ( jac_b_weight - evts_in_bin*mean_mc ) / pow(sigma_mc,2);
		  //error_beta = pow( pow(error_mean_mc,2) + 4*pow(jac_b_weight/evts_in_bin - mean_mc,2)*pow(error_sigma_mc,2)/pow(sigma_mc,2) + pow(error_jac_b_weight,2)/pow(evts_in_bin,2), 0.5) / pow(sigma_mc,2);
		  error_beta = pow( pow(evts_in_bin*error_mean_mc,2) + 4*pow(jac_b_weight - evts_in_bin*mean_mc,2)*pow(error_sigma_mc,2)/pow(sigma_mc,2) + pow(error_jac_b_weight,2), 0.5) / pow(sigma_mc,2);
		  J_beta_control(position_to_fill,1) = value_beta;
		  
		  jac_beta_control->SetBinContent(i, value_beta);
		  jac_beta_control->SetBinError(i, error_beta);

                  position_to_fill++;

                } 
	      }

              //jac_beta_control->Scale(1.0 / jac_beta_control->Integral(1,nbinsmll_diff)); //normalise jacobian to no. of events in MC
              //position_to_fill = 0;
              //for(int i=1; i<=nbinsmll_diff; i++){
	      //if ( multi_hist_proj_diff_smear_beta_val->GetBinContent(i) >= 0 && multi_hist_proj_diff_smear->GetBinContent(i) > 0 ){ //////is this the corrrect way of counting this?????
	      //  J_beta_control(position_to_fill) = jac_beta_control->GetBinContent(i);
	      //  position_to_fill++;
	      //}
              //}
              if (position_to_fill != filled_bins_mll_diff){ std::cout<<"problem counting jac_beta_control size \n"; }

	      //write jacobians
	      f2->WriteObject(jac, ("jac" + name).c_str());
	      f2->WriteObject(jac_control, ("jac_control" + name).c_str());
	      f2->WriteObject(jac_beta, ("jac_beta" + name).c_str());
	      f2->WriteObject(jac_beta_control, ("jac_beta_control" + name).c_str());
		      
	      ////////////////// solve for alpha //////////////////
	      Eigen::MatrixXd A = V_inv_sqrt*J;
	      Eigen::MatrixXd b = V_inv_sqrt*h_smear_minus_smear_vector;
 //!!!!!!!!!!! Attention alpha_vector contains alpha-1, not alpha !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      Eigen::VectorXd alpha_vector = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	      alpha->SetBinError(alpha->Fill(name.c_str(), alpha_vector(1)+1), 0.01); //write alpha
	      
	      ////////////////// solve for beta //////////////////
	      Eigen::MatrixXd A_beta = V_inv_sqrt*J_beta;
 //!!!!!!!!!!! Attention beta_vector contains beta-1, not beta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      Eigen::VectorXd beta_vector = A_beta.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

	      ////////////////// error on beta //////////////////
	      Eigen::MatrixXd V_beta = (J_beta.col(1).transpose()*(V_inv_sqrt*V_inv_sqrt)*J_beta.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      Eigen::MatrixXd V_nu = (J_beta.col(0).transpose()*(V_inv_sqrt*V_inv_sqrt)*J_beta.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      fit_beta_error = pow(V_beta(0,0),0.5);
	      beta->SetBinError(beta->Fill(name.c_str(), beta_vector(1)+1), fit_beta_error); //write beta
	      nu->SetBinError(nu->Fill(name.c_str(), beta_vector(0)+1), pow(V_nu(0,0),0.5));

	      ////////////////// solve for beta //////////////////
	      Eigen::MatrixXd A_beta_control = V_inv_sqrt_control*J_beta_control;
	      Eigen::MatrixXd b_control = V_inv_sqrt_control*h_smear_diff_minus_smear_diff_vector;
	      // for beta_vector_control, we actually need just beta, instead of beta-1 
	      Eigen::VectorXd beta_vector_control = A_beta_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_control);
              ////////////////// error on beta //////////////////
	      Eigen::MatrixXd V_beta_control = (J_beta_control.col(1).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_beta_control.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      Eigen::MatrixXd V_nu_control = (J_beta_control.col(0).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_beta_control.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      nu_control->SetBinError(nu_control->Fill(name.c_str(), beta_vector_control(0)+1), pow(V_nu_control(0,0),0.5)); 

	      //////////////// Finish drawing mll ////////////////////////////////
	      leg_entry = "Fit #beta=" + to_string(beta_vector(1)+1).substr(0, 6) + "#pm" + to_string(fit_beta_error).substr(0, 6); 
	      leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #nu=" + to_string(beta_vector(0)+1).substr(0, 6) + "#pm" + to_string(pow(V_nu(0,0),0.5)).substr(0, 6);
	      leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      
	      leg2->Draw("");

              c1->cd();
	      c1->cd(1);
	      leg_entry = "Fit #varepsilon=" + to_string(beta_vector_control(1)).substr(0, 6) + "#pm" + to_string(pow(V_beta_control(0,0),0.5)).substr(0, 6); //write beta, no need to add one
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #nu=" + to_string(beta_vector_control(0)+1).substr(0, 6) + "#pm" + to_string(pow(V_nu_control(0,0),0.5)).substr(0, 6); 
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      
	      // apply fitted epsilon and nu correction
	      delete gROOT->FindObject("corrected_diff_smear");
	      TH1D *corrected_diff_smear = new TH1D("corrected_diff_smear", "corrected_diff_smear", nbinsmll_diff, -5.0, 5.0);
	      corrected_diff_smear->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      corrected_diff_smear->GetXaxis()->SetTitle("mll_diff [GeV]");
              corrected_diff_smear->GetYaxis()->SetTitle("Events");
	      
	      position_to_fill = 0;
              for(int i=1; i<=nbinsmll_diff; i++){
		if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		  value_corrected = multi_hist_proj_diff_smear->GetBinContent(i) + beta_vector_control(1)*J_beta_control(position_to_fill,1) + beta_vector_control(0)*J_beta_control(position_to_fill,0);
		  error = pow( pow(J_beta_control(position_to_fill,1),2)*V_beta_control(0,0) + pow(beta_vector_control(1)*jac_beta_control->GetBinErrorLow(i),2) + pow(multi_hist_proj_diff_smear->GetBinContent(i),2)*V_nu_control(0,0) , 0.5);
		  corrected_diff_smear->SetBinContent(i, value_corrected);
		  corrected_diff_smear->SetBinError(i, error);

		  position_to_fill++;
		}
	      }
	      if (position_to_fill != filled_bins_mll_diff){ std::cout<<"problem counting vector size \n"; }
	      
	      //draw it, add to leg
	      fitresult = fitHisto(corrected_diff_smear, 1, 1);
	      // save for beta_control
	      mean_corrected_diff_smear = fitresult[0];
	      error_mean_corrected_diff_smear = fitresult[2];
              corrected_diff_smear->SetLineColor(kBlack);
	      corrected_diff_smear->Draw("SAME");

	      leg1->AddEntry(corrected_diff_smear, "corrected smeared-gen, #varepsilon from 0 to -0.15", "l");
	      leg1->Draw("");
	      
	      c1->cd();
              f2->WriteObject(c1, name.c_str());

	      beta_control->SetBinError(beta_control->Fill(name.c_str(), mean_diff_smear_beta_val - mean_corrected_diff_smear), pow( pow(error_mean_diff_smear_beta_val,2) + pow(error_mean_corrected_diff_smear,2) ,0.5) );
	      pull_beta_control->Fill( (mean_diff_smear_beta_val - mean_mc - beta_vector_control(1)) / pow(V_beta_control(0,0),0.5) );

	    }
	  }
	  
	} 
      }
    }
  }
  f3.close();

  hfrac = empty_histos_count / all_histos_count;
  efrac = remaining_nevents / total_nevents;
  std::cout<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";

  ////////////////// Write remaining histograms ///////////////////////////

  mean_reco->SetStats(0);
  mean_reco->LabelsDeflate();
  f1->WriteObject(mean_reco, "mean_diff_reco");
  
  sigma_reco->SetStats(0);
  sigma_reco->LabelsDeflate();
  f1->WriteObject(sigma_reco, "sigma_diff_reco");

  mean_smear->SetStats(0);
  mean_smear->LabelsDeflate();
  f1->WriteObject(mean_smear, "mean_diff_smear");

  sigma_smear->SetStats(0);
  sigma_smear->LabelsDeflate();
  f1->WriteObject(sigma_smear, "sigma_diff_smear");

  f1->WriteObject(alpha, "alpha");
  f1->WriteObject(jac_inclusive, "jacobian_inclusive");

  beta->Fit("pol0");
  f1->WriteObject(beta, "beta");

  std::cout<<"beta control: \n";
  beta_control->Fit("pol0");
  f1->WriteObject(beta_control, "beta_control");

  fitresult = fitHisto(pull_beta_control, 1, 1);
  std::cout<<"Pull epsilon control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
  f1->WriteObject(pull_beta_control, "pull_beta_control");

  f1->WriteObject(nu, "nu");
  f1->WriteObject(nu_control, "nu_control");

  f1->WriteObject(jac_beta_inclusive, "jacobian_beta_inclusive");

  // Superimposed histograms

  TCanvas *c2 = new TCanvas("c2","c2",800,600);

  auto leg3 = new TLegend(0.58, 0.68, 0.90, 0.90);

  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.035);
  leg3->SetFillColor(10);
  leg3->SetNColumns(1);
  leg3->SetHeader("");

  mean_reco->Draw();
  mean_smear->SetMarkerColor(kRed);
  mean_smear->SetLineColor(kRed);
  mean_smear->Draw("SAME");

  leg3->AddEntry(mean_reco, "reco", "l");
  leg3->AddEntry(mean_smear, "smear", "l");
  leg3->Draw("");

  f1->WriteObject(c2, "mean_superimposed");

  gaus_integral->SetStats(0);
  gaus_integral->LabelsDeflate();
  f1->WriteObject(gaus_integral, "gaus_integral");

  occupancy->SetStats(0);
  occupancy->LabelsDeflate();
  f1->WriteObject(occupancy, "bin_occupancy");

  return 0; 
}
