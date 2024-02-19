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

//--------------------------------------------------------------------------------------
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
vector<double> fitHisto(TH1* histogram, int draw_option, int color){

  vector<double> fitresult;

  double mean = histogram->GetMean();
  double sigma = histogram->GetRMS();
  double mean_err, sigma_err, integral;

  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", mean - 3 * sigma, mean + 3 * sigma);
  // first fit
  if(0 == histogram->Fit(gaussianFunc, "QNR")){
    mean = gaussianFunc->GetParameter(1);
    sigma = gaussianFunc->GetParameter(2);

    // second fit: few sigma of first fit around mean of first fit
    gaussianFunc->SetRange(mean - 5 * sigma, mean + 5 * sigma);
    if (0 == histogram->Fit(gaussianFunc, "Q0R")) { // don't draw yet
      if (histogram->GetFunction(gaussianFunc->GetName())) { 
	histogram->GetFunction(gaussianFunc->GetName())->SetLineColor(color);
        if (draw_option == 1){ histogram->GetFunction(gaussianFunc->GetName())->ResetBit(TF1::kNotDraw);} // draw fit
      }
      mean = gaussianFunc->GetParameter(1);
      sigma = gaussianFunc->GetParameter(2);
      mean_err = gaussianFunc->GetParError(1);
      sigma_err = gaussianFunc->GetParError(2);

      TF1 *gaussianNewFunc = new TF1("gaussianNewFunc", "gaus", mean - 5 * sigma, mean + 5 * sigma);
      gaussianNewFunc->SetParameter(0, 1.0/(sigma*pow(2*M_PI,0.5)));
      gaussianNewFunc->SetParameter(1, mean);
      gaussianNewFunc->SetParameter(2, sigma);

      // integral of fitted gaussian function in range of the histogram
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

//--------------------------------------------------------------------------------------

int frame(){

  ROOT::EnableImplicitMT(128);

  //--------------------------------------------------------------------------------------

  // Input files 

  // Reco
  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo_reco.root") );
  std::unique_ptr<THnD> mDh_reco(myFile->Get<THnD>("multi_data_histo_reco"));
  std::unique_ptr<THnD> mDh_gen(myFile->Get<THnD>("multi_data_histo_gen"));
  std::unique_ptr<THnD> mDh_diff_reco(myFile->Get<THnD>("multi_data_histo_diff_reco")); // it's reco - gen

  // Jacobian terms reco
  std::unique_ptr<THnD> mDh_jac_diff_squared_reco_mll(myFile->Get<THnD>("multi_data_histo_jac_diff_squared_reco_mll")); // sum of mll_diff_squared in mll bin reco
  std::unique_ptr<THnD> mDh_jac_diff_squared_reco_mll_diff(myFile->Get<THnD>("multi_data_histo_jac_diff_squared_reco_mll_diff")); // sum of mll_diff_squared in mll_diff bin reco
  std::unique_ptr<THnD> mDh_jac_diff_reco_mll(myFile->Get<THnD>("multi_data_histo_jac_diff_reco_mll")); // sum of mll_diff in mll bin reco
  std::unique_ptr<THnD> mDh_jac_diff_reco_mll_diff(myFile->Get<THnD>("multi_data_histo_jac_diff_reco_mll_diff")); // sum of mll_diff in mll_diff bin reco

  // Data

  //--------------------------------------------------------------------------------------

  // Smear
  std::unique_ptr<TFile> myFile2( TFile::Open("multiD_histo_smear.root") );
  std::unique_ptr<THnD> mDh_smear(myFile2->Get<THnD>("multi_data_histo_smear"));
  std::unique_ptr<THnD> mDh_diff_smear(myFile2->Get<THnD>("multi_data_histo_diff_smear")); //it's smeared - gen
  // Shift weights
  std::unique_ptr<THnD> mDh_smear_beta_val(myFile2->Get<THnD>("multi_data_histo_smear_beta_val"));
  std::unique_ptr<THnD> mDh_diff_smear_beta_val(myFile2->Get<THnD>("multi_data_histo_diff_smear_beta_val")); 
    
  // Jacobian terms smear
  // analytical jacobians
  std::unique_ptr<THnD> mDh_jac_diff_squared_smear_mll(myFile2->Get<THnD>("multi_data_histo_jac_diff_squared_smear_mll")); // sum of mll_diff_squared in mll bin smear
  std::unique_ptr<THnD> mDh_jac_diff_squared_smear_mll_diff(myFile2->Get<THnD>("multi_data_histo_jac_diff_squared_smear_mll_diff")); // sum of mll_diff_squared in mll_diff bin smear
  std::unique_ptr<THnD> mDh_jac_diff_smear_mll(myFile2->Get<THnD>("multi_data_histo_jac_diff_smear_mll")); // sum of mll_diff in mll bin smear
  std::unique_ptr<THnD> mDh_jac_diff_smear_mll_diff(myFile2->Get<THnD>("multi_data_histo_jac_diff_smear_mll_diff")); // sum of mll_diff in mll_diff bin smear
  // numerical jacobians
  std::unique_ptr<THnD> mDh_diff_smear_plus_offset(myFile2->Get<THnD>("multi_data_histo_diff_smear_plus_offset")); 
  std::unique_ptr<THnD> mDh_diff_smear_minus_offset(myFile2->Get<THnD>("multi_data_histo_diff_smear_minus_offset"));
  float mll_offset = -0.1;

  // Smear easy way
  std::unique_ptr<TFile> myFile3( TFile::Open("multiD_histo_smear_beta_val_easy.root") );
  std::unique_ptr<THnD> mDh_diff_smear_beta_val_easy(myFile3->Get<THnD>("multi_data_histo_diff_smear_beta_val_easy")); //it's smeared_beta_val_easy - gen  

  //--------------------------------------------------------------------------------------

  // Binning must match with 5D histogram
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=8, nbinsmll=32, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, mllbinranges, ptbinranges{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};
  
  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  std::cout<<"\n mllbinranges = [";
  for (int i=0; i<=nbinsmll; i++){mllbinranges.push_back(75.0 + i * (105.0 - 75.0)/nbinsmll); std::cout<<mllbinranges[i]<<", ";}
  std::cout<<"] \n";

  //--------------------------------------------------------------------------------------

  // Variables declaration 

  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0;
  double max_hist_mll_diff, max_hist_mll, value, error, value_epsilon, error_epsilon, diff_squared, diff, jac_e_weight, jac_b_weight, evts_in_bin, mean_mc, error_mean_mc, sigma_mc=0.0, error_sigma_mc, error_diff_squared, error_diff, error_jac_e_weight, error_jac_b_weight, fit_epsilon_error, value_corrected, mean_diff_smear_epsilon_val, error_mean_diff_smear_epsilon_val, sigma_diff_smear_epsilon_val, error_sigma_diff_smear_epsilon_val, mean_corrected_diff_smear, error_mean_corrected_diff_smear, sigma_corrected_diff_smear, error_sigma_corrected_diff_smear, integral_diff_smear_epsilon_val, integral_mc;
  int filled_bins_mll, position_to_fill, filled_bins_mll_diff, unchecked_fits=0;
  string name, leg_entry;
  vector<double> fitresult;

  // Histograms for mll_diff distribution properties

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

  TH1D *occupancy = new TH1D("bin_occupancy", "Events in mll_diff_reco", 3, 0, 3);
  occupancy->SetCanExtend(TH1::kAllAxes);
  occupancy->SetMarkerStyle(kPlus);
  occupancy->SetMarkerColor(kBlue);
  occupancy->SetLineColor(kBlue);
  occupancy->GetXaxis()->SetTitle("Bin number");
  occupancy->GetYaxis()->SetTitle("Events");

  // Histograms for mll distribution properties

  TH1D *gaus_integral = new TH1D("gaus_integral", "reco mll integral 75.5-105.0", 3, 0, 3);
  gaus_integral->SetCanExtend(TH1::kAllAxes);
  gaus_integral->SetMarkerStyle(kPlus);
  gaus_integral->SetMarkerColor(kBlue);
  gaus_integral->SetLineColor(kBlue);
  gaus_integral->GetXaxis()->SetTitle("Bin number");
  gaus_integral->GetYaxis()->SetTitle("integral [-]");

  // Histograms for jacobians inclusive in eta, pt

  TH1D *jac_inclusive = new TH1D("jacobian_inclusive", "jacobian alpha inclusive in eta, pt", nbinsmll, 75.0, 105.0);
  jac_inclusive->GetXaxis()->SetTitle("mll smear");
  jac_inclusive->GetYaxis()->SetTitle("jacobian alpha");

  TH1D *jac_epsilon_inclusive = new TH1D("jacobian_epsilon_inclusive", "jacobian epsilon inclusive in eta, pt", nbinsmll, 75.0, 105.0);
  jac_epsilon_inclusive->GetXaxis()->SetTitle("mll smear");
  jac_epsilon_inclusive->GetYaxis()->SetTitle("jacobian epsilon");

  // Histograms for pull of fitted variables per k bin

  TH1D *alpha = new TH1D("alpha", "alpha", 3, 0, 3);
  alpha->SetCanExtend(TH1::kAllAxes);
  alpha->SetMarkerStyle(kPlus);
  alpha->SetMarkerColor(kBlue);
  alpha->GetXaxis()->SetTitle("Bin number");
  alpha->GetYaxis()->SetTitle("alpha");

  TH1D *alpha_control = new TH1D("alpha_control", "((sigma_yellow - sigma_green) - sigma_green*Delta_alpha) / error_alpha", 3, 0, 3);
  alpha_control->SetCanExtend(TH1::kAllAxes);
  alpha_control->SetMarkerStyle(kPlus);
  alpha_control->SetMarkerColor(kBlue);
  alpha_control->GetXaxis()->SetTitle("Bin number");
  alpha_control->GetYaxis()->SetTitle("Pull");

  TH1D *epsilon = new TH1D("epsilon", "epsilon", 3, 0, 3);
  epsilon->SetCanExtend(TH1::kAllAxes);
  epsilon->SetMarkerStyle(kPlus);
  epsilon->SetMarkerColor(kBlue);
  epsilon->GetXaxis()->SetTitle("Bin number");
  epsilon->GetYaxis()->SetTitle("epsilon");

  TH1D *epsilon_control = new TH1D("epsilon_control", "((mean_yellow - mean_green) - epsilon) / error_epsilon", 3, 0, 3);
  epsilon_control->SetCanExtend(TH1::kAllAxes);
  epsilon_control->SetMarkerStyle(kPlus);
  epsilon_control->SetMarkerColor(kBlack);
  epsilon_control->GetXaxis()->SetTitle("Bin number");
  epsilon_control->GetYaxis()->SetTitle("Pull");

  TH1D *epsilon_control2 = new TH1D("epsilon_control2", "((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze #nu #alpha)", 3, 0, 3);
  epsilon_control2->SetCanExtend(TH1::kAllAxes);
  epsilon_control2->SetMarkerStyle(kPlus);
  epsilon_control2->SetMarkerColor(kGreen);
  epsilon_control2->GetXaxis()->SetTitle("Bin number");
  epsilon_control2->GetYaxis()->SetTitle("Pull");

  TH1D *epsilon_test1 = new TH1D("epsilon_test1", " (mean_yellow - mean_green) / sigma_green", 3, 0, 3);
  epsilon_test1->SetCanExtend(TH1::kAllAxes);
  epsilon_test1->SetMarkerStyle(kPlus);
  epsilon_test1->SetMarkerColor(kBlue);
  epsilon_test1->GetXaxis()->SetTitle("Bin number");
  epsilon_test1->GetYaxis()->SetTitle("Pull");

  TH2D *epsilon_test2 = new TH2D("epsilon_test2", "epsilon_test2", 8, 0.1, 0.4, 8, -3.0, 3.0);
  epsilon_test2->GetXaxis()->SetTitle("abs((mean_yellow - mean_green) / sigma_green)");
  epsilon_test2->GetYaxis()->SetTitle("((mean_yellow - mean_green) - epsilon) / error_epsilon");

  TH1D *nu = new TH1D("nu", "nu", 3, 0, 3);
  nu->SetCanExtend(TH1::kAllAxes);
  nu->SetMarkerStyle(kPlus);
  nu->SetMarkerColor(kBlue);
  nu->GetXaxis()->SetTitle("Bin number");
  nu->GetYaxis()->SetTitle("nu");

  TH1D *nu_control = new TH1D("nu_control", "((integral_yellow - integral_green) - integral_green*delta_nu) / error_nu", 3, 0, 3);
  nu_control->SetCanExtend(TH1::kAllAxes);
  nu_control->SetMarkerStyle(kPlus);
  nu_control->SetMarkerColor(kBlue);
  nu_control->GetXaxis()->SetTitle("Bin number");
  nu_control->GetYaxis()->SetTitle("Pull");

  // Histograms for pull distributions

  TH1D *pull_nu_control = new TH1D("pull_nu_control", " ((integral_yellow - integral_green) - integral_green*Delta_nu) / error_nu ", 30, -5.0, 5.0);
  pull_nu_control->GetXaxis()->SetTitle("Pull");
  pull_nu_control->GetYaxis()->SetTitle("Events");

  string title_offset = " (input epsilon(=" + to_string(mll_offset) + ") - fitted epsilon) / error_fitted_epsilon (smear mass additively) ";
  TH1D *pull_epsilon_control = new TH1D("pull_epsilon_control", title_offset.c_str(), 30, -5.0, 5.0);
  pull_epsilon_control->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control->GetYaxis()->SetTitle("Events");

  TH1D *pull_epsilon_control1 = new TH1D("pull_epsilon_control1", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze nu) ", 30, -5.0, 5.0);
  pull_epsilon_control1->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control1->GetYaxis()->SetTitle("Events");

  TH1D *pull_epsilon_control2 = new TH1D("pull_epsilon_control2", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze nu, alpha) ", 30, -5.0, 5.0);
  pull_epsilon_control2->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control2->GetYaxis()->SetTitle("Events");

  TH1D *pull_epsilon_control3 = new TH1D("pull_epsilon_control3", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (|d_nu|, |d_alpha| <0.02) ", 30, -5.0, 5.0);
  pull_epsilon_control3->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control3->GetYaxis()->SetTitle("Events");

  TH1D *pull_alpha_control = new TH1D("pull_alpha_control", " ((sigma_yellow - sigma_green) - sigma_green*Delta_alpha) / error_alpha ", 30, -5.0, 5.0);
  pull_alpha_control->GetXaxis()->SetTitle("Pull");
  pull_alpha_control->GetYaxis()->SetTitle("Events");

  // Histograms for fit constituents for each k bin

  auto multi_hist_proj_diff_reco = mDh_diff_reco->Projection(4);
  auto multi_hist_proj_reco = mDh_reco->Projection(4);
  auto multi_hist_proj_gen = mDh_gen->Projection(4);
  auto multi_hist_proj_smear = mDh_smear->Projection(4);
  auto multi_hist_proj_diff_smear = mDh_diff_smear->Projection(4);
  auto multi_hist_proj_jac_diff_squared_smear_mll = mDh_jac_diff_squared_smear_mll->Projection(4);
  auto multi_hist_proj_jac_diff_squared_smear_mll_diff = mDh_jac_diff_squared_smear_mll_diff->Projection(4);
  auto multi_hist_proj_jac_diff_smear_mll_diff = mDh_jac_diff_smear_mll_diff->Projection(4);
  auto multi_hist_proj_jac_diff_smear_mll = mDh_jac_diff_smear_mll->Projection(4);
  auto multi_hist_proj_smear_beta_val = mDh_smear_beta_val->Projection(4);
  auto multi_hist_proj_diff_smear_beta_val = mDh_diff_smear_beta_val->Projection(4);
  auto multi_hist_proj_diff_smear_beta_val_easy = mDh_diff_smear_beta_val_easy->Projection(4);
  auto multi_hist_proj_diff_smear_plus_offset = mDh_diff_smear_plus_offset->Projection(4);
  auto multi_hist_proj_diff_smear_minus_offset = mDh_diff_smear_minus_offset->Projection(4);

  //--------------------------------------------------------------------------------------

  // Files to write results
  std::unique_ptr<TFile> f_control( TFile::Open("control_bin_histo.root", "RECREATE") ); // histos inclusive in k bins
  std::unique_ptr<TFile> f_fits( TFile::Open("reco_gen_histos.root", "RECREATE") ); // histos per k bin
  ofstream f_pass_reg("passed_regions.txt"); // to check if there are enough k bins to constrain all sagitta parameters

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

  //--------------------------------------------------------------------------------------
  // Prepare counting events in mll_diff_reco

  total_nevents = multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff);
  std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";
  
  all_histos_count = 0.0;
  empty_histos_count = 0.0;
  hfrac = -1.0;  
  efrac = -1.0;
  remaining_nevents = 0.0;

  //--------------------------------------------------------------------------------------
  // Loop over eta+,pt+,eta-,pt- 

  for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
    mDh_diff_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_gen->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_squared_smear_mll->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_squared_smear_mll_diff->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_smear_mll_diff->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_smear_mll->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_smear_beta_val->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear_beta_val->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear_beta_val_easy->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear_plus_offset->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear_minus_offset->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    
    for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){   
      mDh_diff_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_gen->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_squared_smear_mll->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_squared_smear_mll_diff->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_smear_mll_diff->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_smear_mll->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_smear_beta_val->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear_beta_val->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear_beta_val_easy->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear_plus_offset->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear_minus_offset->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
    
      for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	mDh_diff_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_gen->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_diff_squared_smear_mll->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_jac_diff_squared_smear_mll_diff->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_diff_smear_mll_diff->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_diff_smear_mll->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_smear_beta_val->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear_beta_val->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear_beta_val_easy->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear_plus_offset->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear_minus_offset->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	
	for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	  mDh_diff_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  	  
	  all_histos_count++;
	  
	  //--------------------------------------------------------------------------
	  
	  // Require enough stats in mll_diff_reco for sigma_MC fit
	  delete gROOT->FindObject("multi_data_histo_diff_reco_proj_4");
	  multi_hist_proj_diff_reco = mDh_diff_reco->Projection(4);
	  nevents = multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff); 
	  
	  if (nevents < 100.0){ // reject low stats
	    empty_histos_count++;
	  } else {

	    // Require most of mll_reco Gaussian to be in fit window
	    delete gROOT->FindObject("multi_data_histo_reco_proj_4");
	    multi_hist_proj_reco = mDh_reco->Projection(4);
	    multi_hist_proj_reco->GetXaxis()->SetTitle("mll [GeV]");
	    multi_hist_proj_reco->GetYaxis()->SetTitle("Events");
	    fitresult = fitHisto(multi_hist_proj_reco, 0, 4); 
	    
	    if (fitresult[4] < 0.75){ // reject small gaus integral 
	      empty_histos_count++;
	    } else {

	      remaining_nevents += nevents;	      
	      occupancy->SetBinError(occupancy->Fill(name.c_str(), nevents), 100);
	      
	      name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	      //std::cout<< name <<"\n";
	      f_pass_reg << name << "\n"; 
	      
	      //--------------------------------------------------------------------------
	      // Set range of the remaining histograms
	      mDh_gen->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_diff_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_squared_smear_mll->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_squared_smear_mll_diff->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_smear_mll_diff->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_smear_mll->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_smear_beta_val->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_diff_smear_beta_val->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_diff_smear_beta_val_easy->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_diff_smear_plus_offset->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_diff_smear_minus_offset->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      
	      //--------------------------------------------------------------------------
	      // diff_reco histogram

	      // Already projected diff_reco
	      //multi_hist_proj_diff_reco->SetName(name.c_str());
	      multi_hist_proj_diff_reco->GetXaxis()->SetTitle("mll_diff [GeV]");
	      multi_hist_proj_diff_reco->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_diff_reco->SetLineColor(kBlue);
	      
	      // Fit diff_reco
	      fitresult = fitHisto(multi_hist_proj_diff_reco, 1, 4);
	      if (fitresult[0] > -90.0){ mean_reco->SetBinError(mean_reco->Fill(name.c_str(), fitresult[0]), fitresult[2]); }
	      if (fitresult[1] > -5.0){ sigma_reco->SetBinError(sigma_reco->Fill(name.c_str(), fitresult[1]), fitresult[3]); }
	      	      
	      // Prepare to draw mll_diff panel
	      c1->cd(1); 
	      max_hist_mll_diff = -1.0;	      
	      //max_hist_mll_diff = multi_hist_proj_diff_reco->GetBinContent(multi_hist_proj_diff_reco->GetMaximumBin());
	      
	      //--------------------------------------------------------------------------
	      // diff_smear histogram

	      // Project diff_smear
	      delete gROOT->FindObject("multi_data_histo_diff_smear_proj_4");
	      multi_hist_proj_diff_smear = mDh_diff_smear->Projection(4);
	      multi_hist_proj_diff_smear->GetXaxis()->SetTitle("mll_diff [GeV]");
	      multi_hist_proj_diff_smear->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_diff_smear->SetLineColor(kGreen);
	      
	      // Fit diff_smear
	      fitresult = fitHisto(multi_hist_proj_diff_smear, 1, 8);
	      if (fitresult[0] > -90.0){ mean_smear->SetBinError(mean_smear->Fill(name.c_str(), fitresult[0]), fitresult[2]); }
	      if (fitresult[1] > -5.0){ sigma_smear->SetBinError(sigma_smear->Fill(name.c_str(), fitresult[1]), fitresult[3]); }
	      
	      // TODO only carry on if the 2 ifs above are good 

	      // Save for nu, epsilon, alpha fit
	      mean_mc = fitresult[0]; 
	      error_mean_mc = fitresult[2];
	      sigma_mc = fitresult[1]; 
	      error_sigma_mc = fitresult[3]; 

	      if(max_hist_mll_diff < multi_hist_proj_diff_smear->GetBinContent(multi_hist_proj_diff_smear->GetMaximumBin())){
                max_hist_mll_diff = multi_hist_proj_diff_smear->GetBinContent(multi_hist_proj_diff_smear->GetMaximumBin());
              }

	      //--------------------------------------------------------------------------
	      // diff_smear_beta_val histogram 
	     
	      // Project diff_smear_beta_val histogram
	      delete gROOT->FindObject("multi_data_histo_diff_smear_beta_val_proj_4");
	      multi_hist_proj_diff_smear_beta_val = mDh_diff_smear_beta_val->Projection(4);
	      multi_hist_proj_diff_smear_beta_val->GetXaxis()->SetTitle("mll_diff [GeV]");
              multi_hist_proj_diff_smear_beta_val->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_diff_smear_beta_val->SetLineColor(kYellow);

	      // Fit diff_smear_beta_val
	      fitresult = fitHisto(multi_hist_proj_diff_smear_beta_val, 1, 5);
	      
	      // Save these for pull histograms 
	      mean_diff_smear_epsilon_val = fitresult[0]; 
	      error_mean_diff_smear_epsilon_val = fitresult[2];
	      sigma_diff_smear_epsilon_val = fitresult[1];
	      error_sigma_diff_smear_epsilon_val = fitresult[3];

	      if(max_hist_mll_diff < multi_hist_proj_diff_smear_beta_val->GetBinContent(multi_hist_proj_diff_smear_beta_val->GetMaximumBin())){
                max_hist_mll_diff = multi_hist_proj_diff_smear_beta_val->GetBinContent(multi_hist_proj_diff_smear_beta_val->GetMaximumBin());
              }

	      //--------------------------------------------------------------------------
	      // diff_smear_beta_val_easy histogram

	      // Project diff_smear_beta_val_easy
	      delete gROOT->FindObject("multi_data_histo_diff_smear_beta_val_easy_proj_4");
	      multi_hist_proj_diff_smear_beta_val_easy = mDh_diff_smear_beta_val_easy->Projection(4);
	      multi_hist_proj_diff_smear_beta_val_easy->GetXaxis()->SetTitle("mll_diff [GeV]");
              multi_hist_proj_diff_smear_beta_val_easy->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_diff_smear_beta_val_easy->SetLineColor(kPink);
	      
	      // Fit diff_smear_beta_val_easy
	      fitresult = fitHisto(multi_hist_proj_diff_smear_beta_val_easy, 1, 6);
	      
	      //if(max_hist_mll_diff < multi_hist_proj_diff_smear_beta_val_easy->GetBinContent(multi_hist_proj_diff_smear_beta_val_easy->GetMaximumBin())){
	      //max_hist_mll_diff = multi_hist_proj_diff_smear_beta_val_easy->GetBinContent(multi_hist_proj_diff_smear_beta_val_easy->GetMaximumBin());
              //}

	      //--------------------------------------------------------------------------
	      // diff_smear_plus_offset histogram

              // Project diff_smear_plus_offset
              delete gROOT->FindObject("multi_data_histo_diff_smear_plus_offset_proj_4");
              multi_hist_proj_diff_smear_plus_offset = mDh_diff_smear_plus_offset->Projection(4);
	      multi_hist_proj_diff_smear_plus_offset->SetLineColor(kRed);

              // Fit diff_smear_plus_offset
              fitresult = fitHisto(multi_hist_proj_diff_smear_plus_offset, 1, 2);

              if(max_hist_mll_diff < multi_hist_proj_diff_smear_plus_offset->GetBinContent(multi_hist_proj_diff_smear_plus_offset->GetMaximumBin())){
              max_hist_mll_diff = multi_hist_proj_diff_smear_plus_offset->GetBinContent(multi_hist_proj_diff_smear_plus_offset->GetMaximumBin());
              }

	      //--------------------------------------------------------------------------
              // diff_smear_minus_offset histogram

              // Project diff_smear_minus_offset
              delete gROOT->FindObject("multi_data_histo_diff_smear_minus_offset_proj_4");
              multi_hist_proj_diff_smear_minus_offset = mDh_diff_smear_minus_offset->Projection(4);
              multi_hist_proj_diff_smear_minus_offset->SetLineColor(kBlue);

              // Fit diff_smear_minus_offset
              fitresult = fitHisto(multi_hist_proj_diff_smear_minus_offset, 1, 4);

              if(max_hist_mll_diff < multi_hist_proj_diff_smear_minus_offset->GetBinContent(multi_hist_proj_diff_smear_minus_offset->GetMaximumBin())){
		max_hist_mll_diff = multi_hist_proj_diff_smear_minus_offset->GetBinContent(multi_hist_proj_diff_smear_minus_offset->GetMaximumBin());
              }

	      //--------------------------------------------------------------------------
	      // Start drawing mll_diff 
	      	      
	      max_hist_mll_diff = max_hist_mll_diff * 1.3;
              multi_hist_proj_diff_smear->SetMaximum(max_hist_mll_diff);
	      multi_hist_proj_diff_smear->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      multi_hist_proj_diff_smear->Draw();
	      multi_hist_proj_diff_smear_beta_val->Draw("SAME");
	      multi_hist_proj_diff_smear_plus_offset->Draw("SAME");
	      multi_hist_proj_diff_smear_minus_offset->Draw("SAME");
	      //multi_hist_proj_diff_smear_beta_val_easy->Draw("SAME");
	      //multi_hist_proj_diff_reco->Draw("SAME");

	      // Legend

	      leg1->Clear();
	      leg_entry = "Region " + name + ": " + nevents + " events";
	      leg1->SetHeader(leg_entry.c_str(),"C");
              
	      //leg1->AddEntry(multi_hist_proj_diff_reco, "reco-gen", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear, "smeared-gen, #beta_pT=1, #alpha=1", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear_beta_val, "smeared-gen, #beta_pT=0.999, #alpha=1, by weight", "l");
	      //leg1->AddEntry(multi_hist_proj_diff_smear_beta_val_easy, "smeared-gen, #varepsilon_pT=-0.15, #Delta#alpha=0, by sampling", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear_plus_offset, "diff_smear_plus_offset", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear_minus_offset, "diff_smear_minus_offset", "l");
	      leg_entry = "green fit: #mu=" + to_string(mean_mc).substr(0, 5) + ", #sigma=" + to_string(sigma_mc).substr(0, 5);
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      //leg1->Draw("");
	      
	      //--------------------------------------------------------------------------
	      // Prepare to draw mll panel

	      c1->cd();
	      c1->cd(2); 
	      max_hist_mll = -1.0;
	      	      
	      //--------------------------------------------------------------------------
	      // mass reco histogram

	      // Already projected mass reco
	      multi_hist_proj_reco->SetLineColor(kBlue);
	      multi_hist_proj_reco->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      
	      // Fit mass reco
	      fitresult = fitHisto(multi_hist_proj_reco, 0, 4);
	      
	      if (fitresult[4] > -100.0){ gaus_integral->SetBinError(gaus_integral->Fill(name.c_str(), fitresult[4]), 0.01); }

	      max_hist_mll = multi_hist_proj_reco->GetBinContent(multi_hist_proj_reco->GetMaximumBin());

	      //--------------------------------------------------------------------------
              // mass gen histogram

	      // Project mass gen
	      delete gROOT->FindObject("multi_data_histo_gen_proj_4");
	      multi_hist_proj_gen = mDh_gen->Projection(4);
	      multi_hist_proj_gen->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      multi_hist_proj_gen->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_gen->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_gen->SetLineColor(kRed);
	      
	      // Fit mass gen
	      //fitresult = fitHisto(multi_hist_proj_gen, 0, 2);
	      
	      if(max_hist_mll < multi_hist_proj_gen->GetBinContent(multi_hist_proj_gen->GetMaximumBin())){
		max_hist_mll = multi_hist_proj_gen->GetBinContent(multi_hist_proj_gen->GetMaximumBin());
	      }
	      
	      //--------------------------------------------------------------------------
              // mass smear histogram

	      // Project mass smear
	      delete gROOT->FindObject("multi_data_histo_smear_proj_4");
	      multi_hist_proj_smear = mDh_smear->Projection(4);
	      multi_hist_proj_smear->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      multi_hist_proj_smear->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_smear->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_smear->SetLineColor(kGreen);
	    
	      // Fit mass smear
	      //fitresult = fitHisto(multi_hist_proj_smear, 0, 8);

	      if(max_hist_mll < multi_hist_proj_smear->GetBinContent(multi_hist_proj_smear->GetMaximumBin())){
                max_hist_mll = multi_hist_proj_smear->GetBinContent(multi_hist_proj_smear->GetMaximumBin());
	      }

	      //--------------------------------------------------------------------------
              // mass smear_beta_val
	      
	      // Project mass smear_beta_val
	      delete gROOT->FindObject("multi_data_histo_smear_beta_val_proj_4");
              multi_hist_proj_smear_beta_val = mDh_smear_beta_val->Projection(4);
	      multi_hist_proj_smear_beta_val->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_smear_beta_val->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_smear_beta_val->SetLineColor(kYellow);

	      // Fit mass smear_beta_val
	      //fitresult = fitHisto(multi_hist_proj_smear_beta_val, 0, 5);
              
	      if(max_hist_mll < multi_hist_proj_smear_beta_val->GetBinContent(multi_hist_proj_smear_beta_val->GetMaximumBin())){
                max_hist_mll = multi_hist_proj_smear_beta_val->GetBinContent(multi_hist_proj_smear_beta_val->GetMaximumBin());
	      }

	      //--------------------------------------------------------------------------
              // Start drawing mll
	      
	      max_hist_mll = max_hist_mll * 1.3;
	      multi_hist_proj_gen->SetMaximum(max_hist_mll);

	      //multi_hist_proj_reco->Draw();
	      //multi_hist_proj_gen->Draw();
	      multi_hist_proj_smear->Draw("");
	      multi_hist_proj_smear_beta_val->Draw("SAME");

	      // Legend 
	      leg2->Clear();
	      leg_entry = "Region " + name + ": " + nevents + " events";
	      leg2->SetHeader(leg_entry.c_str(),"C"); 

	      //leg2->AddEntry(multi_hist_proj_reco, "reco", "l");
	      //leg2->AddEntry(multi_hist_proj_gen, "gen", "l");
	      leg2->AddEntry(multi_hist_proj_smear, "smeared, #beta=1", "l");
	      leg2->AddEntry(multi_hist_proj_smear_beta_val, "smeared, #beta=0.999", "l"); // insert value of beta here
	      
	      c1->cd();
	      
	      //--------------------------------------------------------------------------
	      // Fill vectors and variance for minimisation
	      //--------------------------------------------------------------------------
	      
	      //--------------------------------------------------------------------------
	      // Find mass bins with defined jacobian
	      
	      filled_bins_mll=0;
	      vector<int> good_indices_mll;
	      for(int i=1; i<=nbinsmll; i++){
		// TODO the requirement for >10 evts in jac might need to be different for numerical jac due to bin migration
		if ( multi_hist_proj_smear->GetBinContent(i) >= 10 && multi_hist_proj_smear_beta_val->GetBinContent(i) > 0 ){ //asking jacobian to be average over at least 10 events, TODO refine the other criterium
		  good_indices_mll.push_back(i);
		  filled_bins_mll++;
		}
	      }
	      if (filled_bins_mll < 3 ){continue;} // Need minimum 3 data points to fit nu, epsilon, alpha TODO refine this continue, e.g maybe i want to ocntinue to mll_diff

	      // Find mll_diff bins with defined jacobian
	      
	      filled_bins_mll_diff=0;
	      vector<int> good_indices_mll_diff;
	      for(int i=1; i<=nbinsmll_diff; i++){
		// TODO the requirement for >10 evts in jac might need to be different for numerical jac due to bin migration
		if( multi_hist_proj_diff_smear->GetBinContent(i) >= 10 && multi_hist_proj_diff_smear_beta_val->GetBinContent(i) > 0 ){ //asking jacobian to be average over at least 10 events, TODO refine the other criterium
		  good_indices_mll_diff.push_back(i);
		  filled_bins_mll_diff++;
		}
	      }
	      if (filled_bins_mll_diff < 3 ){continue;} // Need minimum 3 data points to fit nu, epsilon, alpha TODO refine this continue
	      
	      //--------------------------------------------------------------------------
	      // Declare mass vectors, jacobians and variance
	      
	      VectorXd h_smear_minus_smear_vector(filled_bins_mll);
	      Eigen::MatrixXd V_sqrt(filled_bins_mll, filled_bins_mll), V_inv_sqrt(filled_bins_mll, filled_bins_mll), J(filled_bins_mll, 3); //J.col(0) is nu, (1) is epsilon, (2) is alpha

	      V_sqrt = MatrixXd::Zero(filled_bins_mll, filled_bins_mll); 
	      position_to_fill=0;
              for(int i=1; i<=nbinsmll; i++){
		if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
		  h_smear_minus_smear_vector(position_to_fill) = multi_hist_proj_smear_beta_val->GetBinContent(i) - multi_hist_proj_smear->GetBinContent(i);
		  J(position_to_fill,0) = multi_hist_proj_smear->GetBinContent(i); //J_control.col(0) is nu
		  V_sqrt(position_to_fill,position_to_fill) = pow(pow(multi_hist_proj_smear_beta_val->GetBinError(i),2) + pow(multi_hist_proj_smear->GetBinError(i),2), 0.5); // sqrt(data_stat**2 + mc_stat**2) 
		  //V_inv_sqrt(position_to_fill,position_to_fill) = 1 / multi_hist_proj_smear_beta_val->GetBinError(i);
		  position_to_fill++;
		}
	      }
	      if (position_to_fill != filled_bins_mll){ std::cout<<"problem counting vector size \n"; }
	      // solve for V_inv_sqrt
	      V_inv_sqrt = V_sqrt.completeOrthogonalDecomposition().solve(MatrixXd::Identity(filled_bins_mll,filled_bins_mll));

	      // Declare mll_diff vectors, jacobians and variance
              VectorXd h_smear_diff_minus_smear_diff_vector(filled_bins_mll_diff);
	      Eigen::MatrixXd V_sqrt_control(filled_bins_mll_diff, filled_bins_mll_diff), V_inv_sqrt_control(filled_bins_mll_diff, filled_bins_mll_diff), J_control(filled_bins_mll_diff, 3); //J_control.col(0) is nu, (1) is epsilon, (2) is alpha 

	      V_sqrt_control = MatrixXd::Zero(filled_bins_mll_diff, filled_bins_mll_diff);
	      position_to_fill = 0;
	      integral_diff_smear_epsilon_val = 0.0;
	      integral_mc = 0.0;

	      for(int i=1; i<=nbinsmll_diff; i++){
		if( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		  h_smear_diff_minus_smear_diff_vector(position_to_fill) = multi_hist_proj_diff_smear_beta_val->GetBinContent(i) - multi_hist_proj_diff_smear->GetBinContent(i);
		  J_control(position_to_fill,0) = multi_hist_proj_diff_smear->GetBinContent(i); //J_control.col(0) is nu
		  V_sqrt_control(position_to_fill,position_to_fill) = pow(pow(multi_hist_proj_diff_smear_beta_val->GetBinError(i),2) + pow(multi_hist_proj_diff_smear->GetBinError(i),2), 0.5); // sqrt(data_stat**2 + mc_stat**2)
		  //V_inv_sqrt_control(position_to_fill,position_to_fill) = 1 / multi_hist_proj_diff_smear_beta_val->GetBinError(i);
		  integral_diff_smear_epsilon_val += multi_hist_proj_diff_smear_beta_val->GetBinContent(i);
		  integral_mc += multi_hist_proj_diff_smear->GetBinContent(i);

		  position_to_fill++;
		}
	      }
	      if (position_to_fill != filled_bins_mll_diff){ std::cout<<"problem counting vector size \n"; }
	      // solve for V_inv_sqrt_control
	      V_inv_sqrt_control = V_sqrt_control.completeOrthogonalDecomposition().solve(MatrixXd::Identity(filled_bins_mll_diff,filled_bins_mll_diff));

	      //--------------------------------------------------------------------------
	      // Define jacobian histograms

	      // Jacobian histogram alpha mass smear
	      delete gROOT->FindObject("jacobian");
	      TH1D *jac = new TH1D("jacobian", "jacobian", nbinsmll, 75.0, 105.0);
              jac->GetXaxis()->SetTitle("mll_smear [GeV]");
              jac->GetYaxis()->SetTitle("GeV");

	      // Jacobian histogram alpha diff_smear 
	      delete gROOT->FindObject("jacobian_control");
              TH1D *jac_control = new TH1D("jacobian_control", "jacobian", nbinsmll_diff, -5.0, 5.0);
              jac_control->GetXaxis()->SetTitle("mll_smear(#beta=1) - mll_gen [GeV]");
              jac_control->GetYaxis()->SetTitle("GeV");

	      // Jacobian histogram epsilon mass smear
	      delete gROOT->FindObject("jacobian_epsilon");
	      TH1D *jac_epsilon = new TH1D("jacobian_epsilon", "jacobian epsilon", nbinsmll, 75.0, 105.0);
              jac_epsilon->GetXaxis()->SetTitle("mll_smear [GeV]");
              jac_epsilon->GetYaxis()->SetTitle("GeV");

	      // Jacobian histogram epsilon diff_smear
	      delete gROOT->FindObject("jacobian_epsilon_control");
	      TH1D *jac_epsilon_control = new TH1D("jacobian_epsilon_control", "jacobian epsilon", nbinsmll_diff, -5.0, 5.0);
              jac_epsilon_control->GetXaxis()->SetTitle("mll_smear(#epsilon=0) - mll_gen [GeV]");
              jac_epsilon_control->GetYaxis()->SetTitle("GeV");

	      // Jacobian numerical histogram epsilon diff_smear
	      delete gROOT->FindObject("jacobian_epsilon_control_num_plus");
	      delete gROOT->FindObject("jacobian_epsilon_control_num_minus");
	      delete gROOT->FindObject("jacobian_epsilon_control_num");
	      //TH1D *jac_epsilon_control_num_plus = new TH1D("jacobian_epsilon_control_num_plus", "jacobian epsilon_num_plus", nbinsmll_diff, -5.0, 5.0);
	      //TH1D *jac_epsilon_control_num_minus = new TH1D("jacobian_epsilon_control_num_minus", "jacobian epsilon_num_minus", nbinsmll_diff, -5.0, 5.0);
	      TH1D *jac_epsilon_control_num = new TH1D("jacobian_epsilon_control_num", "jacobian epsilon_num", nbinsmll_diff, -5.0, 5.0);
	      jac_epsilon_control_num->GetXaxis()->SetTitle("mll_smear(#epsilon=0) - mll_gen [GeV]");
              jac_epsilon_control_num->GetYaxis()->SetTitle("GeV");

	      //--------------------------------------------------------------------------
              // Compute jacobians mass smear

	      // Needed for jacobian alpha mass smear
	      delete gROOT->FindObject("multi_data_histo_jac_diff_squared_smear_mll_proj_4");
	      multi_hist_proj_jac_diff_squared_smear_mll = mDh_jac_diff_squared_smear_mll->Projection(4);
	      // Needed for jacobian epsilon mass smear
	      delete gROOT->FindObject("multi_data_histo_jac_diff_smear_mll_proj_4");
	      multi_hist_proj_jac_diff_smear_mll = mDh_jac_diff_smear_mll->Projection(4);
           	      
	      position_to_fill=0;
	      for(int i=1; i<=nbinsmll; i++){
		// alpha mass smear
		diff_squared = multi_hist_proj_jac_diff_squared_smear_mll->GetBinContent(i);
		error_diff_squared = multi_hist_proj_jac_diff_squared_smear_mll->GetBinError(i);
		// epsilon mass smear
		jac_e_weight = multi_hist_proj_jac_diff_smear_mll->GetBinContent(i);
                error_jac_e_weight = multi_hist_proj_jac_diff_smear_mll->GetBinError(i);
		
		if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
		  
		  evts_in_bin = multi_hist_proj_smear->GetBinContent(i);
		  // alpha mass smear
		  value =  diff_squared / (sigma_mc * sigma_mc) - evts_in_bin;
                  error = 1 / (sigma_mc*sigma_mc) * pow((4.0 * diff_squared*diff_squared * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_diff_squared*error_diff_squared) , 0.5);
		  J(position_to_fill,2) = value;
		  
		  jac->SetBinContent(i, value);
		  jac->SetBinError(i, error);
		  jac_inclusive->Fill(mllbinranges[i-1], value);
		  
		  // epsilon mass smear
		  value_epsilon = jac_e_weight / (sigma_mc * sigma_mc);
		  error_epsilon = 1 / (sigma_mc*sigma_mc) * pow((4 * jac_e_weight*jac_e_weight * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_jac_e_weight*error_jac_e_weight) , 0.5);
		  J(position_to_fill,1) = value_epsilon;
		  
		  jac_epsilon->SetBinContent(i, value_epsilon);
		  jac_epsilon->SetBinError(i, error_epsilon);
		  jac_epsilon_inclusive->Fill(mllbinranges[i-1], value_epsilon);
		  
		  position_to_fill++;
		  
		} 
	      }
	      
	      if (position_to_fill != filled_bins_mll){ std::cout<<"PROBLEM: counting jac mass smear size \n"; }
	      
	      //--------------------------------------------------------------------------
              // Compute jacobians diff_smear

	      // Needed for jacobian alpha diff_smear
	      delete gROOT->FindObject("multi_data_histo_jac_diff_squared_smear_mll_diff_proj_4");
	      multi_hist_proj_jac_diff_squared_smear_mll_diff = mDh_jac_diff_squared_smear_mll_diff->Projection(4);
              delete gROOT->FindObject("multi_data_histo_jac_diff_smear_mll_diff_proj_4");
              multi_hist_proj_jac_diff_smear_mll_diff = mDh_jac_diff_smear_mll_diff->Projection(4);
	      // Needed for jacobian epsilon diff_smear
	      delete gROOT->FindObject("multi_data_histo_jac_diff_smear_mll_diff_proj_4");
	      multi_hist_proj_jac_diff_smear_mll_diff = mDh_jac_diff_smear_mll_diff->Projection(4);
	      // Needed for jacobian numerical epsilon diff_smear
	      //delete gROOT->FindObject("multi_data_histo_diff_smear_plus_offset_proj_4");
	      //multi_hist_proj_diff_smear_plus_offset = mDh_diff_smear_plus_offset->Projection(4);
	      //delete gROOT->FindObject("multi_data_histo_diff_smear_minus_offset_proj_4");
              //multi_hist_proj_diff_smear_minus_offset = mDh_diff_smear_minus_offset->Projection(4);

	      position_to_fill = 0;
              for(int i=1; i<=nbinsmll_diff; i++){
		evts_in_bin = multi_hist_proj_diff_smear->GetBinContent(i);
		// alpha diff_smear
		diff_squared = multi_hist_proj_jac_diff_squared_smear_mll_diff->GetBinContent(i);
                error_diff_squared = multi_hist_proj_jac_diff_squared_smear_mll_diff->GetBinError(i);
                diff = multi_hist_proj_jac_diff_smear_mll_diff->GetBinContent(i);
                error_diff = multi_hist_proj_jac_diff_smear_mll_diff->GetBinError(i);
		// epsilon diff_smear
		// TODO just use diff variable
                jac_e_weight = multi_hist_proj_jac_diff_smear_mll_diff->GetBinContent(i);
                error_jac_e_weight = multi_hist_proj_jac_diff_smear_mll_diff->GetBinError(i);
                
                if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		  // alpha diff_smear
		  value =  (diff_squared - 2.0*mean_mc*diff + evts_in_bin*mean_mc*mean_mc) / (sigma_mc * sigma_mc) - evts_in_bin;
		  error = 2.0 / (sigma_mc*sigma_mc) * pow( pow((diff_squared-2.0*mean_mc*diff + evts_in_bin*mean_mc*mean_mc)*error_sigma_mc/sigma_mc,2) + pow(mean_mc*error_diff,2) + pow(error_diff_squared,2)/4.0 + pow((evts_in_bin*mean_mc-diff)*error_mean_mc,2)  ,0.5);

		  J_control(position_to_fill,2) = value; //J_control.col(2) is for alpha
		  jac_control->SetBinContent(i, value);
                  jac_control->SetBinError(i, error);
		  
		  // epsilon diff_smear
		  // TODO check error and formula
		  value_epsilon = ( jac_e_weight - evts_in_bin*mean_mc ) / pow(sigma_mc,2);
		  error_epsilon = pow( pow(evts_in_bin*error_mean_mc,2) + 4.0*pow(jac_e_weight - evts_in_bin*mean_mc,2)*pow(error_sigma_mc,2)/pow(sigma_mc,2) + pow(error_jac_e_weight,2), 0.5) / pow(sigma_mc,2);
	
		  J_control(position_to_fill,1) = value_epsilon; //J_control.col(1) is for epsilon
		  jac_epsilon_control->SetBinContent(i, value_epsilon);
		  jac_epsilon_control->SetBinError(i, error_epsilon);

		  // epsilon diff_smear numerical
		  value_epsilon = (multi_hist_proj_diff_smear_plus_offset->GetBinContent(i) - multi_hist_proj_diff_smear_minus_offset->GetBinContent(i)) / (mll_offset * 2.0);
		  error_epsilon = pow((pow(multi_hist_proj_diff_smear_plus_offset->GetBinError(i),2) + pow(multi_hist_proj_diff_smear_minus_offset->GetBinError(i),2)) / pow(mll_offset*2.0, 2), 0.5);
		  //---------------------------------------------//
		  //---------------------------------------------//
		  // ATTENTION OVERWRITING J_control for the moment
		  J_control(position_to_fill,1) = value_epsilon;
		  //---------------------------------------------//
		  //---------------------------------------------//
		  jac_epsilon_control_num->SetBinContent(i, value_epsilon);
                  jac_epsilon_control_num->SetBinError(i, error_epsilon);

                  position_to_fill++;

                } 
	      }

              if (position_to_fill != filled_bins_mll_diff){ std::cout<<"PROBLEM: counting jac diff smear size \n"; }

	      //--------------------------------------------------------------------------
              // Write jacobian histograms

	      f_fits->WriteObject(jac, ("jac_alpha" + name).c_str());
	      f_fits->WriteObject(jac_control, ("jac_alpha_control" + name).c_str());
	      f_fits->WriteObject(jac_epsilon, ("jac_epsilon" + name).c_str());
	      f_fits->WriteObject(jac_epsilon_control, ("jac_epsilon_control" + name).c_str());
	      f_fits->WriteObject(jac_epsilon_control_num, ("jac_epsilon_control_num" + name).c_str());
		      
	      //--------------------------------------------------------------------------
              // Solve for alpha, epsilon, nu
              //--------------------------------------------------------------------------
	      
	      // Solve for nu, epsilon, alpha mass smear 
	      Eigen::MatrixXd A = V_inv_sqrt*J;
	      Eigen::MatrixXd b = V_inv_sqrt*h_smear_minus_smear_vector;
	      // ATTENTION: n_e_a_vector contains nu, epsilon, alpha
	      Eigen::VectorXd n_e_a_vector = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	      	      
	      // Error on nu, epsilon, alpha mass smear 
	      Eigen::MatrixXd V_nu = (J.col(0).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      Eigen::MatrixXd V_epsilon = (J.col(1).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      fit_epsilon_error = pow(V_epsilon(0,0),0.5); // save for closure test
	      Eigen::MatrixXd V_alpha = (J.col(2).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(2)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      
	      // Write nu, epsilon, alpha mass smear
	      nu->SetBinError(nu->Fill(name.c_str(), n_e_a_vector(0)), pow(V_nu(0,0),0.5)); 
	      epsilon->SetBinError(epsilon->Fill(name.c_str(), n_e_a_vector(1)), pow(V_epsilon(0,0),0.5));
	      alpha->SetBinError(alpha->Fill(name.c_str(), n_e_a_vector(2)), pow(V_alpha(0,0),0.5));

	      //--------------------------------------------------------------------------

	      // Solve for nu, epsilon, alpha diff smear simultaneously
	      Eigen::MatrixXd A_control = V_inv_sqrt_control*J_control;
	      Eigen::MatrixXd b_control = V_inv_sqrt_control*h_smear_diff_minus_smear_diff_vector;
	      // ATTENTION: n_e_a_vector_control contains nu, epsilon, alpha
	      Eigen::VectorXd n_e_a_vector_control = A_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_control);
              
	      // Solve for epsilon, alpha diff smear only
	      Eigen::MatrixXd A1_control = V_inv_sqrt_control*J_control.rightCols(2); 
	      Eigen::MatrixXd b1_control = V_inv_sqrt_control*h_smear_diff_minus_smear_diff_vector;
              // ATTENTION: e_a_vector_control contains epsilon, alpha
	      Eigen::VectorXd e_a_vector_control = A1_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b1_control);

	      // Solve for epsilon diff smear only
	      Eigen::MatrixXd A2_control = V_inv_sqrt_control*J_control.col(1);
	      Eigen::MatrixXd b2_control = V_inv_sqrt_control*h_smear_diff_minus_smear_diff_vector;
              // ATTENTION: e_vector_control contains epsilon
	      Eigen::VectorXd e_vector_control = A2_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b2_control);
	      
	      // Error on epsilon diff smear
	      Eigen::MatrixXd V_epsilon_control(1,1); 
	      V_epsilon_control = (J_control.col(1).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      // Error on nu diff smear
	      Eigen::MatrixXd V_nu_control(1,1); 
	      V_nu_control = (J_control.col(0).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      // TODO change to pull on nu
	      nu_control->SetBinError(nu_control->Fill(name.c_str(), n_e_a_vector_control(0)), pow(V_nu_control(0,0),0.5)); 
	      // Error on alpha diff smear
	      Eigen::MatrixXd V_alpha_control(1,1);
	      V_alpha_control = (J_control.col(2).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(2)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      	      
	      //--------------------------------------------------------------------------
              // Add fitted parameters to mass fit panel 

	      c1->cd(2);
	      leg_entry = "Fit #varepsilon=" + to_string(n_e_a_vector(1)).substr(0, 6) + "#pm" + to_string(pow(V_epsilon(0,0),0.5)).substr(0, 6); 
	      leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #Delta#nu=" + to_string(n_e_a_vector(0)).substr(0, 6) + "#pm" + to_string(pow(V_nu(0,0),0.5)).substr(0, 6); 
	      leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #Delta#alpha=" + to_string(n_e_a_vector(2)).substr(0, 6) + "#pm" + to_string(pow(V_alpha(0,0),0.5)).substr(0, 6);
              leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");

	      // Add fitted parameters to diff fit panel 
              c1->cd();
	      c1->cd(1);

	      leg_entry = "Fit #varepsilon=" + to_string(n_e_a_vector_control(1)).substr(0, 6) + "#pm" + to_string(pow(V_epsilon_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(1) is epsilon
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #varepsilon=" + to_string(e_vector_control(0)).substr(0, 6) + " if freezed #nu, #alpha"; 
              leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #Delta#nu=" + to_string(n_e_a_vector_control(0)).substr(0, 6) + "#pm" + to_string(pow(V_nu_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(0) is nu 
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      leg_entry = "Fit #Delta#alpha=" + to_string(n_e_a_vector_control(2)).substr(0, 6) + "#pm" + to_string(pow(V_alpha_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(2) is alpha
              leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");

	      //--------------------------------------------------------------------------
              // Apply fitted epsilon and nu and alpha correction

	      //--------------------------------------------------------------------------
	      // Define corrected_smear histogram
	      delete gROOT->FindObject("corrected_smear");
              TH1D *corrected_smear = new TH1D("corrected_smear", "corrected_smear", nbinsmll, 75.0, 105.0);
              corrected_smear->SetTitle(("mll " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      corrected_smear->GetXaxis()->SetTitle("mll [GeV]");
              corrected_smear->GetYaxis()->SetTitle("Events");
	      
	      position_to_fill = 0;
              for(int i=1; i<=nbinsmll; i++){
                if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
		  value_corrected = multi_hist_proj_smear->GetBinContent(i) + n_e_a_vector(1)*J(position_to_fill,1) + n_e_a_vector(0)*J(position_to_fill,0) + n_e_a_vector(2)*J(position_to_fill,2);
		  error = pow(
			      pow(multi_hist_proj_smear->GetBinError(i),2) +
			      pow(n_e_a_vector(1)*jac_epsilon->GetBinError(i),2) +
			      pow(J(position_to_fill,1),2)*V_epsilon(0,0) +
			      pow(n_e_a_vector(0)*multi_hist_proj_smear->GetBinError(i),2) +
			      pow(J(position_to_fill,0),2)*V_nu(0,0) +
			      pow(n_e_a_vector(2)*jac->GetBinError(i),2) +
			      pow(J(position_to_fill,2),2)*V_alpha(0,0)
			      ,0.5);
		  
                  corrected_smear->SetBinContent(i, value_corrected);
                  corrected_smear->SetBinError(i, error);

                  position_to_fill++;
                }
              }
              if (position_to_fill != filled_bins_mll){ std::cout<<"PROBLEM: counting corrected_smear bins \n"; }

	      // Draw corrected_smear in mass fit panel
              c1->cd();
              c1->cd(2);
              corrected_smear->Draw("SAME");

	      leg2->AddEntry(corrected_smear, "corrected smeared-gen", "l");
	      leg2->Draw("");

	      //--------------------------------------------------------------------------
	      // Define corrected_diff_smear histogram
	      delete gROOT->FindObject("corrected_diff_smear");
	      TH1D *corrected_diff_smear = new TH1D("corrected_diff_smear", "corrected_diff_smear", nbinsmll_diff, -5.0, 5.0);
	      corrected_diff_smear->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      corrected_diff_smear->GetXaxis()->SetTitle("mll_diff [GeV]");
              corrected_diff_smear->GetYaxis()->SetTitle("Events");
	      
	      position_to_fill = 0;
              for(int i=1; i<=nbinsmll_diff; i++){
		if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		  value_corrected = multi_hist_proj_diff_smear->GetBinContent(i) + n_e_a_vector_control(1)*J_control(position_to_fill,1) + n_e_a_vector_control(0)*J_control(position_to_fill,0) + n_e_a_vector_control(2)*J_control(position_to_fill,2);
		  error = pow(
		            pow(multi_hist_proj_diff_smear->GetBinError(i),2) +
		            pow(n_e_a_vector_control(1)*jac_epsilon_control->GetBinError(i),2) +
		            pow(J_control(position_to_fill,1),2)*V_epsilon_control(0,0) +
		            pow(n_e_a_vector_control(0)*multi_hist_proj_diff_smear->GetBinError(i),2) + 
		            pow(J_control(position_to_fill,0),2)*V_nu_control(0,0) +
		            pow(n_e_a_vector_control(2)*jac_control->GetBinError(i),2) + 
		            pow(J_control(position_to_fill,2),2)*V_alpha_control(0,0)
			    ,0.5);
		  
		  corrected_diff_smear->SetBinContent(i, value_corrected);
		  corrected_diff_smear->SetBinError(i, error);

		  position_to_fill++;
		}
	      }
	      if (position_to_fill != filled_bins_mll_diff){ std::cout<<"PROBLEM: counting corrected_diff_smear bins \n"; }
	      
	      // Draw corrected_diff_smear in diff fit panel
              c1->cd();
              c1->cd(1);
              
	      fitresult = fitHisto(corrected_diff_smear, 1, 1);
	      
	      // Save for pull distribution epsilon_control
	      mean_corrected_diff_smear = fitresult[0];
	      sigma_corrected_diff_smear = fitresult[1];
	      error_mean_corrected_diff_smear = fitresult[2];
	      error_sigma_corrected_diff_smear = fitresult[3];
              corrected_diff_smear->SetLineColor(kBlack);
	      corrected_diff_smear->Draw("SAME");

	      //--------------------------------------------------------------------------
	      // Finish drawing diff panel
	      leg1->AddEntry(corrected_diff_smear, "corrected smeared-gen", "l");
	      leg1->Draw("");
	      
	      c1->cd();
	      // Write plots per k bin
              f_fits->WriteObject(c1, name.c_str());

	      // Draw pull distributions

	      // nu diff
	      nu_control->SetBinError(nu_control->Fill(name.c_str(), (integral_diff_smear_epsilon_val - integral_mc - integral_mc*n_e_a_vector_control(0)) / integral_mc / pow(V_nu_control(0,0),0.5) ), 0.001);
	      pull_nu_control->Fill( (integral_diff_smear_epsilon_val - integral_mc - integral_mc*n_e_a_vector_control(0)) / integral_mc / pow(V_nu_control(0,0),0.5) );
	      
	      // epsilon diff
	      epsilon_control->SetBinError(epsilon_control->Fill(name.c_str(), (mean_diff_smear_epsilon_val - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) ), 0.001);
	      //pull_epsilon_control->Fill( (mean_diff_smear_epsilon_val - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) );
	      pull_epsilon_control->Fill( (mll_offset - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) );
	      
	      if (abs(n_e_a_vector_control(0)) < 0.02 && abs(n_e_a_vector_control(2)) < 0.02){
		pull_epsilon_control3->Fill( (mll_offset - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) );
	      }
	      
	      pull_epsilon_control1->Fill( (mll_offset - e_a_vector_control(0))/pow(V_epsilon_control(0,0),0.5) );
	      pull_epsilon_control2->Fill( (mll_offset - e_vector_control(0))/pow(V_epsilon_control(0,0),0.5) );

	      epsilon_control2->SetBinError(epsilon_control2->Fill(name.c_str(),(mean_diff_smear_epsilon_val - mean_mc - e_vector_control(0))/pow(V_epsilon_control(0,0),0.5)), 0.001);

	      // epsilon test
	      epsilon_test1->SetBinError(epsilon_test1->Fill(name.c_str(), (mean_diff_smear_epsilon_val - mean_mc) / sigma_mc), 0.001);
	      if ( abs((mean_diff_smear_epsilon_val - mean_mc) / sigma_mc) > 0.05 ){
		epsilon_test2->Fill(abs((mean_diff_smear_epsilon_val - mean_mc) / sigma_mc), (mean_diff_smear_epsilon_val - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5));
	      }
	      
	      // alpha diff
	      alpha_control->SetBinError(alpha_control->Fill(name.c_str(), (sigma_diff_smear_epsilon_val - sigma_mc - sigma_mc*n_e_a_vector_control(2))/pow(V_alpha_control(0,0),0.5) ), 0.001);
	      pull_alpha_control->Fill( (sigma_diff_smear_epsilon_val - sigma_mc - sigma_mc*n_e_a_vector_control(2))/pow(V_alpha_control(0,0),0.5) );
	      
	      
	    }
	  }
	  
	} 
      }
    }
  }
  f_pass_reg.close();

  hfrac = empty_histos_count / all_histos_count;
  efrac = remaining_nevents / total_nevents;
  std::cout<<"CHECKPOINT: "<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";

  //--------------------------------------------------------------------------
  // Write remaining histograms 

  mean_reco->SetStats(0);
  mean_reco->LabelsDeflate();
  f_control->WriteObject(mean_reco, "mean_diff_reco");
  
  sigma_reco->SetStats(0);
  sigma_reco->LabelsDeflate();
  f_control->WriteObject(sigma_reco, "sigma_diff_reco");

  mean_smear->SetStats(0);
  mean_smear->LabelsDeflate();
  f_control->WriteObject(mean_smear, "mean_diff_smear");

  sigma_smear->SetStats(0);
  sigma_smear->LabelsDeflate();
  f_control->WriteObject(sigma_smear, "sigma_diff_smear");

  f_control->WriteObject(alpha, "alpha");
  f_control->WriteObject(jac_inclusive, "jacobian_inclusive");

  f_control->WriteObject(epsilon, "epsilon");
  f_control->WriteObject(jac_epsilon_inclusive, "jacobian_epsilon_inclusive");

  f_control->WriteObject(nu, "nu");

  f_control->WriteObject(epsilon_control, "epsilon_control");
  fitresult = fitHisto(pull_epsilon_control, 1, 1);
  std::cout<<"CHECKPOINT: "<<"Pull epsilon control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
  f_control->WriteObject(pull_epsilon_control, "pull_epsilon_control");

  f_control->WriteObject(epsilon_test1, "epsilon_test1");

  f_control->WriteObject(alpha_control, "alpha_control");
  fitresult = fitHisto(pull_alpha_control, 1, 1);
  std::cout<<"CHECKPOINT: "<<"Pull alpha control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
  f_control->WriteObject(pull_alpha_control, "pull_alpha_control");

  f_control->WriteObject(nu_control, "nu_control");
  fitresult = fitHisto(pull_nu_control, 1, 1);
  std::cout<<"CHECKPOINT: "<<"Pull nu control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
  f_control->WriteObject(pull_nu_control, "pull_nu_control");

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

  f_control->WriteObject(c2, "mean_superimposed");

  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  epsilon_test2->SetStats(0);
  epsilon_test2->Draw("COLZ");
  f_control->WriteObject(c3, "epsilon_test2");

  //-----------------------------------------------------------
  // pull_epsilon_many

  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  auto leg4 = new TLegend(0.68, 0.78, 0.90, 0.90);

  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.025);
  leg4->SetFillColor(10);
  leg4->SetNColumns(1);
  leg4->SetHeader("");

  pull_epsilon_control->SetMaximum( pull_epsilon_control->GetBinContent(pull_epsilon_control->GetMaximumBin())*1.3 );
  pull_epsilon_control->SetStats(0);
  fitresult = fitHisto(pull_epsilon_control, 1, 1);
  leg_entry = "freeze none, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
  leg4->AddEntry(pull_epsilon_control, leg_entry.c_str(), "l");
  pull_epsilon_control->SetLineColor(kBlack);
  pull_epsilon_control->Draw();

  fitresult = fitHisto(pull_epsilon_control1, 1, 2);
  leg_entry = "freeze nu, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
  leg4->AddEntry(pull_epsilon_control1, leg_entry.c_str(), "l");
  pull_epsilon_control1->SetLineColor(kRed);
  pull_epsilon_control1->Draw("SAME");

  fitresult = fitHisto(pull_epsilon_control2, 1, 3);
  leg_entry = "freeze nu, alpha, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
  leg4->AddEntry(pull_epsilon_control2, leg_entry.c_str(), "l");
  pull_epsilon_control2->SetLineColor(kGreen);
  pull_epsilon_control2->Draw("SAME");

  fitresult = fitHisto(pull_epsilon_control3, 1, 4);
  leg_entry = "freeze none, |d_nu|, |d_alpha| <0.02, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
  leg4->AddEntry(pull_epsilon_control3, leg_entry.c_str(), "l");
  pull_epsilon_control3->SetLineColor(kBlue);
  pull_epsilon_control3->Draw("SAME");
  
  leg4->Draw("");
  
  f_control->WriteObject(c4, "pull_epsilon_many");
  //-----------------------------------------------------------

  TCanvas *c5 = new TCanvas("c5","c5",800,600);
  auto leg5 = new TLegend(0.58, 0.68, 0.90, 0.90);

  leg5->SetFillStyle(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.035);
  leg5->SetFillColor(10);
  leg5->SetNColumns(1);
  leg5->SetHeader("");

  epsilon_control->Draw();
  epsilon_control2->Draw("SAME");

  leg5->AddEntry(epsilon_control, "freeze none", "l");
  leg5->AddEntry(epsilon_control2, "freeze nu, alpha", "l");
  leg5->Draw("");
  
  f_control->WriteObject(c5, "epsilon_many");

  gaus_integral->SetStats(0);
  gaus_integral->LabelsDeflate();
  f_control->WriteObject(gaus_integral, "gaus_integral");

  occupancy->SetStats(0);
  occupancy->LabelsDeflate();
  f_control->WriteObject(occupancy, "bin_occupancy");

  return 0; 
}
