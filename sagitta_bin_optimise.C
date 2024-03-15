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

//--------------------------------------------------------------------------------------

int frame(){

  ROOT::EnableImplicitMT(128);

  //--------------------------------------------------------------------------------------

  // Choose validation/analysis mode
  string mode_option("validation"), mc_name_root, data_name_root;
  
  // Input files 
  if (mode_option.compare("analysis") == 0) {
    mc_name_root = "reco";
    data_name_root = "data2016";
  } else if (mode_option.compare("validation") == 0){
    mc_name_root = "smear";
    data_name_root = "smear_beta_val";
  }
     
  // MC 
  std::unique_ptr<TFile> myFile( TFile::Open( ("multiD_histo_"+ mc_name_root +".root").c_str() ) );
  std::unique_ptr<THnD> mDh_mll_mc(myFile->Get<THnD>( ("multi_data_histo_mll_" + mc_name_root).c_str() ) ); 
  std::unique_ptr<THnD> mDh_diff_mc(myFile->Get<THnD>( ("multi_data_histo_diff_" + mc_name_root).c_str() ) ); // reco-gen or smear-gen
  
  // Analytical jacobian terms 
  std::unique_ptr<THnD> mDh_jac_diff_squared_mll(myFile->Get<THnD>( ("multi_data_histo_jac_diff_squared_" + mc_name_root + "_mll").c_str() ) ); // sum of mll_diff_squared in mll bin
  std::unique_ptr<THnD> mDh_jac_diff_squared_mll_diff(myFile->Get<THnD>( ("multi_data_histo_jac_diff_squared_" + mc_name_root + "_mll_diff").c_str() ) ); // sum of mll_diff_squared in mll_diff bin
  std::unique_ptr<THnD> mDh_jac_diff_mll(myFile->Get<THnD>( ("multi_data_histo_jac_diff_times_gen_" + mc_name_root + "_mll").c_str() ) ); // sum of mll_diff*m_gen in mll bin //TODO make changed in var name too
  std::unique_ptr<THnD> mDh_jac_diff_mll_diff(myFile->Get<THnD>( ("multi_data_histo_jac_diff_" + mc_name_root + "_mll_diff").c_str() ) ); // sum of mll_diff in mll_diff bin
 
  // Numerical jacobian terms
  std::unique_ptr<THnD> mDh_diff_plus_offset(myFile->Get<THnD>( ("multi_data_histo_diff_" + mc_name_root + "_plus_offset").c_str() ) );
  std::unique_ptr<THnD> mDh_diff_minus_offset(myFile->Get<THnD>( ("multi_data_histo_diff_" + mc_name_root + "_minus_offset").c_str() ) );
  float mll_offset = 0.1;
  
  // Data
  
  THnD *mDh_mll_data=nullptr, *mDh_diff_data=nullptr;
  if (mode_option.compare("analysis") == 0) { 
    //std::unique_ptr<TFile> myFile2( TFile::Open( ("multiD_histo_"+ data_name_root +".root").c_str() ) );
    std::unique_ptr<TFile> myFile2( TFile::Open( ("multiD_histo_"+ mc_name_root +".root").c_str() ) );
    mDh_mll_data = myFile2->Get<THnD>( ("multi_data_histo_mll_" + data_name_root).c_str() );
  } else if (mode_option.compare("validation") == 0){
    mDh_mll_data = myFile->Get<THnD>( ("multi_data_histo_mll_" + data_name_root).c_str() );
    mDh_diff_data = myFile->Get<THnD>( ("multi_data_histo_diff_" + data_name_root).c_str() );
  }
  
  // Smear easy way
  // std::unique_ptr<TFile> myFile3( TFile::Open("multiD_histo_smear_beta_val_easy.root") );
  // std::unique_ptr<THnD> mDh_diff_smear_beta_val_easy(myFile3->Get<THnD>("multi_data_histo_diff_smear_beta_val_easy")); //it's smeared_beta_val_easy - gen  
  
  //--------------------------------------------------------------------------------------
  
  // Binning must match with 5D histogram
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=16, nbinsmll=32, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, mllbinranges, ptbinranges{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};
  
  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  std::cout<<"\n mllbinranges = [";
  for (int i=0; i<=nbinsmll; i++){mllbinranges.push_back(75.0 + i * (105.0 - 75.0)/nbinsmll); std::cout<<mllbinranges[i]<<", ";}
  std::cout<<"] \n";

  //--------------------------------------------------------------------------------------

  // Variables declaration 

  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0; // to check number of passing events
  double value, error, value_epsilon, error_epsilon, diff_squared, diff, evts_in_bin, error_diff_squared, error_diff, jac_e_weight, error_jac_e_weight; // for jac calculation
  double mean_mc, error_mean_mc, sigma_mc=0.0, error_sigma_mc, integral_mc, chi2_mc;
  double mean_diff_data, error_mean_diff_data, sigma_diff_data, error_sigma_diff_data, integral_diff_data;
  double fit_epsilon_error, value_corrected, mean_corrected_diff, error_mean_corrected_diff, sigma_corrected_diff, error_sigma_corrected_diff;
  double max_hist_mll_diff, max_hist_mll;
  int filled_bins_mll, filled_bins_mll_diff, position_to_fill;
  string name, leg_entry, title;
  vector<double> fitresult;

  // Histograms for mll_diff distribution properties

  title = mc_name_root + "-gen mll mean";
  TH1D *mean_diff = new TH1D("mean_diff", title.c_str(), 3, 0, 3);
  mean_diff->SetCanExtend(TH1::kAllAxes);
  //mean_diff->SetMarkerStyle(kPlus);
  mean_diff->GetXaxis()->SetTitle("Bin number");
  title = mc_name_root + "-gen mll mean [GeV]";
  mean_diff->GetYaxis()->SetTitle(title.c_str());

  title = mc_name_root + "-gen mll sigma";
  TH1D *sigma_diff = new TH1D("sigma_smear", title.c_str(), 3, 0, 3);
  sigma_diff->SetCanExtend(TH1::kAllAxes);
  sigma_diff->GetXaxis()->SetTitle("Bin number");
  title = mc_name_root + "-gen mll sigma [GeV]";
  sigma_diff->GetYaxis()->SetTitle(title.c_str());

  TH1D *occupancy = new TH1D("bin_occupancy", "Events in mll_diff", 3, 0, 3);
  occupancy->SetCanExtend(TH1::kAllAxes);
  occupancy->SetMarkerStyle(kPlus);
  occupancy->SetMarkerColor(kBlue);
  occupancy->SetLineColor(kBlue);
  occupancy->GetXaxis()->SetTitle("Bin number");
  occupancy->GetYaxis()->SetTitle("Events");

  // Histograms for mll distribution properties

  TH1D *gaus_integral = new TH1D("gaus_integral", "mll integral 75.5-105.0", 3, 0, 3);
  gaus_integral->SetCanExtend(TH1::kAllAxes);
  gaus_integral->SetMarkerStyle(kPlus);
  gaus_integral->SetMarkerColor(kBlue);
  gaus_integral->SetLineColor(kBlue);
  gaus_integral->GetXaxis()->SetTitle("Bin number");
  gaus_integral->GetYaxis()->SetTitle("integral [evts]");

  // Histograms for jacobians inclusive in eta, pt
  
  TH1D *jac_inclusive = new TH1D("jacobian_inclusive", "jacobian alpha inclusive in eta, pt", nbinsmll, 75.0, 105.0);
  jac_inclusive->GetXaxis()->SetTitle("mll");
  jac_inclusive->GetYaxis()->SetTitle("jacobian alpha");

  TH1D *jac_epsilon_inclusive = new TH1D("jacobian_epsilon_inclusive", "jacobian epsilon inclusive in eta, pt", nbinsmll, 75.0, 105.0);
  jac_epsilon_inclusive->GetXaxis()->SetTitle("mll");
  jac_epsilon_inclusive->GetYaxis()->SetTitle("jacobian epsilon");

  // Histograms for pull of fitted variables per k bin

  TH1D *alpha = new TH1D("alpha", "alpha", 3, 0, 3);
  alpha->SetCanExtend(TH1::kAllAxes);
  alpha->SetMarkerStyle(kPlus);
  alpha->SetMarkerColor(kBlue);
  alpha->GetXaxis()->SetTitle("Bin number");
  alpha->GetYaxis()->SetTitle("#alpha");

  TH1D *epsilon = new TH1D("epsilon", "epsilon", 3, 0, 3);
  epsilon->SetCanExtend(TH1::kAllAxes);
  epsilon->SetMarkerStyle(kPlus);
  epsilon->SetMarkerColor(kBlue);
  epsilon->GetXaxis()->SetTitle("Bin number");
  epsilon->GetYaxis()->SetTitle("#epsilon");

  TH1D *nu = new TH1D("nu", "nu", 3, 0, 3);
  nu->SetCanExtend(TH1::kAllAxes);
  nu->SetMarkerStyle(kPlus);
  nu->SetMarkerColor(kBlue);
  nu->GetXaxis()->SetTitle("Bin number");
  nu->GetYaxis()->SetTitle("#nu");

  // Validation mode
  TH1D *alpha_control=nullptr, *epsilon_control=nullptr, *nu_control=nullptr, *epsilon_test1=nullptr;
  TH2D *epsilon_test2=nullptr;

  if (mode_option.compare("validation") == 0){
  
    alpha_control = new TH1D("alpha_control", "alpha from mll_diff", 3, 0, 3);
    alpha_control->SetCanExtend(TH1::kAllAxes);
    alpha_control->SetMarkerStyle(kPlus);
    alpha_control->SetMarkerColor(kBlue);
    alpha_control->GetXaxis()->SetTitle("Bin number");
    alpha_control->GetYaxis()->SetTitle("#alpha");

    epsilon_control = new TH1D("epsilon_control", "epsilon from mll_diff", 3, 0, 3);
    epsilon_control->SetCanExtend(TH1::kAllAxes);
    epsilon_control->SetMarkerStyle(kPlus);
    epsilon_control->SetMarkerColor(kBlack);
    epsilon_control->GetXaxis()->SetTitle("Bin number");
    epsilon_control->GetYaxis()->SetTitle("#epsilon");

    nu_control = new TH1D("nu_control", "nu from mll_diff", 3, 0, 3);
    nu_control->SetCanExtend(TH1::kAllAxes);
    nu_control->SetMarkerStyle(kPlus);
    nu_control->SetMarkerColor(kBlue);
    nu_control->GetXaxis()->SetTitle("Bin number");
    nu_control->GetYaxis()->SetTitle("#nu");

    epsilon_test1 = new TH1D("epsilon_test1", " (mean_yellow - mean_green) / sigma_green", 3, 0, 3);
    epsilon_test1->SetCanExtend(TH1::kAllAxes);
    epsilon_test1->SetMarkerStyle(kPlus);
    epsilon_test1->SetMarkerColor(kBlue);
    epsilon_test1->GetXaxis()->SetTitle("Bin number");
    epsilon_test1->GetYaxis()->SetTitle("Pull");

    epsilon_test2 = new TH2D("epsilon_test2", "epsilon_test2", 8, 0.1, 0.4, 8, -3.0, 3.0);
    epsilon_test2->GetXaxis()->SetTitle("abs((mean_yellow - mean_green) / sigma_green)");
    epsilon_test2->GetYaxis()->SetTitle("((mean_yellow - mean_green) - epsilon) / error_epsilon");
  
  }
  
  // Histograms for pull distributions

  double limit_p = 3.0, limit_m = -3.0;
  int bins = (limit_p - limit_m)/0.1;

  // Validation mode
  TH1D *pull_alpha_control=nullptr, *pull_nu_control=nullptr, *pull_epsilon_control=nullptr, *pull_epsilon_control1=nullptr, *pull_epsilon_control2=nullptr;
  
  if (mode_option.compare("validation") == 0){

    pull_alpha_control = new TH1D("pull_alpha_control", " ((sigma_yellow - sigma_green) - sigma_green*Delta_alpha) / error_alpha ", bins, limit_m, limit_p);
    pull_alpha_control->GetXaxis()->SetTitle("Pull");
    pull_alpha_control->GetYaxis()->SetTitle("Events");
    
    pull_nu_control = new TH1D("pull_nu_control", " ((integral_yellow - integral_green) - integral_green*Delta_nu) / error_nu ", bins, limit_m, limit_p);
    pull_nu_control->GetXaxis()->SetTitle("Pull");
    pull_nu_control->GetYaxis()->SetTitle("Events");
    
    pull_epsilon_control = new TH1D("pull_epsilon_control", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze none) ", bins, limit_m, limit_p);
    pull_epsilon_control->GetXaxis()->SetTitle("Pull");
    pull_epsilon_control->GetYaxis()->SetTitle("Events");
    
    pull_epsilon_control1 = new TH1D("pull_epsilon_control1", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze nu) ", bins, limit_m, limit_p);
    pull_epsilon_control1->GetXaxis()->SetTitle("Pull");
    pull_epsilon_control1->GetYaxis()->SetTitle("Events");
    
    pull_epsilon_control2 = new TH1D("pull_epsilon_control2", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze nu, alpha) ", bins, limit_m, limit_p);
    pull_epsilon_control2->GetXaxis()->SetTitle("Pull");
    pull_epsilon_control2->GetYaxis()->SetTitle("Events");
    
  }

  // Histograms for fit constituents for each k bin
  auto multi_hist_proj_mll_mc = mDh_mll_mc->Projection(4);
  auto multi_hist_proj_diff_mc = mDh_diff_mc->Projection(4);
  auto multi_hist_proj_jac_diff_squared_mll = mDh_jac_diff_squared_mll->Projection(4);
  auto multi_hist_proj_jac_diff_squared_mll_diff = mDh_jac_diff_squared_mll_diff->Projection(4);
  auto multi_hist_proj_jac_diff_mll_diff = mDh_jac_diff_mll_diff->Projection(4);
  auto multi_hist_proj_jac_diff_mll = mDh_jac_diff_mll->Projection(4);
  auto multi_hist_proj_mll_data = mDh_mll_data->Projection(4);

  TH1D *multi_hist_proj_diff_data=nullptr, *multi_hist_proj_diff_plus_offset=nullptr, *multi_hist_proj_diff_minus_offset=nullptr;
  if (mode_option.compare("validation") == 0){
    multi_hist_proj_diff_data = mDh_diff_data->Projection(4);
    multi_hist_proj_diff_plus_offset = mDh_diff_plus_offset->Projection(4);
    multi_hist_proj_diff_minus_offset = mDh_diff_minus_offset->Projection(4);
  }
  //--------------------------------------------------------------------------------------

  // Files to write results
  title = "mass_fits_control_histos_" + data_name_root + ".root";
  std::unique_ptr<TFile> f_control( TFile::Open(title.c_str(), "RECREATE") ); // histos inclusive in k bins
  title = "mass_fits_" + data_name_root + ".root";
  std::unique_ptr<TFile> f_fits( TFile::Open(title.c_str(), "RECREATE") ); // histos per k bin
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
  // Prepare counting events in mll_diff

  total_nevents = multi_hist_proj_diff_mc->Integral(1,nbinsmll_diff);
  std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";
  
  all_histos_count = 0.0;
  empty_histos_count = 0.0;
  hfrac = -1.0;  
  efrac = -1.0;
  remaining_nevents = 0.0;

  //--------------------------------------------------------------------------------------
  // Loop over eta+,pt+,eta-,pt- 

  for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
    mDh_mll_mc->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_mc->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_squared_mll->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_squared_mll_diff->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_mll_diff->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_jac_diff_mll->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_mll_data->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    if (mode_option.compare("validation") == 0){
      mDh_diff_data->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
      mDh_diff_plus_offset->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
      mDh_diff_minus_offset->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    }
    for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){   
      mDh_mll_mc->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_mc->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_squared_mll->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_squared_mll_diff->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_mll_diff->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_jac_diff_mll->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_mll_data->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      if (mode_option.compare("validation") == 0){
	mDh_diff_data->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
	mDh_diff_plus_offset->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
	mDh_diff_minus_offset->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      }
      for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	mDh_mll_mc->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_mc->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_diff_squared_mll->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_jac_diff_squared_mll_diff->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_diff_mll_diff->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_jac_diff_mll->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_mll_data->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	if (mode_option.compare("validation") == 0){
	  mDh_diff_data->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	  mDh_diff_plus_offset->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	  mDh_diff_minus_offset->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	}
	for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	  mDh_mll_mc->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_mc->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	    
	  all_histos_count++;
	    
	  //--------------------------------------------------------------------------
	    
	  // Require enough stats in diff_mc for sigma_MC fit
	  delete gROOT->FindObject( ("multi_data_histo_diff_" + mc_name_root + "_proj_4").c_str() );
	  multi_hist_proj_diff_mc = mDh_diff_mc->Projection(4);
	  nevents = multi_hist_proj_diff_mc->Integral(1,nbinsmll_diff); 
	    
	  if (nevents < 300.0){ // reject low stats
	    empty_histos_count++;
	  } else {

	    // Require most of mll_mc Gaussian to be in fit window
	    // TODO and require same from mll_data
	    delete gROOT->FindObject( ("multi_data_histo_mll_" + mc_name_root + "_proj_4").c_str() );
	    multi_hist_proj_mll_mc = mDh_mll_mc->Projection(4);
	    multi_hist_proj_mll_mc->GetXaxis()->SetTitle("mll [GeV]");
	    multi_hist_proj_mll_mc->GetYaxis()->SetTitle("Events");
	    fitresult = fitHisto(multi_hist_proj_mll_mc, 0, 8, 5);
	       
	    if (fitresult[4] < 0.75){ // reject small gaus integral 
	      empty_histos_count++;
	    } else {

	      name = stringify_name(pos_eta_bin-1, pos_pt_bin-1, neg_eta_bin-1, neg_pt_bin-1);
	      //std::cout<< name <<"\n";
              f_pass_reg << name << "\n";
	            
	      remaining_nevents += nevents;      
	      occupancy->SetBinError(occupancy->Fill(name.c_str(), nevents), 100);
	      gaus_integral->SetBinError(gaus_integral->Fill(name.c_str(), fitresult[4]), 0.01);
	            
	      //--------------------------------------------------------------------------
	      // Set range of the remaining histograms
	      mDh_jac_diff_squared_mll->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_squared_mll_diff->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_mll_diff->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_jac_diff_mll->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      mDh_mll_data->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      if (mode_option.compare("validation") == 0){
		mDh_diff_data->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
		mDh_diff_plus_offset->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
		mDh_diff_minus_offset->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	      }
	      //--------------------------------------------------------------------------
              // Diff histograms
              //--------------------------------------------------------------------------

	      //--------------------------------------------------------------------------
	      // diff_mc histogram

	      // Already projected diff_mc
	      multi_hist_proj_diff_mc->GetXaxis()->SetTitle("mll_diff [GeV]");
              multi_hist_proj_diff_mc->GetYaxis()->SetTitle("Events");
              multi_hist_proj_diff_mc->SetLineColor(kGreen);
	      // Fit diff_mc
              fitresult = fitHisto(multi_hist_proj_diff_mc, 1, 8, 2);
              if (fitresult[0] == -90.0 || fitresult[1] == -5.0){
		std::cout<<"WARNING bad diff_mc fit "<< name.c_str() <<" \n";
	      } 
	      mean_diff->SetBinError(mean_diff->Fill(name.c_str(), fitresult[0]), fitresult[2]);
	      sigma_diff->SetBinError(sigma_diff->Fill(name.c_str(), fitresult[1]), fitresult[3]);
	            
	      // Save for nu, epsilon, alpha fit
	      mean_mc = fitresult[0];
	      error_mean_mc = fitresult[2];
	      sigma_mc = fitresult[1];
	      error_sigma_mc = fitresult[3];
	      chi2_mc = fitresult[5];
	            
	      //-----------------------------------------------
	      // Prepare to draw mll_diff panel
	            
	      c1->cd(1); 
	      max_hist_mll_diff = -1.0;      
	            
	      if(max_hist_mll_diff < multi_hist_proj_diff_mc->GetBinContent(multi_hist_proj_diff_mc->GetMaximumBin())){
                max_hist_mll_diff = multi_hist_proj_diff_mc->GetBinContent(multi_hist_proj_diff_mc->GetMaximumBin());
              }
	            
	      //--------------------------------------------------------------------------
	      // diff_data histogram 
	            
	      if (mode_option.compare("validation") == 0){
		// Project diff_data histogram
		delete gROOT->FindObject( ("multi_data_histo_diff_" + data_name_root + "_proj_4").c_str() );
		multi_hist_proj_diff_data = mDh_diff_data->Projection(4);
		multi_hist_proj_diff_data->GetXaxis()->SetTitle("mll_diff [GeV]");
		multi_hist_proj_diff_data->GetYaxis()->SetTitle("Events");
		multi_hist_proj_diff_data->SetLineColor(kYellow);
		// Fit diff_data
		fitresult = fitHisto(multi_hist_proj_diff_data, 1, 5, 2);
		      
		// Save these for pull histograms 
		mean_diff_data = fitresult[0]; 
		error_mean_diff_data = fitresult[2];
		sigma_diff_data = fitresult[1];
		error_sigma_diff_data = fitresult[3];
		
		if(max_hist_mll_diff < multi_hist_proj_diff_data->GetBinContent(multi_hist_proj_diff_data->GetMaximumBin())){
		  max_hist_mll_diff = multi_hist_proj_diff_data->GetBinContent(multi_hist_proj_diff_data->GetMaximumBin());
		}
	      } else if (mode_option.compare("analysis") == 0){
		mean_diff_data = 0.0;
                error_mean_diff_data = 0.0;
                sigma_diff_data = 0.0;
                error_sigma_diff_data = 0.0;
	      }

	      //--------------------------------------------------------------------------
	      // Start drawing mll_diff 
	            
	      max_hist_mll_diff = max_hist_mll_diff * 1.3;
              multi_hist_proj_diff_mc->SetMaximum(max_hist_mll_diff);
	      multi_hist_proj_diff_mc->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      multi_hist_proj_diff_mc->Draw();
	      if (mode_option.compare("validation") == 0) multi_hist_proj_diff_data->Draw("SAME");
	            
	      // Legend

	      leg1->Clear();
	      leg_entry = "Region " + name + ": " + nevents + " events";
	      leg1->SetHeader(leg_entry.c_str(),"C");
                    
	      leg1->AddEntry(multi_hist_proj_diff_mc, "smeared-gen, #beta_pT=1, #alpha=1", "l");
	      if (mode_option.compare("validation") == 0) leg1->AddEntry(multi_hist_proj_diff_data, "smeared-gen, #beta_pT=0.999, #alpha=1, by weight", "l");
	      leg_entry = "green fit: #mu=" + to_string(mean_mc).substr(0, 5) + ", #sigma=" + to_string(sigma_mc).substr(0, 5);
	      leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	            
	      //--------------------------------------------------------------------------
	      // Mll histograms
	      //--------------------------------------------------------------------------

	      //--------------------------------------------------------------------------
	      // Prepare to draw mll panel

	      c1->cd();
	      c1->cd(2); 
	      max_hist_mll = -1.0;
	            
	      //--------------------------------------------------------------------------
	      // mll mc histogram

	      // Already projected and fitted mll mc
	      multi_hist_proj_mll_mc->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
              multi_hist_proj_mll_mc->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_mll_mc->GetYaxis()->SetTitle("Events");
              multi_hist_proj_mll_mc->SetLineColor(kGreen);
	            
	      if(max_hist_mll < multi_hist_proj_mll_mc->GetBinContent(multi_hist_proj_mll_mc->GetMaximumBin())){
                max_hist_mll = multi_hist_proj_mll_mc->GetBinContent(multi_hist_proj_mll_mc->GetMaximumBin());
              }

	      //--------------------------------------------------------------------------
	      // mll data
	            
	      // Project mll data
	      delete gROOT->FindObject( ("multi_data_histo_mll_" + data_name_root + "_proj_4").c_str() );
	      multi_hist_proj_mll_data = mDh_mll_data->Projection(4);
	      multi_hist_proj_mll_data->GetXaxis()->SetTitle("mll [GeV]");
              multi_hist_proj_mll_data->GetYaxis()->SetTitle("Events");
	      multi_hist_proj_mll_data->SetLineColor(kYellow);
	      // Fit mll data
	      //fitresult = fitHisto(multi_hist_proj_mll_data, 0, 5, 5);
              
	      if(max_hist_mll < multi_hist_proj_mll_data->GetBinContent(multi_hist_proj_mll_data->GetMaximumBin())){
                max_hist_mll = multi_hist_proj_mll_data->GetBinContent(multi_hist_proj_mll_data->GetMaximumBin());
	      }

	      //--------------------------------------------------------------------------
              // Start drawing mll
	            
	      max_hist_mll = max_hist_mll * 1.3;
	            
	      multi_hist_proj_mll_mc->Draw("");
	      multi_hist_proj_mll_data->Draw("SAME");

	      // Legend 
	      leg2->Clear();
	      leg_entry = "Region " + name + ": " + nevents + " events";
	      leg2->SetHeader(leg_entry.c_str(),"C"); 

	      leg2->AddEntry(multi_hist_proj_mll_mc, "smeared, #beta=1", "l");
	      leg2->AddEntry(multi_hist_proj_mll_data, "smeared, #beta=0.999", "l"); // insert value of beta here
	            
	      c1->cd();
	            
	      //--------------------------------------------------------------------------
	      // Start mass fit
	      //--------------------------------------------------------------------------

	      //--------------------------------------------------------------------------
	      // Fill vectors and variance for minimisation

	      // Find mll bins with defined jacobian
	            
	      filled_bins_mll=0;
	      vector<int> good_indices_mll;
	      for(int i=1; i<=nbinsmll; i++){
		//asking jacobian to be average over at least 10 events
		if ( multi_hist_proj_mll_mc->GetBinContent(i) >= 10 && multi_hist_proj_mll_data->GetBinContent(i) > 0 ){ //TODO refine the other criterium
		  good_indices_mll.push_back(i);
		  filled_bins_mll++;
		}
	      }
	      if (filled_bins_mll < 3 ){
		std::cout<< "WARNING not enough points in mll to fit nu, epsilon, alpha in "<< name.c_str() <<" \n";
		continue;
	      } 

	      // Declare mass vectors, jacobians and variance

              VectorXd h_data_minus_mc_mll_vector(filled_bins_mll);
	      Eigen::MatrixXd V_sqrt(filled_bins_mll, filled_bins_mll), V_inv_sqrt(filled_bins_mll, filled_bins_mll), J(filled_bins_mll, 3); //J.col(0) is nu, (1) is epsilon, (2) is alpha

              V_sqrt = MatrixXd::Zero(filled_bins_mll, filled_bins_mll);
              position_to_fill=0;
              for(int i=1; i<=nbinsmll; i++){
                if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
                  h_data_minus_mc_mll_vector(position_to_fill) = multi_hist_proj_mll_data->GetBinContent(i) - multi_hist_proj_mll_mc->GetBinContent(i);
                  J(position_to_fill,0) = multi_hist_proj_mll_mc->GetBinContent(i); //J_control.col(0) is nu
                  V_sqrt(position_to_fill,position_to_fill) = pow(pow(multi_hist_proj_mll_data->GetBinError(i),2) + pow(multi_hist_proj_mll_mc->GetBinError(i),2), 0.5); // sqrt(data_stat**2 + mc_stat**2)
                  position_to_fill++;
                }
              }
              if (position_to_fill != filled_bins_mll){ std::cout<<"problem counting vector size \n"; }
              // solve for V_inv_sqrt
              V_inv_sqrt = V_sqrt.completeOrthogonalDecomposition().solve(MatrixXd::Identity(filled_bins_mll,filled_bins_mll));

              //--------------------------------------------------------------------------
              // Define jacobian histograms

              // Jacobian histogram alpha mass
              delete gROOT->FindObject("jacobian");
              TH1D *jac = new TH1D("jacobian", "jacobian", nbinsmll, 75.0, 105.0);
              jac->GetXaxis()->SetTitle("mll [GeV]");
              jac->GetYaxis()->SetTitle("jacobian [GeV]");

              // Jacobian histogram epsilon mass
              delete gROOT->FindObject("jacobian_epsilon");
              TH1D *jac_epsilon = new TH1D("jacobian_epsilon", "jacobian epsilon", nbinsmll, 75.0, 105.0);
              jac_epsilon->GetXaxis()->SetTitle("mll [GeV]");
              jac_epsilon->GetYaxis()->SetTitle("jacobian [GeV]");

	      //--------------------------------------------------------------------------
              // Compute jacobians mass

              // Needed for jacobian alpha mass
              title = "multi_data_histo_jac_diff_squared_" + mc_name_root + "_mll_proj_4";
              delete gROOT->FindObject(title.c_str());
              multi_hist_proj_jac_diff_squared_mll = mDh_jac_diff_squared_mll->Projection(4);
              // Needed for jacobian epsilon mass
              title = "multi_data_histo_jac_diff_" + mc_name_root + "_mll_proj_4";
              delete gROOT->FindObject(title.c_str());
              multi_hist_proj_jac_diff_mll = mDh_jac_diff_mll->Projection(4);

              position_to_fill=0;
              for(int i=1; i<=nbinsmll; i++){
		// alpha mass
                diff_squared = multi_hist_proj_jac_diff_squared_mll->GetBinContent(i);
                error_diff_squared = multi_hist_proj_jac_diff_squared_mll->GetBinError(i);
                // epsilon mass
                jac_e_weight = multi_hist_proj_jac_diff_mll->GetBinContent(i);
                error_jac_e_weight = multi_hist_proj_jac_diff_mll->GetBinError(i);
		
                if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
		    
                  evts_in_bin = multi_hist_proj_mll_mc->GetBinContent(i);
                  // alpha mass
                  value =  diff_squared / (sigma_mc * sigma_mc) - evts_in_bin;
                  error = 1 / (sigma_mc*sigma_mc) * pow((4.0 * diff_squared*diff_squared * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_diff_squared*error_diff_squared) , 0.5);
                  J(position_to_fill,2) = value;

                  jac->SetBinContent(i, value);
                  jac->SetBinError(i, error);
                  jac_inclusive->Fill(mllbinranges[i-1], value);

                  // epsilon mass
                  value_epsilon = jac_e_weight / (sigma_mc * sigma_mc);
                  error_epsilon = 1 / (sigma_mc*sigma_mc) * pow((4 * jac_e_weight*jac_e_weight * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_jac_e_weight*error_jac_e_weight) , 0.5);
                  J(position_to_fill,1) = value_epsilon;

                  jac_epsilon->SetBinContent(i, value_epsilon);
                  jac_epsilon->SetBinError(i, error_epsilon);
                  jac_epsilon_inclusive->Fill(mllbinranges[i-1], value_epsilon);
		    
                  position_to_fill++;

                }
              }
	            
              if (position_to_fill != filled_bins_mll){ std::cout<<"PROBLEM: counting jac mass size \n"; }

	      //--------------------------------------------------------------------------
              // Write jacobian histograms

              f_fits->WriteObject(jac, ("jac_alpha" + name).c_str());
              f_fits->WriteObject(jac_epsilon, ("jac_epsilon" + name).c_str());

              //--------------------------------------------------------------------------
              // Solve for alpha, epsilon, nu
              //--------------------------------------------------------------------------

              // Solve for nu, epsilon, alpha mass
	      Eigen::MatrixXd A = V_inv_sqrt*J;
	      Eigen::MatrixXd b = V_inv_sqrt*h_data_minus_mc_mll_vector;
              // ATTENTION: n_e_a_vector contains nu, epsilon, alpha
	      Eigen::VectorXd n_e_a_vector = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

              // Error on nu, epsilon, alpha mass
	      Eigen::MatrixXd V_nu = (J.col(0).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	      Eigen::MatrixXd V_epsilon = (J.col(1).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
              fit_epsilon_error = pow(V_epsilon(0,0),0.5); // save for closure test
	      Eigen::MatrixXd V_alpha = (J.col(2).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(2)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));

              // Write nu, epsilon, alpha mass
              nu->SetBinError(nu->Fill(name.c_str(), n_e_a_vector(0)), pow(V_nu(0,0),0.5));
              epsilon->SetBinError(epsilon->Fill(name.c_str(), n_e_a_vector(1)), pow(V_epsilon(0,0),0.5));
              alpha->SetBinError(alpha->Fill(name.c_str(), n_e_a_vector(2)), pow(V_alpha(0,0),0.5));

	      //--------------------------------------------------------------------------
              // Add fitted parameters to mass fit panel

              c1->cd(2);
              leg_entry = "Fit #varepsilon=" + to_string(n_e_a_vector(1)).substr(0, 6) + "#pm" + to_string(pow(V_epsilon(0,0),0.5)).substr(0, 6);
              leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
              leg_entry = "Fit #Delta#nu=" + to_string(n_e_a_vector(0)).substr(0, 6) + "#pm" + to_string(pow(V_nu(0,0),0.5)).substr(0, 6);
              leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
              leg_entry = "Fit #Delta#alpha=" + to_string(n_e_a_vector(2)).substr(0, 6) + "#pm" + to_string(pow(V_alpha(0,0),0.5)).substr(0, 6);
              leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	            
	      //--------------------------------------------------------------------------
              // Apply fitted epsilon and nu and alpha correction
              
              // Define corrected_mll histogram
              delete gROOT->FindObject("corrected_mll");
              TH1D *corrected_mll = new TH1D("corrected_mll", "corrected_mll", nbinsmll, 75.0, 105.0);
              corrected_mll->SetTitle(("mll " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
              corrected_mll->GetXaxis()->SetTitle("mll [GeV]");
              corrected_mll->GetYaxis()->SetTitle("Events");

              position_to_fill = 0;
              for(int i=1; i<=nbinsmll; i++){
                if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
                  value_corrected = multi_hist_proj_mll_mc->GetBinContent(i) + n_e_a_vector(1)*J(position_to_fill,1) + n_e_a_vector(0)*J(position_to_fill,0) + n_e_a_vector(2)*J(position_to_fill,2);
                  error = pow(
                              pow(multi_hist_proj_mll_mc->GetBinError(i),2) +
                              pow(n_e_a_vector(1)*jac_epsilon->GetBinError(i),2) +
                              pow(J(position_to_fill,1),2)*V_epsilon(0,0) +
                              pow(n_e_a_vector(0)*multi_hist_proj_mll_mc->GetBinError(i),2) +
                              pow(J(position_to_fill,0),2)*V_nu(0,0) +
                              pow(n_e_a_vector(2)*jac->GetBinError(i),2) +
                              pow(J(position_to_fill,2),2)*V_alpha(0,0)
                              ,0.5);

                  corrected_mll->SetBinContent(i, value_corrected);
                  corrected_mll->SetBinError(i, error);

                  position_to_fill++;
                }
              }
              if (position_to_fill != filled_bins_mll){ std::cout<<"PROBLEM: counting corrected mll bins \n"; }

	      // Draw corrected_mll in mass fit panel
              c1->cd();
              c1->cd(2);
              corrected_mll->Draw("SAME");

              leg2->AddEntry(corrected_mll, "corrected mll", "l");
              leg2->Draw("");

	      if (mode_option.compare("validation") == 0){

		//--------------------------------------------------------------------------
		// Start diff fit
		//--------------------------------------------------------------------------
		
		//--------------------------------------------------------------------------
		// Find mll_diff bins with defined jacobian
		
		filled_bins_mll_diff=0;
		vector<int> good_indices_mll_diff;
		for(int i=1; i<=nbinsmll_diff; i++){
		  //asking jacobian to be average over at least 10 events
		  if( multi_hist_proj_diff_mc->GetBinContent(i) >= 10 && multi_hist_proj_diff_data->GetBinContent(i) > 0 ){ // TODO refine the other criterium
		    good_indices_mll_diff.push_back(i);
		    filled_bins_mll_diff++;
		  }
		}
		if (filled_bins_mll_diff < 3 ){
		  std::cout<<"WARNING not enough points in mll_diff to fit nu, epsilon, alpha in "<< name.c_str() <<"\n";
		  continue;
		} 
		
		// Declare mll_diff vectors, jacobians and variance
		VectorXd h_data_minus_mc_diff_vector(filled_bins_mll_diff);
		Eigen::MatrixXd V_sqrt_control(filled_bins_mll_diff, filled_bins_mll_diff), V_inv_sqrt_control(filled_bins_mll_diff, filled_bins_mll_diff), J_control(filled_bins_mll_diff, 3); //J_control.col(0) is nu, (1) is epsilon, (2) is alpha 
		
		V_sqrt_control = MatrixXd::Zero(filled_bins_mll_diff, filled_bins_mll_diff);
		position_to_fill = 0;
		integral_diff_data = 0.0;
		integral_mc = 0.0;
		
		for(int i=1; i<=nbinsmll_diff; i++){
		  if( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		    h_data_minus_mc_diff_vector(position_to_fill) = multi_hist_proj_diff_data->GetBinContent(i) - multi_hist_proj_diff_mc->GetBinContent(i);
		    J_control(position_to_fill,0) = multi_hist_proj_diff_mc->GetBinContent(i); //J_control.col(0) is nu
		    V_sqrt_control(position_to_fill,position_to_fill) = pow(pow(multi_hist_proj_diff_data->GetBinError(i),2) + pow(multi_hist_proj_diff_mc->GetBinError(i),2), 0.5); // sqrt(data_stat**2 + mc_stat**2)
		    integral_diff_data += multi_hist_proj_diff_data->GetBinContent(i);
		    integral_mc += multi_hist_proj_diff_mc->GetBinContent(i);
		        
		    position_to_fill++;
		  }
		}
		if (position_to_fill != filled_bins_mll_diff){ std::cout<<"problem counting vector size \n"; }
		// solve for V_inv_sqrt_control
		V_inv_sqrt_control = V_sqrt_control.completeOrthogonalDecomposition().solve(MatrixXd::Identity(filled_bins_mll_diff,filled_bins_mll_diff));
		
		//--------------------------------------------------------------------------
		// Define jacobian histograms
		
		// Jacobian histogram alpha diff 
		delete gROOT->FindObject("jacobian_control");
		TH1D *jac_control = new TH1D("jacobian_control", "jacobian", nbinsmll_diff, -5.0, 5.0);
		jac_control->GetXaxis()->SetTitle("mll_diff [GeV]");
		jac_control->GetYaxis()->SetTitle("jacobian [GeV]");
		
		// Jacobian histogram epsilon diff
		delete gROOT->FindObject("jacobian_epsilon_control");
		TH1D *jac_epsilon_control = new TH1D("jacobian_epsilon_control", "jacobian epsilon", nbinsmll_diff, -5.0, 5.0);
		jac_epsilon_control->GetXaxis()->SetTitle("mll_diff [GeV]");
		jac_epsilon_control->GetYaxis()->SetTitle("jacobian [GeV]");
		
		// Jacobian numerical histogram epsilon diff
		delete gROOT->FindObject("jacobian_epsilon_control_num");
		TH1D *jac_epsilon_control_num = new TH1D("jacobian_epsilon_control_num", "jacobian epsilon_num", nbinsmll_diff, -5.0, 5.0);
		jac_epsilon_control_num->GetXaxis()->SetTitle("mll_diff [GeV]");
		jac_epsilon_control_num->GetYaxis()->SetTitle("jacobian [GeV]");
		
		//--------------------------------------------------------------------------
		// Compute jacobians diff
		
		// Needed for jacobian alpha and epsilon diff
		title = "multi_data_histo_jac_diff_squared_" + mc_name_root + "_mll_diff_proj_4";
		delete gROOT->FindObject(title.c_str());
		multi_hist_proj_jac_diff_squared_mll_diff = mDh_jac_diff_squared_mll_diff->Projection(4);
		title = "multi_data_histo_jac_diff_" + mc_name_root + "_mll_diff_proj_4";
		delete gROOT->FindObject(title.c_str());
		multi_hist_proj_jac_diff_mll_diff = mDh_jac_diff_mll_diff->Projection(4);
		// Needed for jacobian numerical epsilon diff
		title = "multi_data_histo_diff_" + mc_name_root + "_plus_offset_proj_4";
		delete gROOT->FindObject(title.c_str());
		multi_hist_proj_diff_plus_offset = mDh_diff_plus_offset->Projection(4);
		title = "multi_data_histo_diff_" + mc_name_root + "_minus_offset_proj_4";
		delete gROOT->FindObject(title.c_str());
		multi_hist_proj_diff_minus_offset = mDh_diff_minus_offset->Projection(4);
		
		position_to_fill = 0;
		for(int i=1; i<=nbinsmll_diff; i++){
		  evts_in_bin = multi_hist_proj_diff_mc->GetBinContent(i);
		  // alpha and epsilon diff
		  diff_squared = multi_hist_proj_jac_diff_squared_mll_diff->GetBinContent(i);
		  error_diff_squared = multi_hist_proj_jac_diff_squared_mll_diff->GetBinError(i);
		  diff = multi_hist_proj_jac_diff_mll_diff->GetBinContent(i);
		  error_diff = multi_hist_proj_jac_diff_mll_diff->GetBinError(i);
		    
		  if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		    // alpha diff
		    value =  (diff_squared - 2.0*mean_mc*diff + evts_in_bin*mean_mc*mean_mc) / (sigma_mc * sigma_mc) - evts_in_bin;
		    error = 2.0 / (sigma_mc*sigma_mc) * pow( pow((diff_squared-2.0*mean_mc*diff + evts_in_bin*mean_mc*mean_mc)*error_sigma_mc/sigma_mc,2) + pow(mean_mc*error_diff,2) + pow(error_diff_squared,2)/4.0 + pow((evts_in_bin*mean_mc-diff)*error_mean_mc,2)  ,0.5);
		        
		    J_control(position_to_fill,2) = value; //J_control.col(2) is for alpha
		    jac_control->SetBinContent(i, value);
		    jac_control->SetBinError(i, error);
		        
		    // epsilon diff
		    // TODO check error and formula
		    value_epsilon = ( diff - evts_in_bin*mean_mc ) / pow(sigma_mc,2);
		    error_epsilon = pow( pow(evts_in_bin*error_mean_mc,2) + 4.0*pow(diff - evts_in_bin*mean_mc,2)*pow(error_sigma_mc,2)/pow(sigma_mc,2) + pow(error_diff,2), 0.5) / pow(sigma_mc,2);
		        
		    J_control(position_to_fill,1) = value_epsilon; //J_control.col(1) is for epsilon
		    jac_epsilon_control->SetBinContent(i, value_epsilon);
		    jac_epsilon_control->SetBinError(i, error_epsilon);
		      
		    // epsilon diff numerical
		    value_epsilon = (multi_hist_proj_diff_plus_offset->GetBinContent(i) - multi_hist_proj_diff_minus_offset->GetBinContent(i)) / (mll_offset * 2.0);
		    error_epsilon = pow((pow(multi_hist_proj_diff_plus_offset->GetBinError(i),2) + pow(multi_hist_proj_diff_minus_offset->GetBinError(i),2)) / pow(mll_offset*2.0, 2), 0.5);
		    //---------------------------------------------//
		    //---------------------------------------------//
		    // ATTENTION, for the numerical jacobian, OVERWRITE J_control
		    //J_control(position_to_fill,1) = value_epsilon;
		    //---------------------------------------------//
		    //---------------------------------------------//
		    jac_epsilon_control_num->SetBinContent(i, value_epsilon);
		    jac_epsilon_control_num->SetBinError(i, error_epsilon);

		    position_to_fill++;
		        
		  } 
		}
		
		if (position_to_fill != filled_bins_mll_diff){ std::cout<<"PROBLEM: counting jac diff size \n"; }

		//--------------------------------------------------------------------------
		// Write jacobian histograms
		
		f_fits->WriteObject(jac_control, ("jac_alpha_control" + name).c_str());
		f_fits->WriteObject(jac_epsilon_control, ("jac_epsilon_control" + name).c_str());
		f_fits->WriteObject(jac_epsilon_control_num, ("jac_epsilon_control_num" + name).c_str());
		
		//--------------------------------------------------------------------------
		
		// Solve for nu, epsilon, alpha diff 
		Eigen::MatrixXd A_control = V_inv_sqrt_control*J_control;
		Eigen::MatrixXd b_control = V_inv_sqrt_control*h_data_minus_mc_diff_vector;
		// ATTENTION: n_e_a_vector_control contains nu, epsilon, alpha
		Eigen::VectorXd n_e_a_vector_control = A_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_control);
		
		// Solve for epsilon, alpha diff only
		Eigen::MatrixXd A1_control = V_inv_sqrt_control*J_control.rightCols(2); 
		Eigen::MatrixXd b1_control = V_inv_sqrt_control*h_data_minus_mc_diff_vector;
		// ATTENTION: e_a_vector_control contains epsilon, alpha
		Eigen::VectorXd e_a_vector_control = A1_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b1_control);
		
		// Solve for epsilon diff only
		Eigen::MatrixXd A2_control = V_inv_sqrt_control*J_control.col(1);
		Eigen::MatrixXd b2_control = V_inv_sqrt_control*h_data_minus_mc_diff_vector;
		// ATTENTION: e_vector_control contains epsilon
		Eigen::VectorXd e_vector_control = A2_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b2_control);
		      
		// Error on epsilon diff 
		Eigen::MatrixXd V_epsilon_control(1,1); 
		V_epsilon_control = (J_control.col(1).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
		epsilon_control->SetBinError(epsilon_control->Fill(name.c_str(), n_e_a_vector_control(1)), pow(V_epsilon_control(0,0),0.5));
		// Error on nu diff
		Eigen::MatrixXd V_nu_control(1,1); 
		V_nu_control = (J_control.col(0).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
		nu_control->SetBinError(nu_control->Fill(name.c_str(), n_e_a_vector_control(0)), pow(V_nu_control(0,0),0.5)); 
		// Error on alpha diff 
		Eigen::MatrixXd V_alpha_control(1,1);
		V_alpha_control = (J_control.col(2).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(2)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
		alpha_control->SetBinError(alpha_control->Fill(name.c_str(), n_e_a_vector_control(2)), pow(V_alpha_control(0,0),0.5));
		
		//--------------------------------------------------------------------------
		// Add fitted parameters to diff fit panel 
		
		c1->cd();
		c1->cd(1);
		
		leg_entry = "Fit #varepsilon=" + to_string(n_e_a_vector_control(1)).substr(0, 6) + "#pm" + to_string(pow(V_epsilon_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(1) is epsilon
		leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
		leg_entry = "Fit #Delta#nu=" + to_string(n_e_a_vector_control(0)).substr(0, 6) + "#pm" + to_string(pow(V_nu_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(0) is nu 
		leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
		leg_entry = "Fit #Delta#alpha=" + to_string(n_e_a_vector_control(2)).substr(0, 6) + "#pm" + to_string(pow(V_alpha_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(2) is alpha
		leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
		
		//--------------------------------------------------------------------------
		// Define corrected_diff histogram
		delete gROOT->FindObject("corrected_diff");
		TH1D *corrected_diff = new TH1D("corrected_diff", "corrected_diff", nbinsmll_diff, -5.0, 5.0);
		corrected_diff->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
		corrected_diff->GetXaxis()->SetTitle("mll_diff [GeV]");
		corrected_diff->GetYaxis()->SetTitle("Events");
		
		position_to_fill = 0;
		for(int i=1; i<=nbinsmll_diff; i++){
		  if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
		    value_corrected = multi_hist_proj_diff_mc->GetBinContent(i) + n_e_a_vector_control(1)*J_control(position_to_fill,1) + n_e_a_vector_control(0)*J_control(position_to_fill,0) + n_e_a_vector_control(2)*J_control(position_to_fill,2);
		    error = pow(
				pow(multi_hist_proj_diff_mc->GetBinError(i),2) +
				pow(n_e_a_vector_control(1)*jac_epsilon_control->GetBinError(i),2) +
				pow(J_control(position_to_fill,1),2)*V_epsilon_control(0,0) +
				pow(n_e_a_vector_control(0)*multi_hist_proj_diff_mc->GetBinError(i),2) + 
				pow(J_control(position_to_fill,0),2)*V_nu_control(0,0) +
				pow(n_e_a_vector_control(2)*jac_control->GetBinError(i),2) + 
				pow(J_control(position_to_fill,2),2)*V_alpha_control(0,0)
				,0.5);
		        
		    corrected_diff->SetBinContent(i, value_corrected);
		    corrected_diff->SetBinError(i, error);
		        
		    position_to_fill++;
		  }
		}
		if (position_to_fill != filled_bins_mll_diff){ std::cout<<"PROBLEM: counting corrected_diff bins \n"; }
		
		// Draw corrected_diff in diff fit panel
		c1->cd();
		c1->cd(1);
              
		fitresult = fitHisto(corrected_diff, 1, 1, 2);
		      
		// Save for pull distribution epsilon_control
		mean_corrected_diff = fitresult[0];
		sigma_corrected_diff = fitresult[1];
		error_mean_corrected_diff = fitresult[2];
		error_sigma_corrected_diff = fitresult[3];
		corrected_diff->SetLineColor(kBlack);
		corrected_diff->Draw("SAME");
		
		// Draw pull distributions
		
		// nu diff
		pull_nu_control->Fill( (integral_diff_data - integral_mc - integral_mc*n_e_a_vector_control(0)) / integral_mc / pow(V_nu_control(0,0),0.5) );
		
		// epsilon diff
		pull_epsilon_control->Fill( (mean_diff_data - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) );
		pull_epsilon_control1->Fill( (mean_diff_data - mean_mc - e_a_vector_control(0))/pow(V_epsilon_control(0,0),0.5) );
		pull_epsilon_control2->Fill( (mean_diff_data - mean_mc - e_vector_control(0))/pow(V_epsilon_control(0,0),0.5) );
		
		// epsilon test
		epsilon_test1->SetBinError(epsilon_test1->Fill(name.c_str(), (mean_diff_data - mean_mc) / sigma_mc), 0.001);
		if ( abs((mean_diff_data - mean_mc) / sigma_mc) > 0.05 ){
		  epsilon_test2->Fill(abs((mean_diff_data - mean_mc) / sigma_mc), (mean_diff_data - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5));
		}
		
		// alpha diff
		pull_alpha_control->Fill( (sigma_diff_data - sigma_mc - sigma_mc*n_e_a_vector_control(2))/pow(V_alpha_control(0,0),0.5) );
		
		//--------------------------------------------------------------------------
		// Add to diff panel
		leg1->AddEntry(corrected_diff, "corrected diff", "l");
		leg_entry = "Pull #varepsilon=" + to_string( (mean_diff_data - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) ).substr(0, 6);
		leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	      }

	      // Finish drawing diff panel
	      c1->cd();
	      c1->cd(1);
              leg1->Draw("");
	            
              c1->cd();
              // Write plots per k bin
              f_fits->WriteObject(c1, name.c_str());
	            
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

  gaus_integral->SetStats(0);
  gaus_integral->LabelsDeflate();
  f_control->WriteObject(gaus_integral, "gaus_integral");

  occupancy->SetStats(0);
  occupancy->LabelsDeflate();
  f_control->WriteObject(occupancy, "bin_occupancy");

  mean_diff->SetStats(0);
  mean_diff->LabelsDeflate();
  f_control->WriteObject(mean_diff, "mean_diff");

  sigma_diff->SetStats(0);
  sigma_diff->LabelsDeflate();
  f_control->WriteObject(sigma_diff, "sigma_diff");

  f_control->WriteObject(alpha, "alpha");
  f_control->WriteObject(jac_inclusive, "jacobian_inclusive");

  f_control->WriteObject(epsilon, "epsilon");
  f_control->WriteObject(jac_epsilon_inclusive, "jacobian_epsilon_inclusive");

  f_control->WriteObject(nu, "nu");

  // Validation plots
  if (mode_option.compare("validation") == 0){
    f_control->WriteObject(epsilon_control, "epsilon_control");
    fitresult = fitHisto(pull_epsilon_control, 1, 1, 5);
    std::cout<<"CHECKPOINT: "<<"Pull epsilon control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
    f_control->WriteObject(pull_epsilon_control, "pull_epsilon_control");
    
    f_control->WriteObject(epsilon_test1, "epsilon_test1");
    
    f_control->WriteObject(alpha_control, "alpha_control");
    fitresult = fitHisto(pull_alpha_control, 1, 1, 5);
    std::cout<<"CHECKPOINT: "<<"Pull alpha control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
    f_control->WriteObject(pull_alpha_control, "pull_alpha_control");
    
    f_control->WriteObject(nu_control, "nu_control");
    fitresult = fitHisto(pull_nu_control, 1, 1, 5);
    std::cout<<"CHECKPOINT: "<<"Pull nu control distribution fitted mean: "<<fitresult[0]<<" +/- "<<fitresult[2]<<" and sigma "<<fitresult[1]<<" +/- "<<fitresult[3]<<"\n";
    f_control->WriteObject(pull_nu_control, "pull_nu_control");

    //--------------------------------------------------------------------
    // Superimposed histograms

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
    fitresult = fitHisto(pull_epsilon_control, 1, 1, 5);
    leg_entry = "freeze none, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
    leg4->AddEntry(pull_epsilon_control, leg_entry.c_str(), "l");
    pull_epsilon_control->SetLineColor(kBlack);
    pull_epsilon_control->Draw("E");
    /*
      fitresult = fitHisto(pull_epsilon_control1, 1, 2, 5);
      leg_entry = "freeze nu, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
      leg4->AddEntry(pull_epsilon_control1, leg_entry.c_str(), "l");
      pull_epsilon_control1->SetLineColor(kRed);
      pull_epsilon_control1->Draw("SAME");
      
      fitresult = fitHisto(pull_epsilon_control2, 1, 3, 5);
      leg_entry = "freeze nu, alpha, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
      leg4->AddEntry(pull_epsilon_control2, leg_entry.c_str(), "l");
      pull_epsilon_control2->SetLineColor(kGreen);
      pull_epsilon_control2->Draw("SAME");
      
    */
    leg4->Draw("");
  
    f_control->WriteObject(c4, "pull_epsilon_many");
    
    //-----------------------------------------------------------
    // Other histograms
    
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    epsilon_test2->SetStats(0);
    epsilon_test2->Draw("COLZ");
    f_control->WriteObject(c3, "epsilon_test2");
  }

  return 0; 
}
