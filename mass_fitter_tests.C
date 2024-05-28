#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

using namespace ROOT;

// Function for Gaussian fit and integral 
vector<double> fitHisto(ROOT::RDF::RResultPtr<TH1D> histogram, int draw_option, int color, int nsigma){

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
      integral = gaussianNewFunc->Integral(-5.0, 5.0);
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

int mass_fitter_tests(){

  ROOT::EnableImplicitMT(128);

  RDataFrame df("control_tree", "InOutputFiles/control_tree.root");

  auto d1 = df.Define("pull_delta_alpha","return alpha/error_alpha;");
  
  auto pull_alpha_hist = d1.Histo1D({"pull_alpha_hist", "#Delta alpha - 0 / error", 20, -5.0, 5.0},"pull_delta_alpha");
    
  vector<double> fitresult;
  fitresult = fitHisto(pull_alpha_hist, 1, 2, 5);
  
  std::unique_ptr<TFile> f_out( TFile::Open("InOutputFiles/mass_fits_control_histos_more.root", "RECREATE") ); 
  pull_alpha_hist->Write();
  f_out->Close();

  return 0;
}
