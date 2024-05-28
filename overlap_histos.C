#include "TFile.h"
#include "TCanvas.h"

using namespace ROOT;

int overlap_histos(){

  string filename1("mass_fits_control_histos_smear_beta_val_unique_pdf.root"), filename2("mass_fits_control_histos_smear_beta_val_additive_k_bin_dependent_pdf.root"), hist_name1("beta"), hist_name2("beta");

  std::unique_ptr<TFile> file1( TFile::Open( filename1.c_str() ) );
  std::unique_ptr<TFile> file2( TFile::Open( filename2.c_str() ) );

  std::unique_ptr<TH1D> hist1(file1->Get<TH1D>(hist_name1.c_str()));
  std::unique_ptr<TH1D> hist2(file2->Get<TH1D>(hist_name2.c_str()));

  std::unique_ptr<TFile> f_out(TFile::Open("superimposed.root", "RECREATE")); 

  // Canvas
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  
  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.015);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  double max_hist, min_hist;
  
  max_hist = hist1->GetBinContent(hist1->GetMaximumBin());
  min_hist = hist1->GetBinContent(hist1->GetMinimumBin());
  
  if(max_hist < hist2->GetBinContent(hist2->GetMaximumBin())){
    max_hist = hist2->GetBinContent(hist2->GetMaximumBin());
  }

  if(min_hist < hist2->GetBinContent(hist2->GetMinimumBin())){
    min_hist = hist2->GetBinContent(hist2->GetMinimumBin());
  }

  hist1->SetMaximum(max_hist);
  hist1->SetMinimum(min_hist);
  
  hist1->SetLineColor(kRed);
  hist1->SetMarkerColor(kRed);
  hist2->SetLineColor(kBlue);
  hist2->SetMarkerColor(kBlue);
  
  hist1->Draw();
  hist2->Draw("SAME");

  string leg_entry = "red: " + filename1;
  leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");

  leg_entry = "blue: " + filename2;
  leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
  leg1->Draw("SAME");

  f_out->WriteObject(c1, "superimposed");
  
  return 0;
}
