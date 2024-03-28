#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"

using namespace ROOT;
using namespace ROOT::VecOps;

int frame(){

  ROOT::EnableImplicitMT(128);

  TChain chain("track_tree");
  chain.Add("/scratchnvme/alexe/MomentumBiases/track_params_no_cut.root");
  chain.Add("/scratchnvme/alexe/MomentumBiases/track_params_no_cut_2kmore.root");

  //RDataFrame df("track_tree", {"/scratchnvme/alexe/MomentumBiases/track_params_no_cut.root"});
  RDataFrame df(chain);
  RDataFrame df_cut("track_tree", "/scratchnvme/alexe/MomentumBiases/track_params_M_no_cut.root");

  auto red_chi2_edm = df.Histo2D({"red_chi2_edm", "edm vs red chi2", 50, 0.9630, 1.046, 50, 4.361e-13, 0.001951},"red_chi2","edm");
  auto edm_low = df.Histo1D({"edm_low", "edm bin", 30, 5.072e-4, 5.462e-4},"edm");
  
  auto d_edm_high = df.Filter("edm > 5.345e-4 ");
  auto d_edm_low = df.Filter("edm <= 5.345e-6 "); // 90% toys pass

  auto A12_edm_high = d_edm_high.Histo1D({"A12_edm_high", "pull A12: edm > 5.345e-4", 50, -2.998, 3.139},"pull_A12");
  auto A12_edm_low = d_edm_low.Histo1D({"A12_edm_low", "pull A12: edm <= 5.345e-6", 50, -2.998, 3.139},"pull_A12");
  auto M12_edm_high = d_edm_high.Histo1D({"M12_edm_high", "pull M12: edm > 5.345e-4", 50, -4.532, 3.725},"pull_M12");
  auto M12_edm_low = d_edm_low.Histo1D({"M12_edm_low", "pull M12: edm <= 5.345e-6", 50, -4.532, 3.725},"pull_M12");

  auto red_chi2 = df.Histo1D({"red_chi2", "chi2/ndf ", 50, 0.9630, 1.046},"red_chi2");
  auto red_chi2_low = d_edm_low.Histo1D({"red_chi2_low", "chi2/ndf edm <= 5.345e-6", 50, 0.9630, 1.046},"red_chi2");
  auto red_chi2_cut = df_cut.Histo1D({"red_chi2", "chi2/ndf ", 50, 0.9630, 1.046},"red_chi2");

  TH1D *param_mean_A_edm_low = new TH1D("param_mean_A_edm_low", "param_pull_mean_edm_low A", 3, 0.0, 3.0);
  param_mean_A_edm_low->SetCanExtend(TH1::kAllAxes);
  TH1D *param_width_A_edm_low = new TH1D("param_width_A_edm_low", "param_pull_width_edm_low A", 3, 0.0, 3.0);
  param_width_A_edm_low->SetCanExtend(TH1::kAllAxes);

  TH1D *param_mean_A = new TH1D("param_mean_A", "param_pull_mean A", 3, 0.0, 3.0);
  param_mean_A->SetCanExtend(TH1::kAllAxes);
  TH1D *param_width_A = new TH1D("param_width_A", "param_pull_width A", 3, 0.0, 3.0);
  param_width_A->SetCanExtend(TH1::kAllAxes);

  TH1D *param_mean_A_cut = new TH1D("param_mean_A_cut", "param_pull_mean_cut A", 3, 0.0, 3.0);
  param_mean_A_cut->SetCanExtend(TH1::kAllAxes);
  TH1D *param_width_A_cut = new TH1D("param_width_A_cut", "param_pull_width_cut A", 3, 0.0, 3.0);
  param_width_A_cut->SetCanExtend(TH1::kAllAxes);

  TH1D *param_mean_M_edm_low = new TH1D("param_mean_M_edm_low", "param_pull_mean_edm_low M", 3, 0.0, 3.0);
  param_mean_M_edm_low->SetCanExtend(TH1::kAllAxes);
  TH1D *param_width_M_edm_low = new TH1D("param_width_M_edm_low", "param_pull_width_edm_low M", 3, 0.0, 3.0);
  param_width_M_edm_low->SetCanExtend(TH1::kAllAxes);

  TH1D *param_mean_M = new TH1D("param_mean_M", "param_pull_mean M", 3, 0.0, 3.0);
  param_mean_M->SetCanExtend(TH1::kAllAxes);
  TH1D *param_width_M = new TH1D("param_width_M", "param_pull_width M", 3, 0.0, 3.0);
  param_width_M->SetCanExtend(TH1::kAllAxes);

  TH1D *param_mean_M_cut = new TH1D("param_mean_M_cut", "param_pull_mean_cut M", 3, 0.0, 3.0);
  param_mean_M_cut->SetCanExtend(TH1::kAllAxes);
  TH1D *param_width_M_cut = new TH1D("param_width_M_cut", "param_pull_width_cut M", 3, 0.0, 3.0);
  param_width_M_cut->SetCanExtend(TH1::kAllAxes);

  auto minVal = d_edm_low.Min<double>("pull_A0");
  auto maxVal = d_edm_low.Max<double>("pull_A0");
  
  auto temp = df.Histo1D({"pull_A0", "edm > 5.345e-4", 1, 0, 0.01},"pull_A0"); 

  string name;
  for(int i=0; i<=23; i++){ //TODO loop upper bound input by hand
    
    name = Form("pull_A%d",i);
    delete gROOT->FindObject(name.c_str());
    minVal = d_edm_low.Min<double>(name.c_str());
    maxVal = d_edm_low.Max<double>(name.c_str());
    temp = d_edm_low.Histo1D({name.c_str(), "edm > 5.345e-4", 50, *minVal, *maxVal},name.c_str());
    param_mean_A_edm_low->SetBinError(param_mean_A_edm_low->Fill(name.c_str(), temp->GetMean()), temp->GetMeanError());
    param_width_A_edm_low->SetBinError(param_width_A_edm_low->Fill(name.c_str(),temp->GetStdDev()), temp->GetStdDevError());

    delete gROOT->FindObject(name.c_str());
    minVal = df.Min<double>(name.c_str());
    maxVal = df.Max<double>(name.c_str());
    temp = df.Histo1D({name.c_str(), "edm > 5.345e-4", 50, *minVal, *maxVal},name.c_str());
    param_mean_A->SetBinError(param_mean_A->Fill(name.c_str(), temp->GetMean()), temp->GetMeanError());
    param_width_A->SetBinError(param_width_A->Fill(name.c_str(),temp->GetStdDev()), temp->GetStdDevError());
    
    delete gROOT->FindObject(name.c_str());
    minVal = df_cut.Min<double>(name.c_str());
    maxVal = df_cut.Max<double>(name.c_str());
    temp = df_cut.Histo1D({name.c_str(), "edm > 5.345e-4", 50, *minVal, *maxVal},name.c_str());
    param_mean_A_cut->SetBinError(param_mean_A_cut->Fill(name.c_str(), temp->GetMean()), temp->GetMeanError());
    param_width_A_cut->SetBinError(param_width_A_cut->Fill(name.c_str(),temp->GetStdDev()), temp->GetStdDevError());

    name = Form("pull_M%d",i);
    delete gROOT->FindObject(name.c_str());
    minVal = d_edm_low.Min<double>(name.c_str());
    maxVal = d_edm_low.Max<double>(name.c_str());
    temp = d_edm_low.Histo1D({name.c_str(), "edm > 5.345e-4", 50, *minVal, *maxVal},name.c_str());
    param_mean_M_edm_low->SetBinError(param_mean_M_edm_low->Fill(name.c_str(), temp->GetMean()), temp->GetMeanError());
    param_width_M_edm_low->SetBinError(param_width_M_edm_low->Fill(name.c_str(),temp->GetStdDev()), temp->GetStdDevError());

    delete gROOT->FindObject(name.c_str());
    minVal = df.Min<double>(name.c_str());
    maxVal = df.Max<double>(name.c_str());
    temp = df.Histo1D({name.c_str(), "edm > 5.345e-4", 50, *minVal, *maxVal},name.c_str());
    param_mean_M->SetBinError(param_mean_M->Fill(name.c_str(), temp->GetMean()), temp->GetMeanError());
    param_width_M->SetBinError(param_width_M->Fill(name.c_str(),temp->GetStdDev()), temp->GetStdDevError());

    delete gROOT->FindObject(name.c_str());
    minVal = df_cut.Min<double>(name.c_str());
    maxVal = df_cut.Max<double>(name.c_str());
    temp = df_cut.Histo1D({name.c_str(), "edm > 5.345e-4", 50, *minVal, *maxVal},name.c_str());
    param_mean_M_cut->SetBinError(param_mean_M_cut->Fill(name.c_str(), temp->GetMean()), temp->GetMeanError());
    param_width_M_cut->SetBinError(param_width_M_cut->Fill(name.c_str(),temp->GetStdDev()), temp->GetStdDevError());
  }

  TFile f("constants_fitter_test.root","recreate");

  TCanvas *c = new TCanvas("c","c",800,600);
  red_chi2_edm->GetXaxis()->SetTitle("#chi^{2}/ndf");
  red_chi2_edm->GetYaxis()->SetTitle("EDM");
  red_chi2_edm->Draw("COLZ");
  f.WriteObject(c, "red_chi2_edm");
  f.WriteObject(edm_low.GetPtr(), "edm_bin");
  f.WriteObject(A12_edm_high.GetPtr(), "A12_edm_high");
  f.WriteObject(A12_edm_low.GetPtr(), "A12_edm_low");
  f.WriteObject(M12_edm_high.GetPtr(), "M12_edm_high");
  f.WriteObject(M12_edm_low.GetPtr(), "M12_edm_low");

  TCanvas *c0 = new TCanvas("c0","c0",800,600);

  auto leg0 = new TLegend(0.68, 0.78, 0.90, 0.90);

  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->SetTextSize(0.025);
  leg0->SetFillColor(10);
  leg0->SetNColumns(1);
  leg0->SetHeader("");

  red_chi2_cut->SetLineColor(kGreen);
  red_chi2_cut->GetXaxis()->SetTitle("#chi^{2}/ndf");
  red_chi2_cut->Draw();
  leg0->AddEntry("red_chi2_cut", "GREEN: no cut");
  
  red_chi2->SetLineColor(kBlue);
  red_chi2->Draw("SAME");
  leg0->AddEntry("red_chi2", "BLUE: no edm cut");

  red_chi2_low->SetLineColor(kRed);
  red_chi2_low->Draw("SAME");
  leg0->AddEntry("red_chi2_low", "RED: select edm <= 5.345e-6");

  leg0->Draw("SAME");
  f.WriteObject(c0, "red_chi2");

  TCanvas *c2 = new TCanvas("c2","c2",800,600);

  auto leg = new TLegend(0.68, 0.78, 0.90, 0.90);
    
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->SetFillColor(10);
  leg->SetNColumns(1);
  leg->SetHeader("");
  /*
  param_mean_A_cut->Fit("pol0");
  param_mean_A_cut->GetFunction("pol0")->SetLineColor(kGreen);
  param_mean_A->Fit("pol0");
  param_mean_A->GetFunction("pol0")->SetLineColor(kBlue);
  param_mean_A_edm_low->Fit("pol0");
  param_mean_A_edm_low->GetFunction("pol0")->SetLineColor(kRed);
  */
  param_mean_A_cut->SetStats(0);
  param_mean_A_cut->SetLineColor(kGreen);
  leg->AddEntry(param_mean_A_cut, "cut edm <= 5.345e-4", "l");
  param_mean_A_cut->Draw();

  param_mean_A->SetStats(0);
  param_mean_A->SetLineColor(kBlue);
  leg->AddEntry(param_mean_A, "no edm cut", "l");
  param_mean_A->Draw("SAME");

  param_mean_A_edm_low->SetStats(0);
  param_mean_A_edm_low->SetLineColor(kRed);
  leg->AddEntry(param_mean_A_edm_low, "select edm <= 5.345e-4", "l");
  param_mean_A_edm_low->GetXaxis()->SetTitle("Parameter name");
  param_mean_A_edm_low->GetYaxis()->SetTitle("#mu");
  param_mean_A_edm_low->Draw("SAME");
  
  leg->Draw("");
  f.WriteObject(c2, "param_mean_A");

  TCanvas *c3 = new TCanvas("c3","c3",800,600);

  auto leg1 = new TLegend(0.68, 0.78, 0.90, 0.90);

  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.025);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  /*
  param_width_A_cut->Fit("pol0");
  param_width_A_cut->GetFunction("pol0")->SetLineColor(kGreen);
  param_width_A->Fit("pol0");
  param_width_A->GetFunction("pol0")->SetLineColor(kBlue);
  param_width_A_edm_low->Fit("pol0");
  param_width_A_edm_low->GetFunction("pol0")->SetLineColor(kRed);
  */
  param_width_A_cut->SetStats(0);
  param_width_A_cut->SetLineColor(kGreen);
  leg1->AddEntry(param_width_A_cut, "cut edm <= 5.345e-4", "l");
  param_width_A_cut->SetMaximum(1.1);
  param_width_A_cut->Draw();

  param_width_A->SetStats(0);
  param_width_A->SetLineColor(kBlue);
  leg1->AddEntry(param_width_A, "no edm cut", "l");
  param_width_A->Draw("SAME");

  param_width_A_edm_low->SetStats(0);
  param_width_A_edm_low->SetLineColor(kRed);
  leg1->AddEntry(param_width_A_edm_low, "select edm <= 5.345e-4", "l");
  param_width_A_edm_low->GetXaxis()->SetTitle("Parameter name");
  param_width_A_edm_low->GetYaxis()->SetTitle("#sigma");
  param_width_A_edm_low->Draw("SAME");

  leg1->Draw("SAME");
  f.WriteObject(c3, "param_width_A");

  TCanvas *c4 = new TCanvas("c4","c4",800,600);

  auto leg2 = new TLegend(0.68, 0.78, 0.90, 0.90);

  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.025);
  leg2->SetFillColor(10);
  leg2->SetNColumns(1);
  leg2->SetHeader("");

  param_mean_M_cut->Fit("pol0");
  param_mean_M->Fit("pol0");
  param_mean_M_edm_low->Fit("pol0");

  param_mean_M_cut->SetStats(0);
  param_mean_M_cut->SetLineColor(kGreen);
  //param_mean_M_cut->Fit("pol0");
  param_mean_M_cut->GetFunction("pol0")->SetLineColor(kGreen);
  leg2->AddEntry(param_mean_M_cut, "no cut ", "l");
  param_mean_M_cut->SetMaximum(0.5);
  //param_mean_M_cut->Draw();

  param_mean_M->SetStats(0);
  param_mean_M->SetLineColor(kBlue);
  //param_mean_M->Fit("pol0");
  param_mean_M->GetFunction("pol0")->SetLineColor(kBlue);
  leg2->AddEntry(param_mean_M, "no edm cut", "l");
  param_mean_M->Draw("");
  
  param_mean_M_edm_low->SetStats(0);
  param_mean_M_edm_low->SetLineColor(kRed);
  //param_mean_M_edm_low->Fit("pol0");
  param_mean_M_edm_low->GetFunction("pol0")->SetLineColor(kRed);
  leg2->AddEntry(param_mean_M_edm_low, "select edm <= 5.345e-6", "l");
  param_mean_M_edm_low->GetXaxis()->SetTitle("Parameter name");
  param_mean_M_edm_low->GetYaxis()->SetTitle("#mu");
  param_mean_M_edm_low->Draw("SAME");
  
  leg2->Draw("SAME");
  f.WriteObject(c4, "param_mean_M");

  TCanvas *c5 = new TCanvas("c5","c5",800,600);

  auto leg3 = new TLegend(0.68, 0.78, 0.90, 0.90);

  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.025);
  leg3->SetFillColor(10);
  leg3->SetNColumns(1);
  leg3->SetHeader("");

  param_width_M_cut->Fit("pol0");
  param_width_M->Fit("pol0");
  param_width_M_edm_low->Fit("pol0");

  param_width_M_cut->SetStats(0);
  param_width_M_cut->SetLineColor(kGreen);
  //param_width_M_cut->Fit("pol0");
  param_width_M_cut->GetFunction("pol0")->SetLineColor(kGreen);
  leg3->AddEntry(param_width_M_cut, "no cut ", "l");
  param_width_M_cut->SetMaximum(1.5);
  //param_width_M_cut->Draw();
  
  param_width_M->SetStats(0);
  param_width_M->SetLineColor(kBlue);
  //param_width_M->Fit("pol0");
  param_width_M->GetFunction("pol0")->SetLineColor(kBlue);
  leg3->AddEntry(param_width_M, "no edm cut", "l");
  param_width_M->Draw("");
  
  param_width_M_edm_low->SetStats(0);
  param_width_M_edm_low->SetLineColor(kRed);
  //param_width_M_edm_low->Fit("pol0");
  param_width_M_edm_low->GetFunction("pol0")->SetLineColor(kRed);
  leg3->AddEntry(param_width_M_edm_low, "select edm <= 5.345e-6", "l");
  param_width_M_edm_low->GetXaxis()->SetTitle("Parameter name");
  param_width_M_edm_low->GetYaxis()->SetTitle("#sigma");
  param_width_M_edm_low->Draw("SAME");
  
  leg3->Draw("SAME");
  f.WriteObject(c5, "param_width_M");

  f.Close();

  return 0;

}
