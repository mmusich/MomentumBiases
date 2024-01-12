// Standalone code to fit for sagitta bias corrections

#include "TFile.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TVectorT.h"
//#include <TMatrixD.h>
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
  string txt1 = "eta+ in [", txt2="] pt+ in [", txt3="] eta- in [", txt4="] pt- in [", txt5 = "]";
  return txt1+to_string(e[a-1]).substr(0,4)+comma+to_string(e[a]).substr(0,4)+txt2+to_string(p[b-1]).substr(0,4)+comma+to_string(p[b]).substr(0,4)+txt3+to_string(e[c-1]).substr(0,4)+comma+to_string(e[c]).substr(0,4)+txt4+to_string(p[d-1]).substr(0,4)+comma+to_string(p[d]).substr(0,4)+txt5;
}

vector<double> fitHisto(TH1* histogram, int color){

  vector<double> fitresult;

  double mean = histogram->GetMean();
  double sigma = histogram->GetRMS();
  double mean_err, sigma_err, integral;

  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", mean - 1.5 * sigma, mean + 1.5 * sigma); 
  if(0 == histogram->Fit(gaussianFunc, "QNR")){
    mean = gaussianFunc->GetParameter(1);
    sigma = gaussianFunc->GetParameter(2);

    // second fit: two sigma of first fit around mean of first fit
    gaussianFunc->SetRange(mean - 2 * sigma, mean + 2 * sigma);
    if (0 == histogram->Fit(gaussianFunc, "Q0R")) {
      if (histogram->GetFunction(gaussianFunc->GetName())) { // Take care that it is later on drawn:
	histogram->GetFunction(gaussianFunc->GetName())->SetLineColor(color);
        histogram->GetFunction(gaussianFunc->GetName())->ResetBit(TF1::kNotDraw);
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

  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo_reco.root") );
  std::unique_ptr<THnD> mDh_reco(myFile->Get<THnD>("multi_data_histo_reco"));
  std::unique_ptr<THnD> mDh_gen(myFile->Get<THnD>("multi_data_histo_gen"));
  std::unique_ptr<THnD> mDh_diff_reco(myFile->Get<THnD>("multi_data_histo_diff_reco")); // it's reco - gen

  std::unique_ptr<TFile> myFile2( TFile::Open("multiD_histo_smear.root") );
  std::unique_ptr<THnD> mDh_smear(myFile2->Get<THnD>("multi_data_histo_smear"));
  std::unique_ptr<THnD> mDh_diff_smear(myFile2->Get<THnD>("multi_data_histo_diff_smear")); //it's smeared - gen
  std::unique_ptr<THnD> mDh_diff_squared_smear(myFile2->Get<THnD>("multi_data_histo_diff_squared_smear")); // it's smeared - gen  

  //these must match how the 5D histo was produced
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=6, nbinsmll=5, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, ptbinranges{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  //Initialise variables
  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0;
  double max_hist_mll_diff, max_hist_mll, value, error, diff_squared, evts_in_bin, sigma_mc=0.0, error_sigma_mc, error_diff_squared; 
  int middle_flag;
  string name;

  std::map<string, vector<double>> GenRecoFit;
  vector<double> fitresult;

  TH1D *mean_reco = new TH1D("mean", "reco-gen mll mean", 3, 0, 3);
  mean_reco->SetCanExtend(TH1::kAllAxes);

  TH1D *sigma_reco = new TH1D("sigma_reco", "reco-gen mll sigma", 3, 0, 3);
  sigma_reco->SetCanExtend(TH1::kAllAxes);

  TH1D *mean_smear = new TH1D("mean_smear", "smear-gen mll mean", 3, 0, 3);
  mean_smear->SetCanExtend(TH1::kAllAxes);
  mean_smear->SetMarkerStyle(kPlus);

  TH1D *sigma_smear = new TH1D("sigma_smear", "smear-gen mll sigma", 3, 0, 3);
  sigma_smear->SetCanExtend(TH1::kAllAxes);

  TH1D *gaus_integral = new TH1D("gaus_integral", "reco mll integral 75.5-105.0", 3, 0, 3);
  gaus_integral->SetCanExtend(TH1::kAllAxes);
  gaus_integral->SetMarkerStyle(kPlus);

  TH1D *occupancy = new TH1D("bin_occupancy", "bin occupancy", 3, 0, 3);
  occupancy->SetCanExtend(TH1::kAllAxes);
  occupancy->SetMarkerStyle(kPlus);

  VectorXd h_smear_minus_gen_vector(nbinsmll);
  Eigen::MatrixXd J(nbinsmll, nbinsmll), V_inv_sqrt(nbinsmll, nbinsmll);
  
  std::unique_ptr<TFile> f1( TFile::Open("control_bin_histo.root", "RECREATE") );
  std::unique_ptr<TFile> f2( TFile::Open("reco_gen_histos.root", "RECREATE") );
  ofstream f3("passed_regions.txt");

  auto multi_hist_proj_diff_reco = mDh_diff_reco->Projection(4);
  auto multi_hist_proj_reco = mDh_reco->Projection(4);
  auto multi_hist_proj_gen = mDh_gen->Projection(4);
  auto multi_hist_proj_smear = mDh_smear->Projection(4);
  auto multi_hist_proj_diff_smear = mDh_diff_smear->Projection(4);
  auto multi_hist_proj_diff_squared_smear = mDh_diff_squared_smear->Projection(4);
  f1->WriteObject(multi_hist_proj_diff_smear, "multi_hist_proj_diff_smear"); 
  f1->WriteObject(multi_hist_proj_diff_squared_smear, "multi_hist_proj_diff_squared_smear");

  total_nevents = multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff);
  std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";
  
  all_histos_count = 0.0;
  empty_histos_count = 0.0;
  hfrac = -1.0;  
  efrac = -1.0;
  remaining_nevents = 0.0;

  TCanvas *c1 = new TCanvas("c1","c1",800,600);  
  c1->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  auto leg2 = new TLegend(0.58, 0.68, 0.90, 0.90);

  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.035);
  leg2->SetFillColor(10);
  leg2->SetNColumns(1);
  leg2->SetHeader("");

  for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
    mDh_diff_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_gen->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_diff_squared_smear->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){   
      mDh_diff_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_gen->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_diff_squared_smear->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	mDh_diff_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_gen->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_diff_squared_smear->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	  mDh_diff_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_gen->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  mDh_diff_squared_smear->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);

	  all_histos_count++;
	  
	  delete gROOT->FindObject("multi_data_histo_diff_reco_proj_4");
	  multi_hist_proj_diff_reco = mDh_diff_reco->Projection(4);
	  nevents = multi_hist_proj_diff_reco->Integral(1,nbinsmll_diff);
	  if (nevents < 100.0){ // reject low stats
	    empty_histos_count++;
	  } else {
	    delete gROOT->FindObject("multi_data_histo_reco_proj_4");
	    multi_hist_proj_reco = mDh_reco->Projection(4);
	    fitresult = fitHisto(multi_hist_proj_reco, 4);
	    
	    if (fitresult[4] < 0.75){ // reject small gaus integral 
	      empty_histos_count++;
	    } else {
	      remaining_nevents += nevents;
	      
	      occupancy->SetBinError(occupancy->Fill(name.c_str(), nevents), 100);
	      
	      name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	      f3 << name << "\n";
	      
	      fitresult = fitHisto(multi_hist_proj_diff_reco, 4);
	      
	      GenRecoFit[name] = fitresult;
	      if (fitresult[0] > -90.0){ mean_reco->SetBinError(mean_reco->Fill(name.c_str(), fitresult[0]), fitresult[2]); }
	      if (fitresult[1] > -5.0){ sigma_reco->SetBinError(sigma_reco->Fill(name.c_str(), fitresult[1]), fitresult[3]); }
	      
	      multi_hist_proj_diff_reco->SetName(name.c_str());
	      multi_hist_proj_diff_reco->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      
	      c1->cd(1);
	      max_hist_mll_diff = -1.0;
	      multi_hist_proj_diff_reco->SetLineColor(kBlue);
	      
	      max_hist_mll_diff = multi_hist_proj_diff_reco->GetBinContent(multi_hist_proj_diff_reco->GetMaximumBin());
	      
	      delete gROOT->FindObject("multi_data_histo_diff_smear_proj_4");
	      multi_hist_proj_diff_smear = mDh_diff_smear->Projection(4);
	      fitresult = fitHisto(multi_hist_proj_diff_smear, 8);
	      if (fitresult[0] > -90.0){ mean_smear->SetBinError(mean_smear->Fill(name.c_str(), fitresult[0]), fitresult[2]); }
	      sigma_mc = fitresult[1];
	      error_sigma_mc = fitresult[3];
	      if (fitresult[1] > -5.0){ sigma_smear->SetBinError(sigma_smear->Fill(name.c_str(), fitresult[1]), fitresult[3]); }
	      multi_hist_proj_diff_smear->SetLineColor(kGreen);
	      if(max_hist_mll_diff < multi_hist_proj_diff_smear->GetBinContent(multi_hist_proj_diff_smear->GetMaximumBin())){
		multi_hist_proj_diff_smear->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
		multi_hist_proj_diff_smear->Draw();
		multi_hist_proj_diff_reco->Draw("SAME");
	      } else {
		multi_hist_proj_diff_reco->Draw();
		multi_hist_proj_diff_smear->Draw("SAME");
	      }
	      
	      leg1->Clear();
	      leg1->AddEntry(multi_hist_proj_diff_reco, "reco-gen", "l");
	      leg1->AddEntry(multi_hist_proj_diff_smear, "smeared-gen", "l");
	      leg1->Draw("");
	      
	      c1->cd(2);
	      max_hist_mll = -1.0;
	      middle_flag = 0;
	      
	      multi_hist_proj_reco->SetLineColor(kBlue);
	      fitresult = fitHisto(multi_hist_proj_reco, 4);
	      multi_hist_proj_reco->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	      max_hist_mll = multi_hist_proj_reco->GetBinContent(multi_hist_proj_reco->GetMaximumBin());
	      if (fitresult[4] > -100.0){ gaus_integral->SetBinError(gaus_integral->Fill(name.c_str(), fitresult[4]), 0.01); }
	      delete gROOT->FindObject("multi_data_histo_gen_proj_4");
	      multi_hist_proj_gen = mDh_gen->Projection(4);
	      fitresult = fitHisto(multi_hist_proj_gen, 2);
	      multi_hist_proj_gen->SetLineColor(kRed);
	      if(max_hist_mll < multi_hist_proj_gen->GetBinContent(multi_hist_proj_gen->GetMaximumBin())){
		max_hist_mll = multi_hist_proj_gen->GetBinContent(multi_hist_proj_gen->GetMaximumBin());
		middle_flag = 1;
	      }
	      delete gROOT->FindObject("multi_data_histo_smear_proj_4");
	      multi_hist_proj_smear = mDh_smear->Projection(4);
	      //fill vectors and variance for minimisation
	      V_inv_sqrt=MatrixXd::Zero(nbinsmll, nbinsmll);	      
	      for(int i=1; i<=nbinsmll; i++){
		h_smear_minus_gen_vector(i-1)=multi_hist_proj_smear->GetBinContent(i) - multi_hist_proj_gen->GetBinContent(i);
		V_inv_sqrt(i-1,i-1)=1/(multi_hist_proj_smear->GetBinErrorLow(i));
	      }
	      fitresult = fitHisto(multi_hist_proj_smear, 8);
	      multi_hist_proj_smear->SetLineColor(kGreen);
	      if(max_hist_mll < multi_hist_proj_smear->GetBinContent(multi_hist_proj_smear->GetMaximumBin())){
		multi_hist_proj_smear->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
		multi_hist_proj_smear->Draw();
		multi_hist_proj_gen->Draw("SAME");
		multi_hist_proj_reco->Draw("SAME");
	      } else if(middle_flag==1){
		multi_hist_proj_gen->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
		multi_hist_proj_gen->Draw();
		multi_hist_proj_smear->Draw("SAME");
		multi_hist_proj_reco->Draw("SAME");
	      } else {
		multi_hist_proj_reco->Draw();
		multi_hist_proj_gen->Draw("SAME");
		multi_hist_proj_smear->Draw("SAME");
	      }
	      
	      leg2->Clear();
	      leg2->AddEntry(multi_hist_proj_reco, "reco", "l");
	      leg2->AddEntry(multi_hist_proj_gen, "gen", "l");
	      leg2->AddEntry(multi_hist_proj_smear, "smeared gen", "l");
	      leg2->Draw("");
	      
	      c1->cd();
	      
	      f2->WriteObject(c1, name.c_str());
	      
	      //jacobian
	      delete gROOT->FindObject("jacobian");
	      delete gROOT->FindObject("diff_squared");
	      delete gROOT->FindObject("evts_in_bin");
	      TH1D *jac = new TH1D("jacobian", "jacobian", nbinsmll, 75.0, 105.0);
	      TH1D *dif_sq = new TH1D("diff_squared", "diff_squared", nbinsmll, 75.0, 105.0);
	      TH1D *ev_b = new TH1D("evts_in_bin", "evts_in_bin", nbinsmll, 75.0, 105.0);
	      //fill jacobian bin by bin
	      delete gROOT->FindObject("multi_data_histo_diff_squared_smear_proj_4");
	      multi_hist_proj_diff_squared_smear = mDh_diff_squared_smear->Projection(4);
	      J=MatrixXd::Zero(nbinsmll, nbinsmll);
	      for(int i=1; i<=nbinsmll; i++){
		diff_squared = multi_hist_proj_diff_squared_smear->GetBinContent(i);
		error_diff_squared = multi_hist_proj_diff_squared_smear->GetBinErrorLow(i);
		dif_sq->SetBinContent(i, diff_squared);
		evts_in_bin = multi_hist_proj_smear->GetBinContent(i);
		ev_b->SetBinContent(i, evts_in_bin);
		if (evts_in_bin > 0){
		  value =  diff_squared /  (evts_in_bin * sigma_mc * sigma_mc) - 1;
		  error = 1 / (evts_in_bin * sigma_mc*sigma_mc) * pow((4 * diff_squared*diff_squared * error_sigma_mc*error_sigma_mc / (sigma_mc*sigma_mc) + error_diff_squared*error_diff_squared) , 0.5);
		  //std::cout<< diff_squared <<" / "<< evts_in_bin<<" = "<<diff_squared / evts_in_bin<<"\n";
		} else {
		  value = 0.0; ///// Attention, value must be changed 
		  error = 20;
		}
		jac->SetBinContent(i, value);
		jac->SetBinError(i, error);
		J(i-1,i-1)=value;
	      }	
	      //write jacobian
	      f2->WriteObject(jac, ("jac" + name).c_str());
	      f2->WriteObject(dif_sq, ("dif_sq" + name).c_str());
	      f2->WriteObject(ev_b, ("ev_b" + name).c_str());
	      f2->WriteObject(multi_hist_proj_diff_squared_smear, ("dif_sq_proj" + name).c_str());
	      //solve for alpha
	      //Eigen::MatrixXd A = V_inv_sqrt*J;
	      //Eigen::MatrixXd b = V_inv_sqrt*h_smear_minus_gen_vector;
	      //Attention this is alpha, not alpha -1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      //Eigen::VectorXd alpha_vector = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	      //std::cout<<name.c_str()<<": A= "<<A<<", b= "<<b<<", alpha = "<< alpha_vector <<"\n";
	      //write alpha
	      //delete gROOT->FindObject("alpha");
	      //TH1D *alpha = new TH1D("alpha", "alpha", nbinsmll, 75.0, 105.0);
	      //std::cout<<"alpha"<<"\n";
	      //for(int i=1; i<=nbinsmll; i++){
	      //alpha->SetBinContent(i, alpha_vector(i-1));
	      //}
	      //f2->WriteObject(alpha, ("alpha" + name).c_str());
	      
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
