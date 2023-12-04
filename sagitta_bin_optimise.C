// Standalone code to fit for sagitta bias corrections

#include "TFile.h"
#include "TCanvas.h"

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

vector<float> fitHisto(TH1* histogram, int color){

  vector<float> fitresult;

  float mean = histogram->GetMean();
  float sigma = histogram->GetRMS();

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
    }
  }

  mean = gaussianFunc->GetParameter(1); 
  sigma = gaussianFunc->GetParameter(2); 
  fitresult.push_back(mean);
  fitresult.push_back(sigma);

  float mean_err = gaussianFunc->GetParError(1);
  float sigma_err = gaussianFunc->GetParError(2);

  fitresult.push_back(mean_err);
  fitresult.push_back(sigma_err);

  return fitresult;
}

int frame(){

  ROOT::EnableImplicitMT(128);

  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo.root") );
  std::unique_ptr<THnD> mDh_reco(myFile->Get<THnD>("multi_data_histo_reco"));
  std::unique_ptr<THnD> mDh_gen(myFile->Get<THnD>("multi_data_histo_gen"));
  std::unique_ptr<THnD> mDh_diff(myFile->Get<THnD>("multi_data_histo_diff"));

  //these must match how the 5D histo was produced
  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=18, nbinsmll=18, nbinseta=24, nbinspt=5;
  vector<double> etabinranges, ptbinranges{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  //Initialise variables
  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0;
  string name;

  std::map<string, vector<float>> GenRecoFit;
  vector<float> fitresult;

  TH1F *mean = new TH1F("mean", "gen-reco mll mean", 3, 0, 3);
  mean->SetCanExtend(TH1::kAllAxes);

  TH1F *sigma = new TH1F("sigma", "gen-reco mll sigma", 3, 0, 3);
  sigma->SetCanExtend(TH1::kAllAxes);

  std::unique_ptr<TFile> f1( TFile::Open("control_bin_histo.root", "RECREATE") );
  std::unique_ptr<TFile> f2( TFile::Open("reco_gen_histos.root", "RECREATE") );

  auto multi_hist_proj_diff = mDh_diff->Projection(4);
  auto multi_hist_proj_reco = mDh_reco->Projection(4);
  auto multi_hist_proj_gen = mDh_gen->Projection(4);

  total_nevents = multi_hist_proj_diff->Integral(1,nbinsmll_diff);
  std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";
  
  all_histos_count = 0.0;
  empty_histos_count = 0.0;
  hfrac = -1.0;  
  efrac = -1.0;
  remaining_nevents = 0.0;

  TCanvas *c1 = new TCanvas("c1","c1",800,600);  
  c1->Divide(2,1);

  for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
    mDh_diff->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_reco->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    mDh_gen->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
    for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){   
      mDh_diff->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_reco->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      mDh_gen->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
      for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	mDh_diff->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	mDh_reco->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
        mDh_gen->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	  mDh_diff->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_reco->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          mDh_gen->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	  all_histos_count++;
	  
	  delete gROOT->FindObject("multi_data_histo_diff_proj_4");
	  multi_hist_proj_diff = mDh_diff->Projection(4);
	  nevents = multi_hist_proj_diff->Integral(1,nbinsmll_diff);
	  //std::cout<<nevents<<"\n";

	  if (nevents < 100.0){ 
	    empty_histos_count++;
	  } else {
	    remaining_nevents += nevents;
	    
	    name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	    
	    fitresult = fitHisto(multi_hist_proj_diff, 2);
	    GenRecoFit[name] = fitresult;
	    mean->SetBinError(mean->Fill(name.c_str(), fitresult[0]), fitresult[2]);
	    sigma->SetBinError(sigma->Fill(name.c_str(), fitresult[1]), fitresult[3]);
	    
	    multi_hist_proj_diff->SetName(name.c_str());
	    multi_hist_proj_diff->SetTitle(("reco - gen mll " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	    //delete gROOT->FindObject("c1");
	    c1->cd(1);
	    multi_hist_proj_diff->Draw();
	    //f2->WriteObject(multi_hist_proj_diff, name.c_str());	    

	    c1->cd(2);
	    delete gROOT->FindObject("multi_data_histo_reco_proj_4");
	    multi_hist_proj_reco = mDh_reco->Projection(4);
	    multi_hist_proj_reco->SetLineColor(kBlue);
	    fitresult = fitHisto(multi_hist_proj_reco, 4);
            multi_hist_proj_reco->SetTitle(("reco and gen mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	    multi_hist_proj_reco->Draw();
	    delete gROOT->FindObject("multi_data_histo_gen_proj_4");
	    multi_hist_proj_gen = mDh_gen->Projection(4);
	    fitresult = fitHisto(multi_hist_proj_gen, 2);
	    multi_hist_proj_gen->SetLineColor(kRed);
            multi_hist_proj_gen->Draw("SAME");

	    auto leg = new TLegend(0.10, 0.68, 0.70, 0.90);
	    leg->SetFillStyle(0);
	    leg->SetBorderSize(0);
	    leg->SetTextSize(0.035);
	    leg->SetFillColor(10);
	    leg->SetNColumns(1);
	    leg->SetHeader("");

	    leg->AddEntry(multi_hist_proj_reco, "reco", "l");
	    leg->AddEntry(multi_hist_proj_gen, "gen", "l");
	    leg->Draw("");

	    c1->cd();

	    f2->WriteObject(c1, name.c_str());

	  }

	}
      } 
    }
  }

  hfrac = empty_histos_count / all_histos_count;
  efrac = remaining_nevents / total_nevents;
  std::cout<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";

  mean->SetStats(0);
  mean->LabelsDeflate();
  f1->WriteObject(mean, "mean_diff");
  
  sigma->SetStats(0);
  sigma->LabelsDeflate();
  f1->WriteObject(sigma, "sigma_diff");

  return 0; 
}
