// Standalone code to fit for sagitta

#include "TFile.h"
#include "TCanvas.h"

using namespace ROOT;
using namespace ROOT::VecOps;

string stringify_name(int a, int b, int c, int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

string stringify_title(int a, int b, int c, int d, vector<float> e, vector<float> p){
  string title_seed = "reco - gen pT ", comma = ",";
  string txt1 = "eta+ in [", txt2="] pt+ in [", txt3="] eta- in [", txt4="] pt- in [", txt5 = "]";
  return title_seed+txt1+to_string(e[a-1]).substr(0,4)+comma+to_string(e[a]).substr(0,4)+txt2+to_string(p[b-1]).substr(0,4)+comma+to_string(p[b]).substr(0,4)+txt3+to_string(e[c-1]).substr(0,4)+comma+to_string(e[c]).substr(0,4)+txt4+to_string(p[d-1]).substr(0,4)+comma+to_string(p[d]).substr(0,4)+txt5;
}

vector<float> fitHisto(TH1* histogram){

  vector<float> fitresult;

  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", -6, 6);
  histogram->Fit(gaussianFunc, "Q");

  float mean = gaussianFunc->GetParameter(1); 
  float sigma = gaussianFunc->GetParameter(2); 
  fitresult.push_back(mean);
  fitresult.push_back(sigma);

  return fitresult;
}

int frame(){

  ROOT::EnableImplicitMT(128);

  int nbinseta=24, nbinspt=30;
  float ptlow=25.0, pthigh=55.0;
  //larger bins while debugging
  //int nbinseta=6, nbinspt=2;
  
  //number of eta pt bins should match or be larger than the one chosen later  
  //auto multi_hist = d5.HistoND<float, float, float, float, RVecF, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, 24}, {-2.4,ptlow,-2.4,ptlow,-6}, {2.4,pthigh,2.4,pthigh,6}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff","genWeight"});

  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo.root") );
  std::unique_ptr<THnD> mDh(myFile->Get<THnD>("multi_data_frame"));

  TH1F *mean = new TH1F("mean", "gen-reco mean", 3, 0, 3);
  mean->SetCanExtend(TH1::kAllAxes);

  TH1F *sigma = new TH1F("sigma", "gen-reco sigma", 3, 0, 3);
  sigma->SetCanExtend(TH1::kAllAxes);

  std::cout << "There are " << mDh->Projection(4)->GetEntries() << " entries in the inclusive projection \n";
 
  TFile f("reco_gen_histos.root","recreate");
  
  auto multi_hist_proj = mDh->Projection(4);
  multi_hist_proj->SetName("inclusive_proj_pt_diff");
  multi_hist_proj->SetTitle("reco - gen pT inclusive");
  //multi_hist_proj->Write();
  
  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", -6, 6);
  multi_hist_proj->Fit(gaussianFunc, "Q");
  multi_hist_proj->Write();

  multi_hist_proj = mDh->Projection(0);
  multi_hist_proj->SetName("inclusive_proj_pos_eta");
  multi_hist_proj->SetTitle("pos eta inclusive");
  multi_hist_proj->Write();

  multi_hist_proj = mDh->Projection(1);
  multi_hist_proj->SetName("inclusive_proj_pos_pt");
  multi_hist_proj->SetTitle("pos pt inclusive");
  multi_hist_proj->Write();

  multi_hist_proj = mDh->Projection(2);
  multi_hist_proj->SetName("inclusive_proj_neg_eta");
  multi_hist_proj->SetTitle("neg eta inclusive");
  multi_hist_proj->Write();

  multi_hist_proj = mDh->Projection(3);
  multi_hist_proj->SetName("inclusive_proj_neg_pt");
  multi_hist_proj->SetTitle("neg pt inclusive");
  multi_hist_proj->Write();
 
  std::map<string, vector<float>> GenRecoFit;
  vector<float> fitresult;

  string name, title, title_seed = "reco - gen pT "; 
  
  TH2F* empty_histos = new TH2F("empty_histos", " fraction empty histos mu+", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh);
  TH2F* total_histos = new TH2F("total_histos", " total histos", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh);
  int entries = 0, all_histos_count = 0, empty_histos_count = 0;

  vector<float> etabinranges, ptbinranges;
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta);}
  //optimisation starts
  nbinspt = 10;
  for (int i=0; i<=nbinspt; i++){ptbinranges.push_back(ptlow + i * (pthigh-ptlow)/nbinspt);}
  
  for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
    for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){
      for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
          mDh->GetAxis(0)->SetRangeUser(etabinranges[pos_eta_bin - 1] + 0.00001, etabinranges[pos_eta_bin] - 0.0001);
	  mDh->GetAxis(1)->SetRangeUser(ptbinranges[pos_pt_bin - 1] + 0.00001, ptbinranges[pos_pt_bin] - 0.0001);
          mDh->GetAxis(2)->SetRangeUser(etabinranges[neg_eta_bin - 1] + 0.00001, etabinranges[neg_eta_bin] - 0.0001);
          mDh->GetAxis(3)->SetRangeUser(ptbinranges[neg_pt_bin - 1] + 0.00001, ptbinranges[neg_pt_bin] - 0.0001);
          multi_hist_proj = mDh->Projection(4);

          entries = multi_hist_proj->GetEntries();
          //total_histos->Fill(etabinranges[pos_eta_bin-1]+0.001, ptbinranges[pos_pt_bin-1]+0.001);
          all_histos_count++;
          if (entries < 50){ 
            //empty_histos->Fill(etabinranges[pos_eta_bin-1]+0.001, ptbinranges[pos_pt_bin-1]+0.001);
            empty_histos_count++;
            continue;
          } else {
            
	    name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
            title = stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges);
            multi_hist_proj->SetName(name.c_str());
            multi_hist_proj->SetTitle(title.c_str());

	    fitresult = fitHisto(multi_hist_proj);
            GenRecoFit[name] = fitresult;
            mean->Fill(name.c_str(), fitresult[0]);
	    sigma->Fill(name.c_str(), fitresult[1]);
            multi_hist_proj->Write(name.c_str());
            
	  }
	}
      } 
    }
  } 
 
  /*
  empty_histos->Divide(total_histos);
  empty_histos->Write("empty_histos");    
  auto c = new TCanvas("c","c");
  c->SetRightMargin(0.20);
  empty_histos->Draw("COLZ");
  c->SaveAs("empty_histos.pdf");
  
  mean->SetStats(0);
  mean->LabelsDeflate();
  mean->Write("mean");
  
  sigma->SetStats(0);
  sigma->LabelsDeflate();
  sigma->Write("sigma");
  */

  f.Close();

  return 0; 

}
