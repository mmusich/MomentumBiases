// Standalone code to fit for sagitta

#include "TFile.h"
#include "TCanvas.h"

using namespace ROOT;
using namespace ROOT::VecOps;

string stringify_name(int a, int b, int c, int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

int frame(){

  ROOT::EnableImplicitMT(128);

  //  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo.root") );
  //  std::unique_ptr<THnD> mDh(myFile->Get<THnD>("multi_data_frame"));

  RDataFrame df("Events", "tree_output.root");

  float ptlow=25.0, pthigh=55.0;
  int nbinseta=24, nbinspt;
  vector<float> etabinranges, ptbinranges;
  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  double total_events_count=0.0, entries = 0.0, all_histos_count = 0.0, remaining_events_count = 0.0, empty_histos_count = 0.0, hfrac = -1.0, efrac = -1.0;
  int bestnbin = 0, bestemptycount = 10000, threshold = 0;
  string name, title, title_seed = "reco - gen mll ";

  //optimisation starts
  for(nbinspt=6; nbinspt<=15; nbinspt++){
    delete gROOT->FindObject("multi_data_frame");
    auto mDh = df.HistoND<float, float, float, float, float, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, 24}, {-2.4,ptlow,-2.4,ptlow,-6.0}, {2.4,pthigh,2.4,pthigh,6.0}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","genWeight"});

    total_events_count = mDh->Projection(4)->Integral(1,24);
    std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_events_count << " entries in the inclusive projection \n";

    ptbinranges.clear();
    std::cout<<"ptbinranges = [";
    for (int i=0; i<=nbinspt; i++){ptbinranges.push_back(ptlow + i * (pthigh-ptlow)/nbinspt); std::cout<<ptbinranges[i]<<", ";}
    std::cout<<"] \n";

    all_histos_count = 0.0;
    empty_histos_count = 0.0;
    hfrac = -1.0;  
    efrac = -1.0;
    remaining_events_count = 0.0;
    
    for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
      mDh->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
      for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){
	mDh->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
	for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	  mDh->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	  for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	    mDh->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	    all_histos_count++;

	    delete gROOT->FindObject("multi_data_frame_proj_4");
	    entries = mDh->Projection(4)->Integral(1,24);
	    //std::cout<<entries<<"\n";

	    if (entries < 30.0){ 
	      empty_histos_count++;
	    } else {
	      remaining_events_count += entries;
	    }
	  }
	} 
      }
    }
    hfrac = empty_histos_count / all_histos_count;
    efrac = remaining_events_count / total_events_count;
    std::cout<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";  
  }

  return 0; 
}
