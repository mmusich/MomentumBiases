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

  //int nbinseta=24, nbinspt=50;
  float ptlow=25.0, pthigh=55.0;
  //larger bins while debugging
  int nbinseta=4, nbinspt=2;
  
  //number of eta pt bins should match or be larger than the one chosen later  
  //auto multi_hist = d5.HistoND<float, float, float, float, RVecF, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, 24}, {-2.4,ptlow,-2.4,ptlow,-6}, {2.4,pthigh,2.4,pthigh,6}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff","genWeight"});

  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo.root") );
  std::unique_ptr<THnD> mDh(myFile->Get<THnD>("multi_data_frame"));

  float total_events_count = mDh->Projection(4)->GetEntries();

  std::cout << "There are " << total_events_count << " entries in the inclusive projection \n";
 
  string name, title, title_seed = "reco - gen pT "; 
  
  float entries = 0.0, all_histos_count = 0.0, remaining_events_count = 0.0, empty_histos_count = 0.0, hfrac = -1.0, efrac = -1.0;

  vector<float> etabinranges, ptbinranges;
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta);}
  //optimisation starts

  int bestnbin = 0, bestemptycount = 10000, threshold = 0;

  for (nbinspt = 18; nbinspt <= 20; nbinspt++){

    ptbinranges.clear();
    for (int i=0; i<=nbinspt; i++){ptbinranges.push_back(ptlow + i * (pthigh-ptlow)/nbinspt);}
    all_histos_count = 0.0;
    empty_histos_count = 0.0;
    hfrac = -1.0;  
    efrac = -1.0;
    remaining_events_count = 0.0;

    for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
      for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){
	for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	  for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	    delete gROOT->FindObject("multi_data_frame_proj_4");
	    mDh->GetAxis(0)->SetRangeUser(etabinranges[pos_eta_bin - 1] + 0.00001, etabinranges[pos_eta_bin] - 0.0001);
	    mDh->GetAxis(1)->SetRangeUser(ptbinranges[pos_pt_bin - 1] + 0.00001, ptbinranges[pos_pt_bin] - 0.0001);
	    mDh->GetAxis(2)->SetRangeUser(etabinranges[neg_eta_bin - 1] + 0.00001, etabinranges[neg_eta_bin] - 0.0001);
	    mDh->GetAxis(3)->SetRangeUser(ptbinranges[neg_pt_bin - 1] + 0.00001, ptbinranges[neg_pt_bin] - 0.0001);
	    all_histos_count++;
	    
	    entries = mDh->Projection(4)->GetEntries();
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
