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

  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=24, nbinseta=2, nbinspt;
  vector<double> etabinranges, ptbinranges, mll_diffbinranges;

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  for (int i=0; i<=nbinsmll_diff; i++){mll_diffbinranges.push_back(-6.0 + i * 12.0/nbinsmll_diff);}

  double total_nevents=0.0, nevents = 0.0, all_histos_count = 0.0, remaining_nevents = 0.0, empty_histos_count = 0.0, hfrac = -1.0, efrac = -1.0;
  int bestnbin = 0, bestemptycount = 10000, threshold = 0;
  string name, title, title_seed = "reco - gen mll ";

  //define TH2 for pt eta 1

  //optimisation starts
  //for(nbinspt=6; nbinspt<=15; nbinspt++){
  for(nbinspt=3; nbinspt<=3; nbinspt++){  
    ptbinranges.clear();
    //Fixed bin size
    std::cout<<"ptbinranges = [";
    for (int i=0; i<=nbinspt; i++){ptbinranges.push_back(ptlow + i * (pthigh-ptlow)/nbinspt); std::cout<<ptbinranges[i]<<", ";}
    std::cout<<"] \n";
    
    delete gROOT->FindObject("multi_data_frame");
    auto mDh = df.HistoND<float, float, float, float, float, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","genWeight"});

    delete gROOT->FindObject("multi_data_frame_proj_4");
    total_nevents = mDh->Projection(4)->Integral(1,24);
    std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";

    all_histos_count = 0.0;
    empty_histos_count = 0.0;
    hfrac = -1.0;  
    efrac = -1.0;
    remaining_nevents = 0.0;
    
    for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
      mDh->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
      for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){
	mDh->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
	// empty the TH2

	for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	  mDh->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	  for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	    mDh->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	    all_histos_count++;

	    delete gROOT->FindObject("multi_data_frame_proj_4");
	    nevents = mDh->Projection(4)->Integral(1,24);
	    //std::cout<<nevents<<"\n";

	    if (nevents < 30.0){ 
	      empty_histos_count++;
	    } else {
	      remaining_nevents += nevents;
	    }
	    //what do i fill the 2d histo with? nevents in this pt eta region

	  }
	} 
	//fill the 2d histo of mu 1 in current pt eta bin

      }
    }
    hfrac = empty_histos_count / all_histos_count;
    efrac = remaining_nevents / total_nevents;
    std::cout<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";

    //plot pt vs eta 1 
    //write in some file

  }

  return 0; 
}
