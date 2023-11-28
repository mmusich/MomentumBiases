// Standalone code to fit for sagitta

#include "TFile.h"
#include "TCanvas.h"

using namespace ROOT;
using namespace ROOT::VecOps;

string stringify_name(int a, int b, int c, int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

string stringify_title(int a, int b, int c, int d, vector<double> e, vector<double> p){
  string title_seed = "reco - gen pT ", comma = ",";
  string txt1 = "eta+ in [", txt2="] pt+ in [", txt3="] eta- in [", txt4="] pt- in [", txt5 = "]";
  return title_seed+txt1+to_string(e[a-1]).substr(0,4)+comma+to_string(e[a]).substr(0,4)+txt2+to_string(p[b-1]).substr(0,4)+comma+to_string(p[b]).substr(0,4)+txt3+to_string(e[c-1]).substr(0,4)+comma+to_string(e[c]).substr(0,4)+txt4+to_string(p[d-1]).substr(0,4)+comma+to_string(p[d]).substr(0,4)+txt5;
}

int frame(){

  ROOT::EnableImplicitMT(128);

  //  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo.root") );
  //  std::unique_ptr<THnD> mDh(myFile->Get<THnD>("multi_data_frame"));

  RDataFrame df("Events", "tree_output.root");

  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=24, nbinseta=24, nbinspt;
  vector<double> etabinranges, ptbinranges, mll_diffbinranges;

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";
  //automatise this
  double myetaboundaries[] = {-2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4};

  for (int i=0; i<=nbinsmll_diff; i++){mll_diffbinranges.push_back(-6.0 + i * 12.0/nbinsmll_diff);}

  double total_nevents=0.0, nevents = 0.0, pt_eta_nevents = 0.0, all_histos_count = 0.0, remaining_nevents = 0.0, empty_histos_count = 0.0, hfrac = -1.0, efrac = -1.0;
  int bestnbin = 0, bestemptycount = 10000, threshold = 0;
  string name, title;

  //optimisation starts
  //for(nbinspt=6; nbinspt<=15; nbinspt++){
  for(nbinspt=6; nbinspt<=6; nbinspt++){  
    ptbinranges.clear();
    //Fixed bin size
    //std::cout<<"ptbinranges = [";
    //for (int i=0; i<=nbinspt; i++){ptbinranges.push_back(ptlow + i * (pthigh-ptlow)/nbinspt); std::cout<<ptbinranges[i]<<", ";}
    //std::cout<<"] \n";
    
    //Variable bin size
    std::cout<<"ptbinranges = [";
    double myptboundaries[] = {25.0, 34.0, 38.0, 40.0, 42.0, 46.0, 55.0};
    ptbinranges.assign (myptboundaries,myptboundaries+nbinspt+1); //CHECK why need myptboundaries+nbinspt+1
    for (int i=0; i<=nbinspt; i++){ std::cout<<ptbinranges[i]<<", ";}
    std::cout<<"] \n";

    delete gROOT->FindObject("pt_eta_pos");
    TH2D* pt_eta_pos = new TH2D("pt_eta_pos", "pt eta mu+", nbinseta, myetaboundaries, nbinspt, myptboundaries); //CHECK LOW EDGES

    delete gROOT->FindObject("multi_data_frame");
    auto mDh = df.HistoND<float, float, float, float, float, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","genWeight"}); //CHECK LOW EDGES

    delete gROOT->FindObject("multi_data_frame_proj_4");
    auto multi_hist_proj = mDh->Projection(4);
    total_nevents = multi_hist_proj->Integral(1,24);
    std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";

    all_histos_count = 0.0;
    empty_histos_count = 0.0;
    hfrac = -1.0;  
    efrac = -1.0;
    remaining_nevents = 0.0;
    
    TFile f1("reco_gen_histos.root","recreate");
    for (int pos_eta_bin=1; pos_eta_bin<=nbinseta; pos_eta_bin++){
      mDh->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
      for (int pos_pt_bin=1; pos_pt_bin<=nbinspt; pos_pt_bin++){
	mDh->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);

	pt_eta_nevents=0.0;

	for (int neg_eta_bin=1; neg_eta_bin<=nbinseta; neg_eta_bin++){
	  mDh->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
	  for (int neg_pt_bin=1; neg_pt_bin<=nbinspt; neg_pt_bin++){
	    mDh->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
	    all_histos_count++;

	    delete gROOT->FindObject("multi_data_frame_proj_4");
	    multi_hist_proj = mDh->Projection(4);
	    nevents = multi_hist_proj->Integral(1,24);
	    //std::cout<<nevents<<"\n";
	    pt_eta_nevents += nevents;

	    if (nevents < 30.0){ 
	      empty_histos_count++;
	    } else {
	      remaining_nevents += nevents;
	      
	      name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	      title = stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges);
	      multi_hist_proj->SetName(name.c_str());
	      multi_hist_proj->SetTitle(title.c_str());
	      multi_hist_proj->Write();
	    }
	  }
	} 
	pt_eta_pos->SetBinContent(pos_eta_bin, pos_pt_bin, pt_eta_nevents);
      }
    }
    hfrac = empty_histos_count / all_histos_count;
    efrac = remaining_nevents / total_nevents;
    std::cout<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";
    f1.Close();

    TFile f2("control_bin_histo.root","recreate");
    pt_eta_pos->Write(); 
    f2.Close();
  }

  return 0; 
}
