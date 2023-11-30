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
  string title_seed = "reco - gen mll ", comma = ",";
  string txt1 = "eta+ in [", txt2="] pt+ in [", txt3="] eta- in [", txt4="] pt- in [", txt5 = "]";
  return title_seed+txt1+to_string(e[a-1]).substr(0,4)+comma+to_string(e[a]).substr(0,4)+txt2+to_string(p[b-1]).substr(0,4)+comma+to_string(p[b]).substr(0,4)+txt3+to_string(e[c-1]).substr(0,4)+comma+to_string(e[c]).substr(0,4)+txt4+to_string(p[d-1]).substr(0,4)+comma+to_string(p[d]).substr(0,4)+txt5;
}

int frame(){

  ROOT::EnableImplicitMT(128);

  //  std::unique_ptr<TFile> myFile( TFile::Open("multiD_histo.root") );
  //  std::unique_ptr<THnD> mDh(myFile->Get<THnD>("multi_data_frame"));

  RDataFrame df("Events", "tree_output.root");

  double ptlow=25.0, pthigh=55.0;
  int nbinsmll_diff=18, nbinseta=24, nbinspt;
  vector<double> etabinranges, ptbinranges, mll_diffbinranges;

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta); std::cout<<etabinranges[i]<<", ";}
  std::cout<<"] \n";

  double myetaboundaries[etabinranges.size()];
  for (int i=0; i<etabinranges.size(); i++){
    myetaboundaries[i] = etabinranges[i];
  }

  for (int i=0; i<=nbinsmll_diff; i++) mll_diffbinranges.push_back(-6.0 + i * 12.0/nbinsmll_diff);

  //Initialise variables
  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0;
  string name, title;

  //pT bin optimisation starts
  for(nbinspt=5; nbinspt<=5; nbinspt++){
    TFile f1("control_bin_histo.root","recreate");

    delete gROOT->FindObject("pt_pos_uni");
    auto pt_pos_uni = df.Histo1D({"pt_pos_uni", "pt mu+", nbinspt*3, ptlow, pthigh},"posTrackPt");
    pt_pos_uni->Write();
    // Get quartiles
    double xq[nbinspt+1], myptboundaries[nbinspt+1];
    for (int i=0;i<=nbinspt;i++) xq[i] = float(i)/nbinspt;
    pt_pos_uni->GetQuantiles(nbinspt+1,myptboundaries, xq);

    ptbinranges.clear();
    std::cout<<"ptbinranges = [";
    for (int i=0; i<=nbinspt; i++){
      ptbinranges.push_back(myptboundaries[i]);
      std::cout<<ptbinranges[i]<<", ";
    }
    std::cout<<"] \n";
        
    //TH1 in pT+ with variable bin size -> should be uniform
    delete gROOT->FindObject("pt_pos");
    auto pt_pos = df.Histo1D({"pt_pos", "pt mu+", nbinspt, myptboundaries},"posTrackPt");
    pt_pos->Write();

    delete gROOT->FindObject("pt_eta_pos");
    auto pt_eta_pos = df.Histo2D({"pt_eta_pos", "pt eta mu+", nbinseta, myetaboundaries, nbinspt, myptboundaries},"posTrackEta", "posTrackPt");
    pt_eta_pos->Write();

    f1.Close();

    delete gROOT->FindObject("multi_data_frame");
    auto mDh = df.HistoND<float, float, float, float, float, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","mll_diff","genWeight"});

    delete gROOT->FindObject("multi_data_frame_proj_4");
    auto multi_hist_proj = mDh->Projection(4);
    total_nevents = multi_hist_proj->Integral(1,nbinsmll_diff);
    std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";

    all_histos_count = 0.0;
    empty_histos_count = 0.0;
    hfrac = -1.0;  
    efrac = -1.0;
    remaining_nevents = 0.0;
    
    TFile f2("reco_gen_histos.root","recreate");
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
	    multi_hist_proj = mDh->Projection(4);
	    nevents = multi_hist_proj->Integral(1,nbinsmll_diff);
	    //std::cout<<nevents<<"\n";

	    if (nevents < 1000.0){ 
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
      }
    }
    f2.Close();
    hfrac = empty_histos_count / all_histos_count;
    efrac = remaining_nevents / total_nevents;
    std::cout<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";

  }

  return 0; 
}
