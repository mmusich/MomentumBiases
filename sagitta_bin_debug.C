// Standalone code to fit for sagitta

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"

using namespace ROOT;
using namespace ROOT::VecOps;

// MuonGenMatchedIndex
// works as it is because there is only one muon pair selected  per event
// TODO try a match by deltaR instead
RVecF MuonGenMatchedIndex(RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt){
  int posindex = 999;
  int negindex = 999;
  RVecF genpt;
  for(int i=0;i<GenPart_status.size();i++){
    if (GenPart_status[i] == 1 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 23){
      if (GenPart_pdgId[i] == -13){ // mu(-) has PDGID 13
        posindex = i;
      } else if(GenPart_pdgId[i] == 13) {
        negindex = i;
      }
    }
  }
  genpt.push_back(GenPart_pt[posindex]);
  genpt.push_back(GenPart_pt[negindex]);
    
  return genpt;
}

//dxy_significance

RVecF dxy_significance(RVecF Muon_dxy, RVecF Muon_dxyErr){
  return abs(Muon_dxy)/Muon_dxyErr;
}

// MuonisGood

RVecB MuonisGood(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
       if (Muon_pt[i] > 10 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
            muonisgood.push_back(1);
       }
  }
  return muonisgood;
}

//Pairs

RVec<std::tuple<int,int,float>> pairs(RVecF Muon_pt, RVecF Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, float rest_mass){

  RVec<std::tuple<int,int,float>> pairs; // <pos_muon_index, neg_muon_index, pair_inv_mass>
  std::tuple<int,int,float> temp;

  for(int i=1;i<Muon_pt.size();i++){
  for(int j=0;j<i;j++){
            if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
                 // Define opening angle
                 float gamma_angle = 0;
                 gamma_angle = acos(sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*cos(Muon_phi[i])*cos(Muon_phi[j]) + sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*sin(Muon_phi[i])*sin(Muon_phi[j]) + cos(2*atan(exp((-1)*Muon_eta[i])))*cos(2*atan(exp((-1)*Muon_eta[j]))));
                 if(gamma_angle > 3.141592653/4){
                      TLorentzVector firstTrack, secondTrack, mother;
                      firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);  
                      secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
                      mother = firstTrack + secondTrack;
                      float mll = mother.M();
                      if(75<mll && mll<105){
                           if(Muon_charge[i]==1){
                                temp=make_tuple(i,j,mll);
                           } else {
                                temp=make_tuple(j,i,mll);
                           }
                           pairs.push_back(temp);
                      }
                 }
            }
       }
  }
  if(pairs.size()>1){
       float diff=100.0;
       int best=0;
       for(int i=0;i<pairs.size();i++){
            if(abs(get<2>(pairs.at(i))-91)<diff){
                 diff=(abs(get<2>(pairs.at(i))-91));
                 best=i;
            }
       }
       RVec<std::tuple<int,int,float>> pairss;
       pairss.push_back(pairs.at(best));
       return pairss;
  }
  return pairs;
}

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

  TFile *f1 = new TFile("snapshot_output.root","RECREATE");
  TFile *f2 = new TFile("reco_gen_histos.root","RECREATE");
  
  TChain chain("Events");
  //chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_773.root");
  //chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_772.root");
  
  string line;

  ifstream file("MCfilenames.txt");
  if (file.is_open()) {
    while (getline(file, line)) {
      chain.Add(line.c_str());
    }
    file.close();
  }  

  RDataFrame df(chain);

  auto d1 = df.Filter("HLT_IsoMu24 == 1")
                .Filter("nMuon >= 2")
                .Filter("PV_npvsGood >= 1");

  auto d2 = d1.Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)")
                .Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)")
                .Define("pairs", "pairs(Muon_pt, Muon_charge, Muon_eta, Muon_phi, MuonisGood, Muon_dxy, Muon_dz, 0.105658)"); // muMass = 0.105658 GeV
  
  auto d3 = d2.Filter("pairs.size()==1");

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  auto d4 = d3.Define("posTrackPt","float posTrackPt; posTrackPt=Muon_pt[get<0>(pairs.at(0))]; return posTrackPt;")
              .Define("negTrackPt","float negTrackPt; negTrackPt=Muon_pt[get<1>(pairs.at(0))]; return negTrackPt;")
              .Define("recoPt","RVecF recoPt; recoPt.push_back(posTrackPt); recoPt.push_back(negTrackPt); return recoPt;")
              .Define("genPt", "MuonGenMatchedIndex(GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)")
              .Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs.at(0))]; return posTrackEta;")
	      .Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs.at(0))]; return negTrackEta;")
              .Define("mll","RVecF temporary; for(int i=0; i<pairs.size(); i++){temporary.push_back(get<2>(pairs.at(i)));} return temporary;");

  auto d5 = d4.Define("ptDiff","RVecF ptDiff; ptDiff.push_back(recoPt[0] - genPt[0]); ptDiff.push_back(recoPt[1] - genPt[1]); return ptDiff;")
              .Define("posPtDiff","float posptDiff; posptDiff=recoPt[0] - genPt[0]; return posptDiff;")
              .Define("negPtDiff","float negptDiff; negptDiff=recoPt[1] - genPt[1]; return negptDiff;");

  d5.Snapshot("Events", "snapshot_output.root", {"Muon_pt", "Muon_charge", "recoPt", "genPt", "posPtDiff", "negPtDiff", "GenPart_pdgId", "GenPart_pt"});

  int nbinseta=1, nbinspt=1;
  float ptlow=25.0, pthigh=55.0;
  //larger bins while debugging
  //int nbinseta=6, nbinspt=2;
  /*
  auto posPt_hist = d5.Histo3D<float, float, float, float>({"posPt_hist", "posPt_hist", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh, 24, -6., 6.},"posTrackEta", "posTrackPt", "posPtDiff", "genWeight");
  auto negPt_hist = d5.Histo3D<float, float, float, float>({"negPt_hist", "negPt_hist", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh, 24, -6., 6.},"negTrackEta", "negTrackPt", "negPtDiff", "genWeight");

  auto multi_hist = d5.HistoND<float, float, float, float, RVecF, float>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, 24}, {-2.4,ptlow,-2.4,ptlow,-6}, {2.4,pthigh,2.4,pthigh,6}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff","genWeight"});
  */
  // don't fill with weights for now
  auto posPt_hist = d5.Histo3D<float, float, float>({"posPt_hist", "posPt_hist", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh, 24, -6., 6.},"posTrackEta", "posTrackPt", "posPtDiff");
  auto negPt_hist = d5.Histo3D<float, float, float>({"negPt_hist", "negPt_hist", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh, 24, -6., 6.},"negTrackEta", "negTrackPt", "negPtDiff");

  auto multi_hist = d5.HistoND<float, float, float, float, RVecF>({"multi_data_frame", "multi_data_frame", 5, {nbinseta, nbinspt, nbinseta, nbinspt, 24}, {-2.4,ptlow,-2.4,ptlow,-6.0}, {2.4,pthigh,2.4,pthigh,6.0}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff"});

  auto posPt_hist2 = d5.Histo2D<float, float>({"posPt_hist2", "posPt_hist2", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh},"posTrackEta", "posTrackPt");
  auto negPt_hist2 = d5.Histo2D<float, float>({"negPt_hist2", "negPt_hist2", nbinseta, -2.4, 2.4, nbinspt, ptlow, pthigh},"negTrackEta", "negTrackPt");

  vector<float> etabinranges, ptbinranges;
  for (int i=0; i<=nbinseta; i++){etabinranges.push_back(-2.4 + i * 4.8/nbinseta);}

  nbinspt = 1;
  for (int i=0; i<=nbinspt; i++){ptbinranges.push_back(ptlow + i * (pthigh-ptlow)/nbinspt);}

  std::cout.precision(9);
  std::cout << std::scientific;

  double initial_entries = 0.0, initial_integral = 0.0;
  auto posPt_hist_proj = posPt_hist->ProjectionZ();
  initial_integral += posPt_hist_proj->Integral(1,24);
  std::cout << "There are " << std::setprecision(9) << posPt_hist_proj->Integral(1,24) << " entries in the pos inclusive projection \n";
 
  auto negPt_hist_proj = negPt_hist->ProjectionZ();
  initial_integral += negPt_hist_proj->Integral(1,24);
  std::cout << "There are " << std::setprecision(9) << negPt_hist_proj->Integral(1,24) << " entries in the neg inclusive projection, total: "<< std::setprecision(9) << initial_integral <<"\n";
  
  auto multi_hist_proj = multi_hist->Projection(4);
  initial_entries = multi_hist_proj->Integral(1,24);
  std::cout << "There are " << std::setprecision(9) << initial_entries << " entries in the inclusive projection \n";

  TFile f("reco_gen_histos.root","recreate");
  
  posPt_hist_proj->SetName("inclusive_proj_posPtDiff");
  posPt_hist_proj->SetTitle("pos reco - gen pT inclusive");
  posPt_hist_proj->Write();

  negPt_hist_proj->SetName("inclusive_proj_negPtDiff");
  negPt_hist_proj->SetTitle("neg reco - gen pT inclusive");
  negPt_hist_proj->Write();

  multi_hist_proj->SetName("inclusive_proj_pt_diff");
  multi_hist_proj->SetTitle("reco - gen pT inclusive");
  multi_hist_proj->Write();

  posPt_hist->Write();
  negPt_hist->Write();
  posPt_hist2->Write();
  negPt_hist2->Write();
  

  string name, title, title_seed = "reco - gen pT "; 

  TH1F *myhist = new TH1F("pt_diff_4Dbin", "pt diff", 24, -6.0, 6.0);
  
  double events_count = 0.0, entries = 0.0, integral = 0.0, supposed_integral = 0.0;
  int global_bin;
  double content;
  //include under overflow
  for (int pos_eta_bin=0; pos_eta_bin<=nbinseta+1; pos_eta_bin++){
    for (int pos_pt_bin=0; pos_pt_bin<=nbinspt+1; pos_pt_bin++){
      for (int neg_eta_bin=0; neg_eta_bin<=nbinseta+1; neg_eta_bin++){
	for (int neg_pt_bin=0; neg_pt_bin<=nbinspt+1; neg_pt_bin++){

          name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	  title = stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges);
	  //std::cout<<"\n"<< name;

	  /////////// 2 3D histograms method ///////////////////////
	 
	  myhist->SetName(name.c_str());
	  myhist->SetTitle(title.c_str());
	  
	  for (int binz=1; binz<=24; binz++){
	    content = 0.0;
	    content += posPt_hist->GetBinContent( posPt_hist->GetBin(pos_eta_bin,pos_pt_bin,binz) );
	    supposed_integral += posPt_hist2->GetBinContent( posPt_hist2->GetBin(pos_eta_bin,pos_pt_bin) );

	    content += negPt_hist->GetBinContent(negPt_hist->GetBin(neg_eta_bin,neg_pt_bin,binz));
	    myhist->SetBinContent(binz, content);
	    
	    supposed_integral += negPt_hist2->GetBinContent( negPt_hist2->GetBin(neg_eta_bin,neg_pt_bin)  );  
	  }

	  myhist->Write(name.c_str());
	  integral += myhist->Integral(1, 24);

	  ///////// 1 5D histogram method /////////////////////////
          
	  delete gROOT->FindObject("multi_data_frame_proj_4");

	  //For historical reasons, SetRange(0,0) resets the range, SetRange(-1, 0) gets underflow only
	  if(pos_eta_bin != 0){
	    multi_hist->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
	  } else {
	    multi_hist->GetAxis(0)->SetRange(-1, 0);
	  }

	  if(pos_pt_bin != 0){
            multi_hist->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
          } else {
            multi_hist->GetAxis(1)->SetRange(-1, 0);
          }

	  if(neg_eta_bin != 0){
            multi_hist->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
          } else {
            multi_hist->GetAxis(2)->SetRange(-1, 0);
          }
	  
	  if(neg_pt_bin != 0){
            multi_hist->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          } else {
            multi_hist->GetAxis(3)->SetRange(-1, 0);
          }

          multi_hist_proj = multi_hist->Projection(4); 
	  
          entries = multi_hist_proj->Integral(1, 24);
	  events_count =  events_count + entries;
	  //std::cout<<"in histo "<<name<<" there are with SetRange "<< entries << " entries and with SetRangeUser there are ";          
	  /*
	  delete gROOT->FindObject("multi_data_frame_proj_4");
          multi_hist->GetAxis(0)->SetRangeUser(etabinranges[pos_eta_bin - 1]+0.00001, etabinranges[pos_eta_bin]-0.0001);
	  multi_hist->GetAxis(1)->SetRangeUser(ptbinranges[pos_pt_bin - 1]+0.00001, ptbinranges[pos_pt_bin]-0.0001);
          multi_hist->GetAxis(2)->SetRangeUser(etabinranges[neg_eta_bin - 1]+0.00001, etabinranges[neg_eta_bin]-0.0001);
          multi_hist->GetAxis(3)->SetRangeUser(ptbinranges[neg_pt_bin - 1]+0.00001, ptbinranges[neg_pt_bin]-0.0001);
          multi_hist_proj = multi_hist->Projection(4);
	  
	  entries = multi_hist_proj->GetEntries();
	  std::cout<<"entries: "<<entries;
	  events_count =  events_count + entries;
          */	  
          multi_hist_proj->SetName(name.c_str());
          multi_hist_proj->SetTitle(title.c_str());
          multi_hist_proj->Write(name.c_str());
	  
	}
      } 
    }
  } 

  std::cout<<"2 3D histo approach: Initial no. of events: "<<initial_integral<<" sum of events in individual histos: "<<integral<<"; ev diff: "<<initial_integral-integral<<"\n"; 
  std::cout<<"1 5D histo approach: Initial no. of events: "<<initial_entries<<" sum of events in individual histos: "<<events_count<<"; ev diff: "<<initial_entries - events_count<<"\n";
  std::cout<<"Supposed integral: "<<supposed_integral;
  
  f.Close();

  return 0; 

}
