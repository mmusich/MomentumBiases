// Standalone code to fit for sagitta

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

using namespace ROOT;
using namespace ROOT::VecOps;

// MuonGenMatchedIndex
// works as it is because there is only one muon pair selected  per event
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

string stringify(int a, int b, int c,int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

int frame(){

  TFile *f1 = new TFile("snapshot_output.root","RECREATE");
  TFile *f2 = new TFile("reco_gen_histos.root","RECREATE");
  
  TChain chain("Events");
  //chain.Add("/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_5.root");
  //chain.Add("/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_6.root");
  chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_773.root");
  chain.Add("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/NanoV9MCPostVFP_772.root");
  
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

  auto d5 = d4.Define("ptDiff","RVecF ptDiff; ptDiff.push_back(recoPt[0] - genPt[0]); ptDiff.push_back(recoPt[1] - genPt[1]); return ptDiff;");

  //d5.Snapshot("Events", "snapshot_output.root", {"Muon_pt", "Muon_charge", "recoPt", "genPt", "ptDiff", "GenPart_pdgId", "GenPart_pt"});
  
  //auto multi_hist = d5.HistoND<float, float, float, float, RVecF, float>({"multi_data_frame", "multi_data_frame", 5, {16,20,16,20,12}, {-2.4,10,-2.4,10,-6}, {2.4,90,2.4,90,6}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff","genWeight"});

  //larger bins while debugging
  auto multi_hist = d5.HistoND<float, float, float, float, RVecF, float>({"multi_data_frame", "multi_data_frame", 5, {2,4,2,4,12}, {-2.4,10,-2.4,10,-6}, {2.4,90,2.4,90,6}}, {"posTrackEta","posTrackPt","negTrackEta","negTrackPt","ptDiff","genWeight"});
  
  std::cout << "There are " << multi_hist->Projection(4)->GetEntries() << " entries in the inclusive projection \n";
  
  auto multi_hist_proj = multi_hist->Projection(4);
  multi_hist_proj->SetName("inclusive_proj_pt_diff");
  multi_hist_proj->SetTitle("reco - gen pT inclusive");
  multi_hist_proj->Write("reco_gen_histos.root");
  
  multi_hist_proj = multi_hist->Projection(0);
  multi_hist_proj->SetName("inclusive_proj_pos_eta");
  multi_hist_proj->SetTitle("pos eta inclusive");
  multi_hist_proj->Write("reco_gen_histos.root");

  multi_hist_proj = multi_hist->Projection(1);
  multi_hist_proj->SetName("inclusive_proj_pos_pt");
  multi_hist_proj->SetTitle("pos pt inclusive");
  multi_hist_proj->Write("reco_gen_histos.root");

  multi_hist_proj = multi_hist->Projection(2);
  multi_hist_proj->SetName("inclusive_proj_neg_eta");
  multi_hist_proj->SetTitle("neg eta inclusive");
  multi_hist_proj->Write("reco_gen_histos.root");

  multi_hist_proj = multi_hist->Projection(3);
  multi_hist_proj->SetName("inclusive_proj_neg_pt");
  multi_hist_proj->SetTitle("neg pt inclusive");
  multi_hist_proj->Write("reco_gen_histos.root");

  string name, title, title_seed = "reco - gen pT "; 

  for (int pos_eta_bin=1; pos_eta_bin<=2; pos_eta_bin++){
    for (int pos_pt_bin=1; pos_pt_bin<=4; pos_pt_bin++){
      for (int neg_eta_bin=1; neg_eta_bin<=2; neg_eta_bin++){
	for (int neg_pt_bin=1; neg_pt_bin<=4; neg_pt_bin++){
          multi_hist->GetAxis(0)->SetRange(pos_eta_bin, pos_eta_bin);
          multi_hist->GetAxis(1)->SetRange(pos_pt_bin, pos_pt_bin);
          multi_hist->GetAxis(2)->SetRange(neg_eta_bin, neg_eta_bin);
          multi_hist->GetAxis(3)->SetRange(neg_pt_bin, neg_pt_bin);
          multi_hist_proj = multi_hist->Projection(4);
	  
	  name = stringify(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
          title = title_seed + name;

          multi_hist_proj->SetName(name.c_str());
          multi_hist_proj->SetTitle(title.c_str());
          multi_hist_proj->Write("reco_gen_histos.root");
	}
      } 
    }
  }
  
  return 0; 

}
