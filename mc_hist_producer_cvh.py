# Stand alone code to produce histograms from MC cvh data per eta, phi region

from ROOT import TFile, TLorentzVector, RDataFrame, TH1D, TH2D, TCanvas, TGraph
from ROOT import gROOT
from ROOT import Math, TMath
#TODO maybe replace TLorentzVector with LorentzVector 
import ROOT
from array import array
import numpy as np
import os
from os.path import exists
import time, math

# this is to supress ROOT garbage collector for histograms
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
# this is to stop ROOT from displaying the plot when running
# (better to save it)
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;

//dxy_significance
ROOT::RVecD dxy_significance(ROOT::RVecD Muon_dxy, ROOT::RVecD Muon_dxyErr){
return abs(Muon_dxy)/Muon_dxyErr;
}

// MuonisGood

ROOT::RVecB MuonisGood(ROOT::RVecD Muon_cvhPt, ROOT::RVecD Muon_cvhEta, ROOT::RVecB Muon_isGlobal, ROOT::RVecB Muon_mediumId, ROOT::RVecD Muon_pfRelIso04_all, ROOT::RVecD dxy_significance){
ROOT::RVecB muonisgood;
for(int i=0;i<Muon_cvhPt.size();i++){
     if (Muon_cvhPt[i] > 10 && abs(Muon_cvhEta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
          muonisgood.push_back(1);
     }
}
return muonisgood;
}

//Pairs

ROOT::RVec<std::tuple<int,int,double>> pairs(ROOT::RVecD Muon_cvhPt, ROOT::RVecD Muon_cvhCharge, ROOT::RVecD Muon_cvhEta, ROOT::RVecD Muon_cvhPhi, ROOT::RVecB MuonisGood, ROOT::RVecD Muon_dxy, ROOT::RVecD Muon_dz){

ROOT::RVec<std::tuple<int,int,double>> pairs;
for(int i=1;i<Muon_cvhPt.size();i++){
     for(int j=0;j<i;j++){
          if(MuonisGood[i] && MuonisGood[j] && Muon_cvhCharge[i]*Muon_cvhCharge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
               // Define opening angle
               double gamma_angle = 0;
               gamma_angle = acos(sin(2*atan(exp((-1)*Muon_cvhEta[i])))*sin(2*atan(exp((-1)*Muon_cvhEta[j])))*cos(Muon_cvhPhi[i])*cos(Muon_cvhPhi[j]) + sin(2*atan(exp((-1)*Muon_cvhEta[i])))*sin(2*atan(exp((-1)*Muon_cvhEta[j])))*sin(Muon_cvhPhi[i])*sin(Muon_cvhPhi[j]) + cos(2*atan(exp((-1)*Muon_cvhEta[i])))*cos(2*atan(exp((-1)*Muon_cvhEta[j]))));
               if(gamma_angle > 3.141592653/4){
                    double mll=pow(2*Muon_cvhPt[i]*Muon_cvhPt[j]*(cosh(Muon_cvhEta[i]-Muon_cvhEta[j])-cos(Muon_cvhPhi[i]-Muon_cvhPhi[j])),0.5);
                    if(75<mll && mll<105){
                         std::tuple<int,int,double> temp;
                         if(Muon_cvhCharge[i]==1){
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
     double diff=100.0;
     int best=0;
     for(int i=0;i<pairs.size();i++){
          if(abs(get<2>(pairs.at(i))-91)<diff){
               diff=(abs(get<2>(pairs.at(i))-91));
               best=i;
          }
     }
     ROOT::RVec<std::tuple<int,int,double>> pairss;
     pairss.push_back(pairs.at(best));
     return pairss;
}
return pairs;
}
""")

c = TCanvas()

def frame():

    ROOT.EnableImplicitMT(128) 

    names = ROOT.std.vector('string')()

    files=[]

    for root, dirnames, filenames in os.walk("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV718_NanoProdv1"):
        for filename in filenames:
            if '.root' in filename:
                files.append(os.path.join(root, filename))

    for name in files:
        names.push_back(name)

    # Run over one file for debugging
#    for n in range(2,101):
#        filename = f"/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_{n}.root"
#        names.push_back(filename)

    frame = RDataFrame("Events", names)

    # Muon selection
    frame = frame.Filter("HLT_IsoMu24 == 1")
    frame = frame.Filter("nMuon >= 2")
    frame = frame.Filter("PV_npvsGood >= 1")

    frame = frame.Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)")
    
    frame = frame.Define("MuonisGood", "MuonisGood(Muon_cvhPt, Muon_cvhEta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)")
    
    # Define invariant mass per mu pairs and filter for opposite sign, MuonisGood and invariant mass closest to Z
    frame = frame.Define("pairs", "pairs(Muon_cvhPt, Muon_cvhCharge, Muon_cvhEta, Muon_cvhPhi, MuonisGood, Muon_dxy, Muon_dz)")
    
    frame = frame.Filter("pairs.size()>=1")

    # This works because actually we kept only one pair per event
    # Define branches pt of leading and subleading muon
    frame = frame.Define("leading_pt","double lead_pt; lead_pt=Muon_cvhPt[get<0>(pairs.at(0))]; if(Muon_cvhPt[get<1>(pairs.at(0))]>lead_pt){lead_pt=Muon_cvhPt[get<1>(pairs.at(0))];} return lead_pt;")
    frame = frame.Define("subleading_pt","double sublead_pt; sublead_pt=Muon_cvhPt[get<0>(pairs.at(0))]; if(Muon_cvhPt[get<1>(pairs.at(0))]<sublead_pt){sublead_pt=Muon_cvhPt[get<1>(pairs.at(0))];} return sublead_pt;")

    frame = frame.Define("good_muon_pt","std::vector<double> temporary; temporary.push_back(Muon_cvhPt[get<0>(pairs.at(0))]) ; temporary.push_back(Muon_cvhPt[get<1>(pairs.at(0))]); return temporary;")
    frame = frame.Define("good_muon_eta","std::vector<double> temporary; temporary.push_back(Muon_cvhEta[get<0>(pairs.at(0))]) ; temporary.push_back(Muon_cvhEta[get<1>(pairs.at(0))]); return temporary;")
    frame = frame.Define("good_muon_phi","std::vector<double> temporary; temporary.push_back(Muon_cvhPhi[get<0>(pairs.at(0))]) ; temporary.push_back(Muon_cvhPhi[get<1>(pairs.at(0))]); return temporary;")

    # Define momentum
    frame = frame.Define("good_muon_p", "std::vector<double> muon_p; muon_p.push_back(Muon_cvhPt[get<0>(pairs.at(0))]*std::cosh(Muon_cvhEta[get<0>(pairs.at(0))])); muon_p.push_back(Muon_cvhPt[get<1>(pairs.at(0))]*std::cosh(Muon_cvhEta[get<1>(pairs.at(0))])); return muon_p;")

    
    # Define mll branch inclusive in all eta,phi bins
    frame = frame.Define("mll","std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){temporary.push_back(get<2>(pairs.at(i)));} return temporary;")
    
    # Binning in eta1 {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}    
    frame = frame.Define("etabin1","int etabin_1; double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int i=0;i<6;i++){if(etabins[i]<Muon_cvhEta[get<0>(pairs.at(0))] && Muon_cvhEta[get<0>(pairs.at(0))]<etabins[i+1]){etabin_1=i+1;}} return etabin_1;")
    # Binning in eta2 {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    frame = frame.Define("etabin2","int etabin_2; double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int i=0;i<6;i++){if(etabins[i]<Muon_cvhEta[get<1>(pairs.at(0))] && Muon_cvhEta[get<1>(pairs.at(0))]<etabins[i+1]){etabin_2=i+1;}} return etabin_2;")
    # Binning in phi1 {-3,-2,-1,0,1,2,3}
    frame = frame.Define("phibin1","int phibin_1; double phibins[7]={-M_PI,-2*M_PI/3,-M_PI/3,0,M_PI/3,2*M_PI/3,M_PI}; for(int i=0;i<6;i++){if(phibins[i]<Muon_cvhPhi[get<0>(pairs.at(0))] && Muon_cvhPhi[get<0>(pairs.at(0))]<phibins[i+1]){phibin_1=i+1;}} return phibin_1;")
    # Binning in phi2 {-3,-2,-1,0,1,2,3}
    frame = frame.Define("phibin2","int phibin_2; double phibins[7]={-M_PI,-2*M_PI/3,-M_PI/3,0,M_PI/3,2*M_PI/3,M_PI}; for(int i=0;i<6;i++){if(phibins[i]<Muon_cvhPhi[get<1>(pairs.at(0))] && Muon_cvhPhi[get<1>(pairs.at(0))]<phibins[i+1]){phibin_2=i+1;}} return phibin_2;")  

    frame = frame.Define("posTrackEta","float posTrackEta; posTrackEta=Muon_cvhEta[get<0>(pairs.at(0))]; return posTrackEta;")
    frame = frame.Define("negTrackEta","float negTrackEta; negTrackEta=Muon_cvhEta[get<1>(pairs.at(0))]; return negTrackEta;")
    frame = frame.Define("posTrackPhi","float posTrackPhi; posTrackPhi=Muon_cvhPhi[get<0>(pairs.at(0))]; return posTrackPhi;")
    frame = frame.Define("negTrackPhi","float negTrackPhi; negTrackPhi=Muon_cvhPhi[get<1>(pairs.at(0))]; return negTrackPhi;")
    frame = frame.Define("posTrackPt","float posTrackPt; posTrackPt=Muon_cvhPt[get<0>(pairs.at(0))]; return posTrackPt;")
    frame = frame.Define("negTrackPt","float negTrackPt; negTrackPt=Muon_cvhPt[get<1>(pairs.at(0))]; return negTrackPt;")
    frame = frame.Define("posTrackDz","float posTrackDz; posTrackDz=Muon_dz[get<0>(pairs.at(0))]; return posTrackDz;")
    frame = frame.Define("negTrackDz","float negTrackDz; negTrackDz=Muon_dz[get<1>(pairs.at(0))]; return negTrackDz;")
    frame = frame.Define("posTrackD0","float posTrackD0; posTrackD0=Muon_dxy[get<0>(pairs.at(0))]; return posTrackD0;")
    frame = frame.Define("negTrackD0","float negTrackD0; negTrackD0=Muon_dxy[get<1>(pairs.at(0))]; return negTrackD0;")
    
    # Binning in eta {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    #frame = frame.Define("etabin","std::vector<int> etabin; for(int j=0;j<nMuon;j++){etabin.push_back(-99999);} double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(etabins[i]<Muon_cvhEta[k] && Muon_cvhEta[k]<etabins[i+1]){etabin[k]=i+1;}}} return etabin;")
           
    # Binning in phi {-3,-2,-1,0,1,2,3}    
    #frame = frame.Define("phibin","std::vector<int> phibin; for(int j=0;j<nMuon;j++){phibin.push_back(-99999);} double phibins[7]={-3,-2,-1,0,1,2,3}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(phibins[i]<Muon_cvhPhi[k] && Muon_cvhPhi[k]<phibins[i+1]){phibin[k]=i+1;}}} return phibin;")
    
    frame.Snapshot("outputTree", "mc_outputFile_cvh_shorter_range_medium.root",{"posTrackEta","negTrackEta","posTrackPhi","negTrackPhi","posTrackPt","negTrackPt","posTrackDz","negTrackDz","posTrackD0","negTrackD0"})

    '''
    # Get average p and pT per eta,phi bin    # Remember to change the ranges for a new binning
    p_profile_1 = frame.Profile2D(("p", "p", 6, 1, 7, 6, 1, 7), "etabin1", "phibin1", "Muon_p") 
    hist_p_1 = p_profile_1.ProjectionXY()

    p_profile_2 = frame.Profile2D(("p", "p", 6, 1, 7, 6, 1, 7), "etabin2", "phibin2", "Muon_p")
    hist_p_2 = p_profile_2.ProjectionXY()
    
    pt_profile_1 = frame.Profile2D(("pt", "pt", 6, 1, 7, 6, 1, 7), "etabin1", "phibin1", "Muon_cvhPt")
    hist_pt_1 = pt_profile_1.ProjectionXY()

    pt_profile_2 = frame.Profile2D(("pt", "pt", 6, 1, 7, 6, 1, 7), "etabin2", "phibin2", "Muon_cvhPt")
    hist_pt_2 = pt_profile_2.ProjectionXY()
    '''
    
    # Book histograms per eta, phi as a 5D histogram

    multi_hist = frame.HistoND(("multi_data_frame", "multi_data_frame", 5, (6,6,6,6,30), (1,1,1,1,75), (7,7,7,7,105)), ("etabin1","phibin1","etabin2","phibin2","mll"))
    hist_entries = multi_hist.Projection(4).GetEntries() # to compare entries in the inclusive profile vs the sum of the individual ones
    
    # Profiles
    outfile_fast = ROOT.TFile.Open("mc_histos_cvh_shorter_range_medium.root","RECREATE")
    outfile_fast.cd()
    
    eta_1_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
    phi_1_edges = np.array([-math.pi,-2*math.pi/3,-math.pi/3,0,math.pi/3,2*math.pi/3,math.pi])
    eta_2_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
    phi_2_edges = np.array([-math.pi,-2*math.pi/3,-math.pi/3,0,math.pi/3,2*math.pi/3,math.pi])

    #eta_2_edges = np.array([-2.4,-0.8,0.8,2.4])
    #phi_2_edges = np.array([-3,-1,1,3])
    
    total_entries=0
    
    for etabin1_idx in range(1,len(eta_1_edges)): 
        for phibin1_idx in range(1,len(phi_1_edges)):
            for etabin2_idx  in range(1,len(eta_2_edges)):
                for phibin2_idx in range(1,len(phi_2_edges)):
                    multi_hist.GetAxis(0).SetRange(etabin1_idx,etabin1_idx)
                    multi_hist.GetAxis(1).SetRange(phibin1_idx,phibin1_idx)
                    multi_hist.GetAxis(2).SetRange(etabin2_idx,etabin2_idx)
                    multi_hist.GetAxis(3).SetRange(phibin2_idx,phibin2_idx)
                    multi_hist_proj = multi_hist.Projection(4)
                    multi_hist_proj.SetName(f"mll_fast_{etabin1_idx}_{phibin1_idx}_{etabin2_idx}_{phibin2_idx}") #this name DOES NOT match the slow histograms
                    multi_hist_proj.SetTitle(f"Inv mass eta1 in [{eta_1_edges[etabin1_idx-1]},{eta_1_edges[etabin1_idx]}], phi1 in [{phi_1_edges[phibin1_idx-1]},{phi_1_edges[phibin1_idx]}], eta2 in [{eta_2_edges[etabin2_idx-1]},{eta_2_edges[etabin2_idx]}], phi2 in [{phi_2_edges[phibin2_idx-1]},{phi_2_edges[phibin2_idx]}]")
                    entries = multi_hist_proj.GetEntries()
                    total_entries = total_entries + entries
                    multi_hist_proj.Write()
    outfile_fast.Close()

    print("The profile inclusive in all eta and phi bins has: ", hist_entries, " entries and the sum from the individual histos (fast) is: ", total_entries, "the summed/initial (fast): ", total_entries/hist_entries)

    
    # Book histograms with data from all eta, phi regions
    h1 = frame.Histo1D(ROOT.RDF.TH1DModel("pt", "cvhPt", 100, 0, 100),"good_muon_pt")
    h2 = frame.Histo1D(ROOT.RDF.TH1DModel("eta", "cvhEta", 60, -3, 3),"good_muon_eta")
    h3 = frame.Histo1D(ROOT.RDF.TH1DModel("phi", "cvhPhi", 80, -4, 4),"good_muon_phi")
    h4 = frame.Histo1D(ROOT.RDF.TH1DModel("mll", "mll", 1000, 70, 110),"mll")
    h5 = frame.Histo1D(ROOT.RDF.TH1DModel("leading_pt", "Leading Pt", 100, 0, 100),"leading_pt")
    h6 = frame.Histo1D(ROOT.RDF.TH1DModel("subleading_pt", "Subleading Pt", 100, 0, 100),"subleading_pt")

    outfile_mll = ROOT.TFile.Open("mc_cvh_shorter_range_medium_mll.root","RECREATE")
    outfile_mll.cd()
    h4.Write()
    outfile_mll.Close()
    
    outfilem = ROOT.TFile.Open("mc_other_histos_cvh_shorter_range_medium.root","RECREATE")
    outfilem.cd()
    '''
    hist_p_1.Write()
    hist_pt_1.Write()
    hist_p_2.Write()
    hist_pt_2.Write()
    '''
    multi_hist.Write()
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    h5.Write()
    h6.Write()
    
    outfilem.Close()

start_time = time.time()
frame ()
stop_time = time.time()
print(stop_time-start_time)
