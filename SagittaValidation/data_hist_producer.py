# Stand alone code to produce histograms from data per eta, phi region

from ROOT import TFile, TLorentzVector, RDataFrame, TH1D, TH2D, TCanvas, TGraph
from ROOT import gROOT
from ROOT import Math, TMath
#TODO maybe replace TLorentzVector with LorentzVector 
import ROOT
from array import array
import numpy as np
import os
from os.path import exists
import time
import math

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

ROOT::RVecB MuonisGood(ROOT::RVecD Muon_pt, ROOT::RVecD Muon_eta, ROOT::RVecB Muon_isGlobal, ROOT::RVecB Muon_mediumId, ROOT::RVecD Muon_pfRelIso04_all, ROOT::RVecD dxy_significance){
ROOT::RVecB muonisgood;
for(int i=0;i<Muon_pt.size();i++){
     if (Muon_pt[i] > 10 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
          muonisgood.push_back(1);
     }
}
return muonisgood;
}

//Pairs

ROOT::RVec<std::tuple<int,int,double>> pairs(ROOT::RVecD Muon_pt, ROOT::RVecD Muon_charge, ROOT::RVecD Muon_eta, ROOT::RVecD Muon_phi, ROOT::RVecB MuonisGood, ROOT::RVecD Muon_dxy, ROOT::RVecD Muon_dz, double rest_mass){

ROOT::RVec<std::tuple<int,int,double>> pairs;
std::tuple<int,int,double> temp;

for(int i=1;i<Muon_pt.size();i++){
for(int j=0;j<i;j++){
          if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
               // Define opening angle
               double gamma_angle = 0;
               gamma_angle = acos(sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*cos(Muon_phi[i])*cos(Muon_phi[j]) + sin(2*atan(exp((-1)*Muon_eta[i])))*sin(2*atan(exp((-1)*Muon_eta[j])))*sin(Muon_phi[i])*sin(Muon_phi[j]) + cos(2*atan(exp((-1)*Muon_eta[i])))*cos(2*atan(exp((-1)*Muon_eta[j]))));
               if(gamma_angle > 3.141592653/4){
                    TLorentzVector firstTrack, secondTrack, mother;
                    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);  
                    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
                    mother = firstTrack + secondTrack;
                    double mll = mother.M();
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

    for root, dirnames, filenames in os.walk("/scratchnvme/wmass/NANOV9/postVFP/SingleMuon"):
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
    
    frame = frame.Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)")
    
    # Define invariant mass per mu pairs and filter for opposite sign, MuonisGood and invariant mass closest to Z
    frame = frame.Define("pairs", "pairs(Muon_pt, Muon_charge, Muon_eta, Muon_phi, MuonisGood, Muon_dxy, Muon_dz, 0.105658)") #muMass = 0.105658 GeV
    
    frame = frame.Filter("pairs.size()>=1")

    # This works because actually we kept only one pair per event
    # Define branches pt of leading and subleading muon
    frame = frame.Define("leading_pt","double lead_pt; lead_pt=Muon_pt[get<0>(pairs.at(0))]; if(Muon_pt[get<1>(pairs.at(0))]>lead_pt){lead_pt=Muon_pt[get<1>(pairs.at(0))];} return lead_pt;")
    frame = frame.Define("subleading_pt","double sublead_pt; sublead_pt=Muon_pt[get<0>(pairs.at(0))]; if(Muon_pt[get<1>(pairs.at(0))]<sublead_pt){sublead_pt=Muon_pt[get<1>(pairs.at(0))];} return sublead_pt;")

    frame = frame.Define("good_muon_pt","std::vector<double> temporary; temporary.push_back(Muon_pt[get<0>(pairs.at(0))]) ; temporary.push_back(Muon_pt[get<1>(pairs.at(0))]); return temporary;")
    frame = frame.Define("good_muon_eta","std::vector<double> temporary; temporary.push_back(Muon_eta[get<0>(pairs.at(0))]) ; temporary.push_back(Muon_eta[get<1>(pairs.at(0))]); return temporary;")
    frame = frame.Define("good_muon_phi","std::vector<double> temporary; temporary.push_back(Muon_phi[get<0>(pairs.at(0))]) ; temporary.push_back(Muon_phi[get<1>(pairs.at(0))]); return temporary;")

    frame = frame.Define("posTrackEta","float posTrackEta; posTrackEta=Muon_eta[get<0>(pairs.at(0))]; return posTrackEta;")
    frame = frame.Define("firstPairsIndex", "float first; first = get<0>(pairs.at(0)); return first;")
    frame = frame.Define("secondPairsIndex", "float second; second = get<1>(pairs.at(0)); return second;")
    frame = frame.Define("negTrackEta","float negTrackEta; negTrackEta=Muon_eta[get<1>(pairs.at(0))]; return negTrackEta;")
    frame = frame.Define("posTrackPhi","float posTrackPhi; posTrackPhi=Muon_phi[get<0>(pairs.at(0))]; return posTrackPhi;")
    frame = frame.Define("negTrackPhi","float negTrackPhi; negTrackPhi=Muon_phi[get<1>(pairs.at(0))]; return negTrackPhi;")
    frame = frame.Define("posTrackPt","float posTrackPt; posTrackPt=Muon_pt[get<0>(pairs.at(0))]; return posTrackPt;")
    frame = frame.Define("negTrackPt","float negTrackPt; negTrackPt=Muon_pt[get<1>(pairs.at(0))]; return negTrackPt;")
    frame = frame.Define("posTrackDz","float posTrackDz; posTrackDz=Muon_dz[get<0>(pairs.at(0))]; return posTrackDz;")
    frame = frame.Define("negTrackDz","float negTrackDz; negTrackDz=Muon_dz[get<1>(pairs.at(0))]; return negTrackDz;")
    frame = frame.Define("posTrackD0","float posTrackD0; posTrackD0=Muon_dxy[get<0>(pairs.at(0))]; return posTrackD0;")
    frame = frame.Define("negTrackD0","float negTrackD0; negTrackD0=Muon_dxy[get<1>(pairs.at(0))]; return negTrackD0;")
    
    # Define momentum
    frame = frame.Define("good_muon_p", "std::vector<double> muon_p; muon_p.push_back(Muon_pt[get<0>(pairs.at(0))]*std::cosh(Muon_eta[get<0>(pairs.at(0))])); muon_p.push_back(Muon_pt[get<1>(pairs.at(0))]*std::cosh(Muon_eta[get<1>(pairs.at(0))])); return muon_p;")

    
    # Define mll branch inclusive in all eta,phi bins
    frame = frame.Define("mll","std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){temporary.push_back(get<2>(pairs.at(i)));} return temporary;")
    
    # Binning in eta1 {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}    
    frame = frame.Define("etabin1","int etabin_1; double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[get<0>(pairs.at(0))] && Muon_eta[get<0>(pairs.at(0))]<etabins[i+1]){etabin_1=i+1;}} return etabin_1;")
    # Binning in eta2 {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    frame = frame.Define("etabin2","int etabin_2; double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[get<1>(pairs.at(0))] && Muon_eta[get<1>(pairs.at(0))]<etabins[i+1]){etabin_2=i+1;}} return etabin_2;")
    # Binning in phi1 {-3,-2,-1,0,1,2,3}
    frame = frame.Define("phibin1","int phibin_1; double phibins[7]={-M_PI,-2*M_PI/3,-M_PI/3,0,M_PI/3,2*M_PI/3,M_PI}; for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[get<0>(pairs.at(0))] && Muon_phi[get<0>(pairs.at(0))]<phibins[i+1]){phibin_1=i+1;}} return phibin_1;")
    # Binning in phi2 {-3,-2,-1,0,1,2,3}
    frame = frame.Define("phibin2","int phibin_2; double phibins[7]={-M_PI,-2*M_PI/3,-M_PI/3,0,M_PI/3,2*M_PI/3,M_PI}; for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[get<1>(pairs.at(0))] && Muon_phi[get<1>(pairs.at(0))]<phibins[i+1]){phibin_2=i+1;}} return phibin_2;")  

    # Binning in eta {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    #frame = frame.Define("etabin","std::vector<int> etabin; for(int j=0;j<nMuon;j++){etabin.push_back(-99999);} double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[k] && Muon_eta[k]<etabins[i+1]){etabin[k]=i+1;}}} return etabin;")
           
    # Binning in phi {-3,-2,-1,0,1,2,3}    
    #frame = frame.Define("phibin","std::vector<int> phibin; for(int j=0;j<nMuon;j++){phibin.push_back(-99999);} double phibins[7]={-3,-2,-1,0,1,2,3}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[k] && Muon_phi[k]<phibins[i+1]){phibin[k]=i+1;}}} return phibin;")
    
    frame.Snapshot("outputTree", "outputFile_std_shorter_range_medium.root",{"MuonisGood","Muon_charge","firstPairsIndex","secondPairsIndex","Muon_eta","posTrackEta","negTrackEta","Muon_phi","posTrackPhi","negTrackPhi","Muon_pt","posTrackPt","negTrackPt","posTrackDz","negTrackDz","posTrackD0","negTrackD0"})

    '''
    # Get average p and pT per eta,phi bin    # Remember to change the ranges for a new binning
    p_profile_1 = frame.Profile2D(("p", "p", 6, 1, 7, 6, 1, 7), "etabin1", "phibin1", "Muon_p") 
    hist_p_1 = p_profile_1.ProjectionXY()

    p_profile_2 = frame.Profile2D(("p", "p", 6, 1, 7, 6, 1, 7), "etabin2", "phibin2", "Muon_p")
    hist_p_2 = p_profile_2.ProjectionXY()
    
    pt_profile_1 = frame.Profile2D(("pt", "pt", 6, 1, 7, 6, 1, 7), "etabin1", "phibin1", "Muon_pt")
    hist_pt_1 = pt_profile_1.ProjectionXY()

    pt_profile_2 = frame.Profile2D(("pt", "pt", 6, 1, 7, 6, 1, 7), "etabin2", "phibin2", "Muon_pt")
    hist_pt_2 = pt_profile_2.ProjectionXY()
    '''
    
    # Book histograms per eta, phi as a 5D histogram

    multi_hist = frame.HistoND(("multi_data_frame", "multi_data_frame", 5, (6,6,6,6,30), (1,1,1,1,75), (7,7,7,7,105)), ("etabin1","phibin1","etabin2","phibin2","mll"))
    hist_entries = multi_hist.Projection(4).GetEntries() # to compare entries in the inclusive profile vs the sum of the individual ones

    region_proj = ROOT.TFile.Open("region_proj.root","RECREATE")
    region_proj.cd()

    etabin1_proj = multi_hist.Projection(0)
    etabin1_proj.Write()
    multi_hist.GetAxis(0).SetRange(1,1)
    etabin1_proj = multi_hist.Projection(0)
    etabin1_proj.Write()
    multi_hist.GetAxis(0).SetRange(3,3)
    etabin1_proj = multi_hist.Projection(0)
    etabin1_proj.Write()
    etabin1_proj = multi_hist.Projection(0)
    phibin1_proj = multi_hist.Projection(1)
    etabin2_proj = multi_hist.Projection(2)
    phibin2_proj = multi_hist.Projection(3)
    
    etabin1_proj.Write()
    phibin1_proj.Write()
    etabin2_proj.Write()
    phibin2_proj.Write()
    region_proj.Close()
    
    # Profiles
    outfile_fast = ROOT.TFile.Open("data_histos_std_shorter_range_medium.root","RECREATE")
    outfile_fast.cd()
    
    eta_1_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
    phi_1_edges = np.array([-math.pi,-2*math.pi/3,-math.pi/3,0,math.pi/3,2*math.pi/3,math.pi])
    eta_2_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
    phi_2_edges = np.array([-math.pi,-2*math.pi/3,-math.pi/3,0,math.pi/3,2*math.pi/3,math.pi])

    #eta_2_edges = np.array([-2.4,-0.8,0.8,2.4])
    #phi_2_edges = np.array([-3,-1,1,3])

    open("events_per_region.txt", "w").close() # this deletes the previous content of the file
    f = open('events_per_region.txt', 'a')
    
    total_entries=0
    total_entries_region = 0 
    
    for etabin1_idx in range(1,len(eta_1_edges)): 
        for phibin1_idx in range(1,len(phi_1_edges)):
            total_entries_region = 0
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
                    total_entries_region = total_entries_region + entries
                    multi_hist_proj.Write()
            #count events in eta,phi region of the first muon
            f.write(f"lower eta bin bound:{eta_1_edges[etabin1_idx-1]}, lower phi bin bound:{phi_1_edges[phibin1_idx-1]}, total_entries_region: ")
            f.write("{0:.0f}".format(total_entries_region))
            f.write("\n")
    outfile_fast.Close()

    print("The profile inclusive in all eta and phi bins has: ", hist_entries, " entries and the sum from the individual histos (fast) is: ", total_entries, "the summed/initial (fast): ", total_entries/hist_entries)

    '''
    # The slow way of booking histograms per eta1,eta2,phi1,phi2
    histo = {}

    for x in range(1,7): # range changes depending on number of bins
        for y in range(1,7): # x,y label eta bins
            for m in range(1,7):
                for n in range(1,7): # m,n label phi bins
                    branch_name = f"mll_{x}_{y}_{m}_{n}"
                    branch_text = f"std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){{if((etabin[get<0>(pairs.at(i))]=={x} && etabin[get<1>(pairs.at(i))]=={y}) && (phibin[get<0>(pairs.at(i))]=={m} && phibin[get<1>(pairs.at(i))]=={n})){{temporary.push_back(get<2>(pairs.at(i)));}}}} return temporary;"
                    frame = frame.Define(branch_name,branch_text)
                    eta_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4]) 
                    phi_edges = np.array([-3,-2,-1,0,1,2,3])
                    plottitle = f"Inv mass bins eta([{eta_edges[x-1]},{eta_edges[x]}],[{eta_edges[y-1]},{eta_edges[y]}]), phi([{phi_edges[m-1]},{phi_edges[m]}],[{phi_edges[n-1]},{phi_edges[n]}])"
                    var_name = f"h_{x}_{y}_{m}_{n}"
                    histo[var_name] = frame.Histo1D(ROOT.RDF.TH1DModel(f"data_frame_{x}_{y}_{m}_{n}", plottitle, 60, 60, 120),f"mll_{x}_{y}_{m}_{n}")
                    histo[var_name].GetXaxis().SetTitle("Inv mass [GeV]")
                    histo[var_name].GetYaxis().SetTitle("Events")

    # Writing the histograms in a root file
    outfile = ROOT.TFile.Open("data_histos_cuts.root","RECREATE")
    outfile.cd()

    total_entries_slow=0
    
    for x in range(1,7): # range changes depending on binning
        for y in range(1,7):
            for m in range(1,7):
                for n in range(1,7):
                    entries_slow = histo[f"h_{x}_{y}_{m}_{n}"].GetEntries()
                    total_entries_slow = total_entries_slow + entries_slow
                    histo[f"h_{x}_{y}_{m}_{n}"].Write()
    outfile.Close()

    print("The profile inclusive in all eta and phi bins has: ", hist_entries, " entries and the sum from the individual histos (slow) is: ", total_entries_slow, "the summed/initial (slow): ", total_entries_slow/hist_entries)
    '''
    
    # Book histograms with data from all eta, phi regions
    h1 = frame.Histo1D(ROOT.RDF.TH1DModel("pt", "Pt", 100, 0, 100),"good_muon_pt")
    h2 = frame.Histo1D(ROOT.RDF.TH1DModel("eta", "Eta", 60, -3, 3),"good_muon_eta")
    h3 = frame.Histo1D(ROOT.RDF.TH1DModel("phi", "Phi", 80, -4, 4),"good_muon_phi")
    h4 = frame.Histo1D(ROOT.RDF.TH1DModel("mll", "mll", 40, 70, 110),"mll")
    h5 = frame.Histo1D(ROOT.RDF.TH1DModel("leading_pt", "Leading Pt", 100, 0, 100),"leading_pt")
    h6 = frame.Histo1D(ROOT.RDF.TH1DModel("subleading_pt", "Subleading Pt", 100, 0, 100),"subleading_pt")
    
    outfilem = ROOT.TFile.Open("other_histos_std_shorter_range_medium.root","RECREATE")
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
