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

# this is to supress ROOT garbage collector for histograms
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
# this is to stop ROOT from displaying the plot when running
# (better to save it)
ROOT.gROOT.SetBatch(True)

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

    frame = frame.Define("MuonisGood", "Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_isGlobal && Muon_mediumId && Muon_pfRelIso04_all < 0.15")
    
    # Define invariant mass per mu pairs and filter for opposite sign, MuonisGood and invariant mass closest to Z
    frame = frame.Define("pairs","std::vector<std::tuple<int,int,double>> pairs; for(int i=1;i<nMuon;i++){for(int j=0;j<i;j++){if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1){double mll=pow(2*Muon_pt[i]*Muon_pt[j]*(cosh(Muon_eta[i]-Muon_eta[j])-cos(Muon_phi[i]-Muon_phi[j])),0.5); if(60<mll && mll<120){std::tuple<int,int,double> temp; if(Muon_charge[i]==1){temp=make_tuple(i,j,mll);} else{temp=make_tuple(j,i,mll);} pairs.push_back(temp);}}}} if(pairs.size()>1){double diff=100.0; int best=0; for(int i=0;i<pairs.size();i++){if(abs(get<2>(pairs.at(i))-91)<diff){diff=(abs(get<2>(pairs.at(i))-91)); best=i;}} std::vector<std::tuple<int,int,double>> pairss; pairss.push_back(pairs.at(best)); return pairss;} return pairs;")

    frame = frame.Filter("pairs.size()>=1")

    # Define mll branch inclusive in all eta,phi bins
    frame = frame.Define("mll","std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){temporary.push_back(get<2>(pairs.at(i)));} return temporary;")

    # Define momentum
    frame = frame.Define("Muon_p", "std::vector<double> muon_p; for(int j=0;j<nMuon;j++){muon_p.push_back(Muon_pt[j]*std::cosh(Muon_eta[j]));} return muon_p;")
    
    # Binning in eta {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    frame = frame.Define("etabin","std::vector<int> etabin; for(int j=0;j<nMuon;j++){etabin.push_back(-99999);} double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[k] && Muon_eta[k]<etabins[i+1]){etabin[k]=i+1;}}} return etabin;")
           
    # Binning in phi {-3,-2,-1,0,1,2,3}    
    frame = frame.Define("phibin","std::vector<int> phibin; for(int j=0;j<nMuon;j++){phibin.push_back(-99999);} double phibins[7]={-3,-2,-1,0,1,2,3}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[k] && Muon_phi[k]<phibins[i+1]){phibin[k]=i+1;}}} return phibin;")

   # frame = frame.Define("etabin1","int etabin_1; double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int i=0;i<8;i++){if(etabins[i]<Muon_eta[0] && Muon_eta[0]<etabins[i+1]){etabin_1=i+1;}} return etabin_1;")
  #  frame = frame.Define("etabin2","int etabin_2; double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[1] && Muon_eta[1]<etabins[i+1]){etabin_2=i+1;}} return etabin_2;")
  #  frame = frame.Define("phibin1","int phibin_1; double phibins[7]={-3,-2,-1,0,1,2,3}; for(int i=0;i<8;i++){if(phibins[i]<Muon_phi[0] && Muon_phi[0]<phibins[i+1]){phibin_1=i+1;}} return phibin_1;")
  #  frame = frame.Define("phibin2","int phibin_2; double phibins[7]={-3,-2,-1,0,1,2,3}; for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[1] && Muon_phi[1]<phibins[i+1]){phibin_2=i+1;}} return phibin_2;")  

    frame = frame.Define("etabin1","int etabin_1; double etabins[3]={-2.4,0,2.4}; for(int i=0;i<8;i++){if(etabins[i]<Muon_eta[0] && Muon_eta[0]<etabins[i+1]){etabin_1=i+1;}} return etabin_1;")
    frame = frame.Define("etabin2","int etabin_2; double etabins[3]={-2.4,0,2.4}; for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[1] && Muon_eta[1]<etabins[i+1]){etabin_2=i+1;}} return etabin_2;")
    frame = frame.Define("phibin1","int phibin_1; double phibins[3]={-3,0,3}; for(int i=0;i<8;i++){if(phibins[i]<Muon_phi[0] && Muon_phi[0]<phibins[i+1]){phibin_1=i+1;}} return phibin_1;")
    frame = frame.Define("phibin2","int phibin_2; double phibins[3]={-3,0,3}; for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[1] && Muon_phi[1]<phibins[i+1]){phibin_2=i+1;}} return phibin_2;")
    
    frame.Snapshot("outputTree", "outputFile.root",{"MuonisGood","etabin","phibin","Muon_pt","mll"})

    # Get average p and pT per eta,phi bin
    p_profile = frame.Profile2D(("p", "p", 6, 1, 7, 6, 1, 7), "etabin", "phibin", "Muon_p")
    hist_p = p_profile.ProjectionXY()
    
    pt_profile = frame.Profile2D(("pt", "pt", 6, 1, 7, 6, 1, 7), "etabin", "phibin", "Muon_pt")
    hist_pt = pt_profile.ProjectionXY()

    # Book histograms per eta, phi as a 5D histogram

    #multi_hist = frame.HistoND(("multi_data_frame", "multi_data_frame", 5, (6,6,6,6,60), (1,1,1,1,60), (7,7,7,7,120)), ("etabin1","phibin1","etabin2","phibin2","mll"))
    multi_hist = frame.HistoND(("multi_data_frame", "multi_data_frame", 5, (2,2,2,2,60), (1,1,1,1,60), (3,3,3,3,120)), ("etabin1","phibin1","etabin2","phibin2","mll"))
    hist_entries = multi_hist.Projection(4).GetEntries() # to compare entries in the inclusive profile vs the sum of the individual ones
    
    # Profiles
    outfile = ROOT.TFile.Open("data_histos_cuts_profiled.root","RECREATE")
    outfile.cd()

    
   # eta_1_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
   # phi_1_edges = np.array([-3,-2,-1,0,1,2,3])
   # eta_2_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
   # phi_2_edges = np.array([-3,-2,-1,0,1,2,3])
    
    eta_1_edges = np.array([-2.4,0,2.4])
    phi_1_edges = np.array([-3,0,3])
    eta_2_edges = np.array([-2.4,0,2.4])
    phi_2_edges = np.array([-3,0,3])
    total_entries=0
    
    for etabin1_idx in range(1,3): # ranges change with number of bins 
        for phibin1_idx in range(1,2): 
            for etabin2_idx  in range(1,2):
                for phibin2_idx in range(1,2): 
                    multi_hist.GetAxis(0).SetRange(etabin1_idx,etabin1_idx+1)
                    multi_hist.GetAxis(1).SetRange(phibin1_idx,phibin1_idx+1)
                    multi_hist.GetAxis(2).SetRange(etabin2_idx,etabin2_idx+1)
                    multi_hist.GetAxis(3).SetRange(phibin2_idx,phibin2_idx+1)
                    multi_hist_proj = multi_hist.Projection(4)
                    entries = multi_hist_proj.GetEntries()
                    total_entries = total_entries + entries
                    multi_hist_proj.SetName(f"mll_{etabin1_idx}_{etabin2_idx}_{phibin1_idx}_{phibin2_idx}") #this name matches the histograms from the other method
                    multi_hist_proj.SetTitle(f"Inv mass bins eta([{eta_1_edges[etabin1_idx-1]},{eta_1_edges[etabin1_idx]}],[{eta_2_edges[etabin2_idx-1]},{eta_2_edges[etabin2_idx]}]), phi([{phi_1_edges[phibin1_idx-1]},{phi_1_edges[phibin1_idx]}],[{phi_2_edges[phibin2_idx-1]},{phi_2_edges[phibin2_idx]}])")
                    multi_hist_proj.Write()
    outfile.Close()

    print("The profile inclusive in all eta and phi bins has: ", hist_entries, " entries and the sum from the individual histos is: ", total_entries, "the summed/initial: ", total_entries/hist_entries)
    
    '''
    # The slow way of booking histograms per eta1,eta2,phi1,phi2
    # Book histograms per eta,phi region            

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
    for x in range(1,7): # range changes depending on binning
        for y in range(1,7):
            for m in range(1,7):
                for n in range(1,7):
                    histo[f"h_{x}_{y}_{m}_{n}"].Write()
    outfile.Close()
    '''

    
    # Book histograms with data from all eta, phi regions
    h1 = frame.Histo1D(ROOT.RDF.TH1DModel("pt", "Pt", 100, 0, 100),"Muon_pt")
    h2 = frame.Histo1D(ROOT.RDF.TH1DModel("eta", "Eta", 60, -3, 3),"Muon_eta")
    h3 = frame.Histo1D(ROOT.RDF.TH1DModel("phi", "Phi", 80, -4, 4),"Muon_phi")
    h4 = frame.Histo1D(ROOT.RDF.TH1DModel("mll", "mll", 60, 60, 120),"mll")
    
    outfilem = ROOT.TFile.Open("other_histos.root","RECREATE")
    outfilem.cd()
    hist_p.Write()
    hist_pt.Write()
    multi_hist.Write()
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    
    outfilem.Close()

start_time = time.time()
frame ()
stop_time = time.time()
print(stop_time-start_time)
