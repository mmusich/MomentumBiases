# Code to produce histograms from the data per eta, phi region

from ROOT import TFile, TLorentzVector, RDataFrame, TH1D, TH2D, TCanvas, TGraph
from ROOT import gROOT
from ROOT import Math
#TODO maybe replace TLorentzVector with LorentzVector 
import ROOT
from array import array
import numpy as np

# this is to supress ROOT garbage collector for histograms
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
# this is to stop ROOT from displaying the plot when running
# (better to save it)
ROOT.gROOT.SetBatch(True)

c = TCanvas()

def frame():

    ROOT.EnableImplicitMT() 

    # TODO get data from all files
    frame = RDataFrame("Events","/scratchnvme/wmass/NANOV9/staging/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_2.root")

    # Muon selection, to be refined later
    frame = frame.Filter("HLT_IsoMu24 == 1")
    frame = frame.Filter("nMuon == 2") 
    frame = frame.Filter("Muon_charge[0]*Muon_charge[1] == -1") 

    # tight pID 
    # pT cut

    frame = frame.Define("mll","pow(pow(Muon_pt[0]+Muon_pt[1],2)-pow(Muon_eta[0]+Muon_eta[1],2)-pow(Muon_phi[0]+Muon_phi[1],2)-pow(0.10566+0.10566,2),0.5)") #check muon masss

    # Binning in eta {-2.47,-1.37,-0.4,0,0.4,1.37,2.47}
    frame = frame.Define("etabin1","double etabins[7]={-2.47,-1.37,-0.4,0,0.4,1.37,2.47}; for(int i=0;i<5;i++){if(etabins[i]<Muon_eta[0] && Muon_eta[0]<etabins[i+1]){return i;}} return -99999;")
    frame = frame.Define("etabin2","double etabins[7]={-2.47,-1.37,-0.4,0,0.4,1.37,2.47}; for(int i=0;i<5;i++){if(etabins[i]<Muon_eta[1] && Muon_eta[1]<etabins[i+1]){return i;}} return -99999;")
    
    # Binning in phi {-2.47,-1.37,-0.4,0,0.4,1.37,2.47}    
    frame = frame.Define("phibin1","double phibins[7]={-2.47,-1.37,-0.4,0,0.4,1.37,2.47}; for(int i=0;i<5;i++){if(phibins[i]<Muon_phi[0] && Muon_phi[0]<phibins[i+1]){return i;}} return -99999;") 
    frame = frame.Define("phibin2","double phibins[7]={-2.47,-1.37,-0.4,0,0.4,1.37,2.47}; for(int i=0;i<5;i++){if(phibins[i]<Muon_phi[1] && Muon_phi[1]<phibins[i+1]){return i;}} return -99999;")
    
    # Book histograms per eta,phi region
    histo = {}

    for x in range(1,3): # range changes depending on binning
        for y in range(1,x+1):
            for m in range(1,3):
                for n in range(1,m+1):
                    histtext = f"((etabin1=={x} && etabin2=={y}) || (etabin1=={y} && etabin2=={x})) && ((phibin1=={m} && phibin2=={n}) || (phibin1=={n} && phibin2=={m}))"
                    plottitle = f"mll Bins eta({x},{y}), phi({m},{n})"
                    var_name = f"h_{x}_{y}_{m}_{n}"
                    histo[var_name] = frame.Filter(histtext).Histo1D(ROOT.RDF.TH1DModel(f"data_frame_{x}_{y}_{m}_{n}", plottitle, 100, 0, 500),"mll")
                    histo[var_name].GetXaxis().SetTitle("mll [GeV]")
                    histo[var_name].GetYaxis().SetTitle("Events")

    # Writing the individual histograms in a root file
    outfile = ROOT.TFile.Open("data_histos.root","RECREATE")
    outfile.cd()
    for x in range(1,3): # range changes depending on binning
        for y in range(1,x+1):
            for m in range(1,3):
                for n in range(1,m+1):
                    histo[f"h_{x}_{y}_{m}_{n}"].Write()
    outfile.Close()

    # Book histogram with data from all eta, phi regions
    h = frame.Histo1D(ROOT.RDF.TH1DModel("mll", "Invariant mass", 120, 0, 120),"mll")
    h1 = frame.Histo1D(ROOT.RDF.TH1DModel("pt", "Pt", 100, 0, 90),"Muon_pt")
    h2 = frame.Histo1D(ROOT.RDF.TH1DModel("eta", "Eta", 60, -3, 3),"Muon_eta")
    h3 = frame.Histo1D(ROOT.RDF.TH1DModel("phi", "Phi", 80, -4, 4),"Muon_phi")
   
    outfilem = ROOT.TFile.Open("other_histos.root","RECREATE")
    outfilem.cd()
    h.Write()
    h1.Write()
    h2.Write()
    h3.Write()
    outfilem.Close()

    
frame ()
