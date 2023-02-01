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

    ROOT.EnableImplicitMT(32) 

    names = ROOT.std.vector('string')()
    #NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1
    for n in range(2,101):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)
        
    for n in range(102,185):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016FDataPostVFP_TrackFitV718_NanoProdv1/221230_011433/0000/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    #NanoV9Run2016GDataPostVFP_TrackFitV718_NanoProdv1    
    for n in range(1,1000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV718_NanoProdv1/221230_011512/0000/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    for n in range(1000,2000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV718_NanoProdv1/221230_011512/0001/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    for n in range(2000,3000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV718_NanoProdv1/221230_011512/0002/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    for n in range(3000,3449):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV718_NanoProdv1/221230_011512/0003/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    #NanoV9Run2016HDataPostVFP_TrackFitV718_NanoProdv1    
    for n in range(1,1000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV718_NanoProdv1/221230_011551/0000/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    for n in range(1000,2000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV718_NanoProdv1/221230_011551/0001/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    for n in range(2000,3000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV718_NanoProdv1/221230_011551/0002/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)

    for n in range(3000,4000):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV718_NanoProdv1/221230_011551/0003/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)
        
    for n in range(4000,4596):
        filename = f"/scratchnvme/wmass/NANOV9/postVFP/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV718_NanoProdv1/221230_011551/0004/NanoV9DataPostVFP_{n}.root"
        names.push_back(filename)
        
    frame = RDataFrame("Events", names)

    # Muon selection, to be refined later
    frame = frame.Filter("HLT_IsoMu24 == 1")
    frame = frame.Filter("nMuon >= 2")
    
    # tight pID 
    # pT cut

    # Define invariant mass per mu pairs and filter through invariant mass
    frame = frame.Define("pairs","std::vector<std::tuple<int,int,double>> pairs; for(int i=1;i<nMuon;i++){for(int j=0;j<i;j++){if(Muon_charge[i]*Muon_charge[j]==-1){double mll=pow(2*Muon_pt[i]*Muon_pt[j]*(cosh(Muon_eta[i]-Muon_eta[j])-cos(Muon_phi[i]-Muon_phi[j])),0.5); if(60<mll && mll<120){std::tuple<int,int,double> temp(i,j,mll); pairs.push_back(temp);}}}} return pairs;")
    frame = frame.Filter("pairs.size()>=1")
    
    # Binning in eta {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    frame = frame.Define("etabin","std::array<double, 2> etabin; for(int j=0;j<2;j++){etabin[j]=-99999;} double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int k=0;k<2;k++){for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[k] && Muon_eta[k]<etabins[i+1]){etabin[k]=i+1;}}} return etabin;")
           
    # Binning in phi {-3,-2,-1,0,1,2,3}    
    frame = frame.Define("phibin","std::array<double, 2> phibin; for(int j=0;j<2;j++){phibin[j]=-99999;} double phibins[7]={-3,-2,-1,0,1,2,3}; for(int k=0;k<2;k++){for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[k] && Muon_phi[k]<phibins[i+1]){phibin[k]=i+1;}}} return phibin;")
    
    # Book histograms per eta,phi region            

    histo = {}

    for x in range(1,7): # range changes depending on number of bins
        for y in range(1,x+1): # x,y label eta bins
            for m in range(1,7):
                for n in range(1,m+1): # m,n label phi bins
                    branch_name = f"mll_{x}_{y}_{m}_{n}"
                    branch_text = f"std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){{if(((etabin[get<0>(pairs.at(i))]=={x} && etabin[get<1>(pairs.at(i))]=={y}) || (etabin[get<0>(pairs.at(i))]=={y} && etabin[get<1>(pairs.at(i))]=={x})) && ((phibin[get<0>(pairs.at(i))]=={m} && phibin[get<1>(pairs.at(i))]=={n}) || (phibin[get<0>(pairs.at(i))]=={n} && phibin[get<1>(pairs.at(i))]=={m}))){{temporary.push_back(get<2>(pairs.at(i)));}}}} return temporary;"
                    frame = frame.Define(branch_name,branch_text)
                    eta_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4]) 
                    phi_edges = np.array([-3,-2,-1,0,1,2,3])
                    plottitle = f"Inv mass bins eta([{eta_edges[x-1]},{eta_edges[x]}],[{eta_edges[y-1]},{eta_edges[y]}]), phi([{phi_edges[m-1]},{phi_edges[m]}],[{phi_edges[n-1]},{phi_edges[n]}])"
                    var_name = f"h_{x}_{y}_{m}_{n}"
                    histo[var_name] = frame.Histo1D(ROOT.RDF.TH1DModel(f"data_frame_{x}_{y}_{m}_{n}", plottitle, 60, 60, 120),f"mll_{x}_{y}_{m}_{n}")
                    histo[var_name].GetXaxis().SetTitle("Inv mass [GeV]")
                    histo[var_name].GetYaxis().SetTitle("Events")

    # Writing the histograms in a root file
    outfile = ROOT.TFile.Open("data_histos.root","RECREATE")
    outfile.cd()
    for x in range(1,7): # range changes depending on binning
        for y in range(1,x+1):
            for m in range(1,7):
                for n in range(1,m+1):
                    histo[f"h_{x}_{y}_{m}_{n}"].Write()
    outfile.Close()

    # Book histograms with data from all eta, phi regions
    h1 = frame.Histo1D(ROOT.RDF.TH1DModel("pt", "Pt", 100, 0, 100),"Muon_pt")
    h2 = frame.Histo1D(ROOT.RDF.TH1DModel("eta", "Eta", 60, -3, 3),"Muon_eta")
    h3 = frame.Histo1D(ROOT.RDF.TH1DModel("phi", "Phi", 80, -4, 4),"Muon_phi")
     
    outfilem = ROOT.TFile.Open("other_histos.root","RECREATE")
    outfilem.cd()
    h1.Write()
    h2.Write()
    h3.Write()
    
    outfilem.Close()

    
frame ()
