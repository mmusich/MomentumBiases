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
    files=[]

    for root, dirnames, filenames in os.walk("/scratchnvme/wmass/NANOV9/postVFP/SingleMuon"):
        for filename in filenames:
            if '.root' in filename:
                files.append(os.path.join(root, filename))

    for name in files:
        names.push_back(name)
        
    frame = RDataFrame("Events", names)

    # Muon selection
    frame = frame.Filter("HLT_IsoMu24 == 1")
    frame = frame.Filter("nMuon >= 2")
    frame = frame.Filter("PV_npvsGood >= 1")

    #Muon_highPurity

    frame = frame.Define("MuonisGood", "Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_isGlobal && Muon_mediumId && Muon_pfRelIso04_all < 0.15")
    
    # Define invariant mass per mu pairs and filter for opposite sign, MuonisGood and invariant mass range
    frame = frame.Define("pairs","std::vector<std::tuple<int,int,double>> pairs; for(int i=1;i<nMuon;i++){for(int j=0;j<i;j++){if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1){double mll=pow(2*Muon_pt[i]*Muon_pt[j]*(cosh(Muon_eta[i]-Muon_eta[j])-cos(Muon_phi[i]-Muon_phi[j])),0.5); if(60<mll && mll<120){std::tuple<int,int,double> temp; if(Muon_charge[i]==1){temp=make_tuple(i,j,mll);} else{temp=make_tuple(j,i,mll);} pairs.push_back(temp);}}}} if(pairs.size()>1){double diff=100.0; int best=0; for(int i=0;i<pairs.size();i++){if(abs(get<2>(pairs.at(i))-91)<diff){diff=(abs(get<2>(pairs.at(i))-91)); best=i;}} std::vector<std::tuple<int,int,double>> pairss; pairss.push_back(pairs.at(best)); return pairss;} return pairs;")

    #frame = frame.Define("pairs","std::vector<std::tuple<int,int,double>> pairs; for(int i=1;i<nMuon;i++){for(int j=0;j<i;j++){if(MuonisGood[i] && MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1){double mll=pow(2*Muon_pt[i]*Muon_pt[j]*(cosh(Muon_eta[i]-Muon_eta[j])-cos(Muon_phi[i]-Muon_phi[j])),0.5); if(60<mll && mll<120){std::tuple<int,int,double> temp; if(Muon_charge[i]==1){temp=make_tuple(i,j,mll);} else{temp=make_tuple(j,i,mll);} pairs.push_back(temp);}}}} return pairs;")
    frame = frame.Filter("pairs.size()>=1")

    # Define mll branch inclusive in all eta,phi bins
    frame = frame.Define("mll","std::vector<double> temporary; for(int i=0; i<pairs.size(); i++){temporary.push_back(get<2>(pairs.at(i)));} return temporary;")

    # Define momentum
    frame = frame.Define("Muon_p", "std::vector<double> muon_p; for(int j=0;j<nMuon;j++){muon_p.push_back(Muon_pt[j]*std::cosh(Muon_eta[j]));} return muon_p;")
    
    # Binning in eta {-2.4,-1.6,-0.8,0,0.8,1.6,2.4}
    frame = frame.Define("etabin","std::vector<int> etabin; for(int j=0;j<nMuon;j++){etabin.push_back(-99999);} double etabins[7]={-2.4,-1.6,-0.8,0,0.8,1.6,2.4}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(etabins[i]<Muon_eta[k] && Muon_eta[k]<etabins[i+1]){etabin[k]=i+1;}}} return etabin;")
           
    # Binning in phi {-3,-2,-1,0,1,2,3}    
    frame = frame.Define("phibin","std::vector<int> phibin; for(int j=0;j<nMuon;j++){phibin.push_back(-99999);} double phibins[7]={-3,-2,-1,0,1,2,3}; for(int k=0;k<nMuon;k++){for(int i=0;i<6;i++){if(phibins[i]<Muon_phi[k] && Muon_phi[k]<phibins[i+1]){phibin[k]=i+1;}}} return phibin;")

    frame.Snapshot("outputTree", "outputFile.root",{"MuonisGood","etabin","phibin","Muon_pt","mll"})

    # Get average p and pT per eta,phi bin
    p_profile = frame.Profile2D(("p", "p", 6, 1, 7, 6, 1, 7), "etabin", "phibin", "Muon_p")
    hist_p = p_profile.ProjectionXY()
    
    pt_profile = frame.Profile2D(("pt", "pt", 6, 1, 7, 6, 1, 7), "etabin", "phibin", "Muon_pt")
    hist_pt = pt_profile.ProjectionXY()
    
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

    # Book histograms with data from all eta, phi regions
    h1 = frame.Histo1D(ROOT.RDF.TH1DModel("pt", "Pt", 100, 0, 100),"Muon_pt")
    h2 = frame.Histo1D(ROOT.RDF.TH1DModel("eta", "Eta", 60, -3, 3),"Muon_eta")
    h3 = frame.Histo1D(ROOT.RDF.TH1DModel("phi", "Phi", 80, -4, 4),"Muon_phi")
    h4 = frame.Histo1D(ROOT.RDF.TH1DModel("mll", "mll", 60, 60, 120),"mll")
    
    outfilem = ROOT.TFile.Open("other_histos.root","RECREATE")
    outfilem.cd()
    hist_p.Write()
    hist_pt.Write()
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    
    outfilem.Close()

    
frame ()
