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
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
ROOT::RVecI bestPairMuons(ROOT::RVecF GoodMuon_pt, ROOT::RVecF GoodMuon_eta, ROOT::RVecF GoodMuon_phi, ROOT::RVecI GoodMuon_charge) {
  double bestMassdiff = 99999.;
  double bestMass = -9999.;
  unsigned int mu1idx = -1;
  unsigned int mu2idx = -1;
  unsigned int nGoodMuon = GoodMuon_pt.size();
  for(int i=1;i<nGoodMuon;i++){
    for(int j=0;j<i;j++){
      if(GoodMuon_charge[i]*GoodMuon_charge[j]!=-1)    continue; 
      double mll=pow(2*GoodMuon_pt[i]*GoodMuon_pt[j]*(cosh(GoodMuon_eta[i]-GoodMuon_eta[j])-cos(GoodMuon_phi[i]-GoodMuon_phi[j])),0.5);
      if(60. > mll && mll > 120. )  continue;
      double massDiff = abs(mll-91);
      if(massDiff > bestMassdiff)  continue;
      bestMassdiff = massDiff;
      bestMass = mll;
      if(GoodMuon_charge[i] == 1) {
        mu1idx = i;
        mu2idx = j;
      } else {
        mu1idx = j;
        mu2idx = i;
      }
    }
  }
  ROOT::RVecI bPair;
  if(mu1idx != -1 && mu2idx != -1) {
    bPair.push_back(mu1idx);//+ve
    bPair.push_back(mu2idx);//-ve
  }
  return bPair;
}
""")

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
ROOT::RVecF getPair(ROOT::RVecF prop, ROOT::RVecI idx) {
  ROOT::RVecF temp;
  temp.push_back(prop[idx[0]]);
  temp.push_back(prop[idx[1]]);
  return temp;
}
""")

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
ROOT::RVecF getPairPtDiff(ROOT::RVecF reco, ROOT::RVecI gen) {
  ROOT::RVecF temp;
  temp.push_back(reco[0] - gen[0]);
  temp.push_back(reco[1] - gen[1]);
  return temp;
}
""")

ROOT.gInterpreter.Declare("""
// MuonGenMatchedIndex
// works as it is because there is only one muon pair selected  per event
// TODO try a match by deltaR instead
using namespace ROOT::VecOps;
ROOT::RVecF MuonGenMatchedIndex(ROOT::RVecI GenPart_status, ROOT::RVecI GenPart_pdgId, ROOT::RVecI GenPart_genPartIdxMother, ROOT::RVecF GenPart_pt){
  int posindex = 999;
  int negindex = 999;
  ROOT::RVecF genpt;
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
""")


c = TCanvas()

def frame():

    ROOT.EnableImplicitMT(128) 

    names = ROOT.std.vector('string')()
    files=[]

    for root, dirnames, filenames in os.walk("/scratchnvme/wmass/NANOV9/postVFP/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv3/231019_193617/0000/"):
        for filename in filenames:
            if '.root' in filename:
                files.append(os.path.join(root, filename))

    for name in files:
        names.push_back(name)
        
    frame = RDataFrame("Events", names)

    # Muon selection
    frame = frame.Filter("HLT_IsoMu24 == 1", "HLT ")
    frame = frame.Filter("nMuon >= 2", "2 muons in the event")
    frame = frame.Filter("PV_npvsGood >= 1", "good PV")

    #Muon_highPurity
    frame = frame.Define("MuonisGood", "Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_isGlobal && Muon_mediumId && Muon_pfRelIso04_all < 0.15")
    frame = frame.Define("GoodMuon_pt", "Muon_pt[MuonisGood]")
    frame = frame.Define("GoodMuon_eta", "Muon_eta[MuonisGood]")
    frame = frame.Define("GoodMuon_phi", "Muon_phi[MuonisGood]")
    frame = frame.Define("GoodMuon_charge", "Muon_charge[MuonisGood]")
    frame = frame.Filter("GoodMuon_pt.size() >= 2", "event has atleast 2 good muons")    
    # Define invariant mass per mu pairs and filter for opposite sign, MuonisGood and invariant mass range
    frame = frame.Define("bestPairM", "bestPairMuons(GoodMuon_pt, GoodMuon_eta, GoodMuon_phi, GoodMuon_charge)")
    frame = frame.Filter("bestPairM.size() == 2", "event has Mu pair(best pair)")    

    frame = frame.Define("pairRecoPt", "getPair(GoodMuon_pt, bestPairM)")
    frame = frame.Define("pairRecoEta", "getPair(GoodMuon_eta, bestPairM)")

    h2d = frame.Histo2D(("reco_eta_pt", "", 24, -2.4, 2.4, 30, 25., 55.), "pairRecoEta", "pairRecoPt" )

    frame = frame.Define("pairGenPt", "MuonGenMatchedIndex(GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)")
    frame = frame.Define("pairPtDiff", "getPairPtDiff(pairRecoPt, pairGenPt)")

    h2d_eta_ptdiff = frame.Histo2D(("reco_eta_ptDiff", "", 24, -2.4, 2.4, 24, -6., 6.), "pairRecoEta", "pairPtDiff" )
    h2d_pt_ptdiff = frame.Histo2D(("reco_eta_ptDiff", "", 30, 25., 55., 24, -6., 6.), "pairRecoPt", "pairPtDiff" )

    h3d_ptdiff = frame.Histo3D(("reco_eta_pt_ptDiff", "", 24, -2.4, 2.4, 30, 25., 55., 24, -6., 6.), "pairRecoEta", "pairRecoPt", "pairPtDiff" )

    #h5d = frame.HistoND(("name","title", 5, (6, 8, 6, 8, 60), (-2.4, -3.2, -2.4, -3.2, 60.), (2.4, 3.2, 2.4, 3.2, 120.)), ("Mu1_eta", "Mu1_phi", "Mu2_eta","Mu2_phi", "mll"))
    fout=TFile.Open("histo_reco_gen_suv.root", "RECREATE")
    fout.cd()
    h2d.Write()
    h2d_eta_ptdiff.Write()
    h2d_pt_ptdiff.Write()
    h3d_ptdiff.Write()
    fout.Save()
    fout.Close()

    res=frame.Report()
    res.Print()
    
fnp=frame()


