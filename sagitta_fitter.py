import ROOT
from ROOT import RooRealVar, RooConstVar, RooFormulaVar, RooVoigtian, RooExponential
from ROOT import RooAddPdf, RooProdPdf, RooDataSet, RooPlot, RooArgList, RooCategory, RooSimultaneous
from ROOT import gROOT, TLegend, TArrayD, TH2F
import numpy as np
from array import array

# stop ROOT from displaying the plot when running
ROOT.gROOT.SetBatch(True)
# supress ROOT garbage collector for histograms
ROOT.TH1.AddDirectory(False)

m = RooRealVar("m", "m", 60, 120)

# Define sagitta correction delta_sag

#delta_sags = {} # empty dictionary
#for eta in range(1,7): # range depends on number of bins
#    for phi in range(1,7):
#        delta_sags[f"({eta},{phi})"] = RooRealVar(f"delta_sag_{eta}_{phi}", f"delta_sag_{eta}_{phi}",  0, -0.1, 0.1)

# Define average muon transverse momentum per eta bin
#infile_pt = ROOT.TFile.Open("other_histos.root","READ")
#average_pt_histo = infile_pt.Get("pt_pxy")

#average_pts = {} # empty dictionary
#for eta in range(1,7): # range depends on number of bins
#    for phi in range(1,7):
#        average_pts[f"({eta},{phi})"] = RooConstVar(f"average_pt_{eta}_{phi}", f"average_pt_{eta}_{phi}",  average_pt_histo.GetBinContent(eta,phi)) 
       #print('\n eta bin:', eta, 'phi bin:', phi, 'val pt:', average_pt_histo.GetBinContent(eta,phi))

#infile_pt.Close()

# Define dictionaries to store data and model, to be used in data+MC fit
models = {}
signals = {}
bkgs = {} 
data_rdh = {} 

# Define dictionaries for fit variables in each region
# Signal Voigtian variables
m0s = {}
widths = {}
sigmas = {}

# Background variables
cs = {}

# Fraction of signal to bkg in PDFs
fracs_sig = {}

# Define quantities for data+MC fit
regions = ROOT.RooCategory("regions", "regions")
simPdf = ROOT.RooSimultaneous("simPdf","simultaneous pdf", regions)

# Data+MC fit
#infile = ROOT.TFile.Open("data_histos_testfile.root","READ")
infile = ROOT.TFile.Open("data_histos_cuts.root","READ")

# Define PDG value as a fixed variable
m0_pdg = RooConstVar("m0_pdg", "m0_pdg", 91.1876)

for eta1 in range(1,7): # range depends on number of eta bins
    for eta2 in range(1, 7):
        for phi1 in range(1,7): # range depends on number of phi bins
            for phi2 in range(1, 7):
                region_name = f"region_{eta1}_{eta2}_{phi1}_{phi2}" 
                # Define data, get eqn and check rooarglist for constant m_mc
                m0s[region_name] = RooRealVar(f"mean_{eta1}_{eta2}_{phi1}_{phi2}", f"mean_{eta1}_{eta2}_{phi1}_{phi2}", 91, 80, 100) 
         
                # Define remaining variables
                widths[region_name] = RooRealVar(f"width_{eta1}_{eta2}_{phi1}_{phi2}", f"width_{eta1}_{eta2}_{phi1}_{phi2}", 2, 0.9, 7)
                sigmas[region_name] = RooRealVar(f"sigma_{eta1}_{eta2}_{phi1}_{phi2}", f"sigma_{eta1}_{eta2}_{phi1}_{phi2}", 2, 0.9, 7)

                # Define exponential bkg variables
                cs[region_name] = RooRealVar(f"c_{eta1}_{eta2}_{phi1}_{phi2}", f"c_{eta1}_{eta2}_{phi1}_{phi2}", -0.00001, -0.1, 0)
 
                # Fraction of signal PDF in model PDF
                fracs_sig[region_name] = RooRealVar(f"frac_sig_{eta1}_{eta2}_{phi1}_{phi2}", f"frac_sig_{eta1}_{eta2}_{phi1}_{phi2}", 0.9, 0.7, 0.999)    

                # Define PDFs
                signals[region_name] = RooVoigtian("signal", "signal", m, m0s[region_name], widths[region_name], sigmas[region_name])
                bkgs[region_name] = RooExponential("bkg", "bkg", m, cs[region_name])

                # Construct a signal and background PDF
                models[region_name] = RooAddPdf(f"model_{eta1}_{eta2}_{phi1}_{phi2}",f"s+b_{eta1}_{eta2}_{phi1}_{phi2}", signals[region_name], bkgs[region_name], fracs_sig[region_name]) 

                # Get data as RooDataHist
                data_histo = infile.Get(f"data_frame_{eta1}_{eta2}_{phi1}_{phi2}")
                #if data_histo: 
                #if not data_histo:
                #    raise RuntimeError(f"Object data_frame_{eta1}_{eta2}_{phi1}_{phi2} does not exist!")

                data_rdh[region_name] = ROOT.RooDataHist(f"data_rdh[region_{eta1}_{eta2}_{phi1}_{phi2}]", f"data_rdh[region_{eta1}_{eta2}_{phi1}_{phi2}]", [m], Import=data_histo)
        
                # Define new region
                regions.defineType(region_name)
                # Add model to simultaneous PDF
                simPdf.addPdf(models[region_name], region_name)

# Combine data histograms
combHist = ROOT.RooDataHist("combHist","combined hist", RooArgList(m), regions, data_rdh)

# Perform fit
simPdf.fitTo(combHist)

# Plot eta, phi distribution of the means - pdg val
m0_eta_phi = TH2F('m0_eta_phi', 'm0-pdg', 6, 1, 7, 6, 1, 7)
sigma_eta_phi = TH2F('sigma_eta_phi', 'sigma', 6, 1, 7, 6, 1, 7)
width_eta_phi = TH2F('width_eta_phi', 'width-pdg', 6, 1, 7, 6, 1, 7)
for i in range(1,7):
    for k in range(1,7):
        m0_weight = 0
        sigma_weight = 0
        width_weight = 0
        m0_count = 0
        sigma_count = 0
        width_count = 0
        for j in range(1,7):
            for n in range(1,7):
                #data_histo = infile.Get(f"data_frame_{i}_{j}_{k}_{n}")
                #if data_histo:
                region_name = f"region_{i}_{j}_{k}_{n}"
                m0_weight = m0_weight + m0s[region_name].getValV()-91.1876
                sigma_weight =  sigma_weight + sigmas[region_name].getValV()
                width_weight = width_weight + widths[region_name].getValV()-2.4952
                m0_count = m0_count+1
                sigma_count = sigma_count+1
                width_count = width_count+1
        if m0_count:
            m0_weight = m0_weight/m0_count
            m0_eta_phi.Fill(i, k, m0_weight)
            sigma_weight = sigma_weight/sigma_count
            sigma_eta_phi.Fill(i, k, sigma_weight)
            width_weight = width_weight/width_count
            width_eta_phi.Fill(i, k, width_weight)
            print(i, " ", k, "m0: ", m0_weight, "sig: ",sigma_weight, "width: ",width_weight, " ")
            
c = ROOT.TCanvas()
m0_eta_phi.SetStats(0)
m0_eta_phi.Draw("COLZ")
c.SaveAs("m0_eta_phi.pdf")

sigma_eta_phi.SetStats(0)
sigma_eta_phi.Draw("COLZ")
c.SaveAs("sigma_eta_phi.pdf")

width_eta_phi.SetStats(0)
width_eta_phi.Draw("COLZ")
c.SaveAs("width_eta_phi.pdf")

# Save plots of the data fits per eta regions in a separate folder
mframes = {}

for eta1 in range(1,7): # range depends on number of eta bins
    for eta2 in range(1, 7):
        for phi1 in range(1,7): # range depends on number of phi bins
            for phi2 in range(1, 7):
                #data_histo = infile.Get(f"data_frame_{eta1}_{eta2}_{phi1}_{phi2}")
                #if data_histo:
                region_name = f"region_{eta1}_{eta2}_{phi1}_{phi2}"        
                mframes[region_name] = m.frame(Title=f"Fit region {eta1} {eta2} {phi1} {phi2}")
                data_rdh[region_name].plotOn(mframes[region_name])
                models[region_name].plotOn(mframes[region_name], LineColor="kGreen", LineWidth=2)
                models[region_name].plotOn(mframes[region_name], Components = {signals[region_name]}, LineColor="kBlue", LineWidth=2)
                models[region_name].plotOn(mframes[region_name], Components = {bkgs[region_name]}, LineColor="kRed", LineWidth=2)
            
                #c = ROOT.TCanvas()
                mframes[region_name].Draw("same")
                c.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fits/sim_fit_region_{eta1}_{eta2}_{phi1}_{phi2}.pdf")

infile.Close()

# TO DO
# compute reduced chi squares

