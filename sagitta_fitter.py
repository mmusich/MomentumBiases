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

# Define dictionaries to store data and model
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

# Get data
#infile = ROOT.TFile.Open("data_histos_testfile.root","READ") #small sample for debugging 
infile = ROOT.TFile.Open("data_histos_cuts.root","READ")

# Define PDG value as a fixed variable
m0_pdg = RooConstVar("m0_pdg", "m0_pdg", 91.1876)

# Save plots of the data fits per eta regions in a separate folder
mframes = {}

# Plot eta, phi distribution of the means - pdg val
m0_eta_phi = TH2F('m0_eta_phi', 'm0-pdg', 6, 1, 7, 6, 1, 7)
sigma_eta_phi = TH2F('sigma_eta_phi', 'sigma', 6, 1, 7, 6, 1, 7)
width_eta_phi = TH2F('width_eta_phi', 'width-pdg', 6, 1, 7, 6, 1, 7)

m0_weight = {}
sigma_weight = {}
width_weight = {}
m0_count = {}
sigma_count = {}
width_count = {}

for eta1 in range(1,7):
    for phi1 in range(1,7):
        key_name = f"key_{eta1}_{phi1}"
        m0_weight[key_name] = 0
        sigma_weight[key_name] =  0
        width_weight[key_name] = 0
        m0_count[key_name] = 0
        sigma_count[key_name] = 0
        width_count[key_name] = 0
        
c = ROOT.TCanvas()

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
                if data_histo: 
                    data_rdh[region_name] = ROOT.RooDataHist(f"data_rdh[region_{eta1}_{eta2}_{phi1}_{phi2}]", f"data_rdh[region_{eta1}_{eta2}_{phi1}_{phi2}]", [m], Import=data_histo)

                #Fit
                    models[region_name].fitTo(data_rdh[region_name])

                #Plot fits
                    mframes[region_name] = m.frame(Title=f"Fit region {eta1} {eta2} {phi1} {phi2}")
                    data_rdh[region_name].plotOn(mframes[region_name])
                    models[region_name].plotOn(mframes[region_name], LineColor="kGreen", LineWidth=2)
                    models[region_name].plotOn(mframes[region_name], Components = {signals[region_name]}, LineColor="kBlue", LineWidth=2)
                    models[region_name].plotOn(mframes[region_name], Components = {bkgs[region_name]}, LineColor="kRed", LineWidth=2)

                    mframes[region_name].Draw("same")
                    c.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fits/sim_fit_region_{eta1}_{eta2}_{phi1}_{phi2}.pdf")
                                
                # Plot eta, phi distribution of the means - pdg val
                    key_name = f"key_{eta1}_{phi1}"
                    m0_weight[key_name] = m0_weight[key_name] + m0s[region_name].getValV()-91.1876
                    sigma_weight[key_name] =  sigma_weight[key_name] + sigmas[region_name].getValV()
                    width_weight[key_name] = width_weight[key_name] + widths[region_name].getValV()-2.4952
                    m0_count[key_name] = m0_count[key_name]+1
                    sigma_count[key_name] = sigma_count[key_name]+1
                    width_count[key_name] = width_count[key_name]+1

for eta1 in range(1,7): #this can be joined to the loop above, check most efficient way to read the data
    for phi1 in range(1,7):
        key_name = f"key_{eta1}_{phi1}"
        if m0_count[key_name]:
            m0_weight[key_name] = m0_weight[key_name]/m0_count[key_name]
            m0_eta_phi.Fill(eta1, phi1, m0_weight[key_name])
            sigma_weight[key_name] = sigma_weight[key_name]/sigma_count[key_name]
            sigma_eta_phi.Fill(eta1, phi1, sigma_weight[key_name])
            width_weight[key_name] = width_weight[key_name]/width_count[key_name]
            width_eta_phi.Fill(eta1, phi1, width_weight[key_name])
            print(eta1, " ", phi1, "m0: ", m0_weight[key_name], "sig: ",sigma_weight[key_name], "width: ",width_weight[key_name], " ")
            
m0_eta_phi.SetStats(0)
m0_eta_phi.Draw("COLZ")
c.SaveAs("m0_eta_phi.pdf")

sigma_eta_phi.SetStats(0)
sigma_eta_phi.Draw("COLZ")
c.SaveAs("sigma_eta_phi.pdf")

width_eta_phi.SetStats(0)
width_eta_phi.Draw("COLZ")
c.SaveAs("width_eta_phi.pdf")

infile.Close()

# TO DO
# compute reduced chi squares

