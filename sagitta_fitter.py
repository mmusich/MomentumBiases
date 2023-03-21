import ROOT
from ROOT import RooRealVar, RooConstVar, RooFormulaVar, RooVoigtian, RooExponential
from ROOT import RooAddPdf, RooProdPdf, RooDataSet, RooPlot, RooArgList, RooCategory, RooSimultaneous
from ROOT import gROOT, TLegend, TArrayD, TH2F
import numpy as np
from array import array
import time

# stop ROOT from displaying the plot when running
ROOT.gROOT.SetBatch(True)
# supress ROOT garbage collector for histograms
ROOT.TH1.AddDirectory(False)

start_time = time.time()

# Get data
#infile = ROOT.TFile.Open("data_histos_testfile.root","READ") #small sample for debugging 
infile = ROOT.TFile.Open("data_histos_cuts_profiled.root","READ")

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
        
ca = ROOT.TCanvas()

#Fit variables
m = RooRealVar("m", "m", 60, 120)

# Define data
m0 = RooRealVar("mean", "mean", 91, 80, 100)

# Define remaining variables
width = RooRealVar("width", "width", 2, 0.9, 7)
sigma = RooRealVar("sigma", "sigma", 2, 0.9, 7)

# Define exponential bkg variables
c = RooRealVar("c", "c", -0.00001, -0.1, 0)

# Fraction of signal PDF in model PDF
frac_sig = RooRealVar("frac_sig", "frac_sig", 0.9, 0.7, 0.999)

# Define PDFs
signal = RooVoigtian("signal", "signal", m, m0, width, sigma)
bkg = RooExponential("bkg", "bkg", m, c)

# Construct a signal and background PDF
model = RooAddPdf("model","s+b", signal, bkg, frac_sig)

for eta1 in range(1,7): # range depends on number of eta bins
    for eta2 in range(1, 7):
        for phi1 in range(1,7): # range depends on number of phi bins
            for phi2 in range(1, 7):
                region_name = f"region_{eta1}_{eta2}_{phi1}_{phi2}" 

                # Get data as RooDataHist
                data_histo = infile.Get(f"mll_fast_{eta1}_{eta2}_{phi1}_{phi2}")
                if data_histo: 
                    data_rdh = ROOT.RooDataHist(f"data_rdh[region_{eta1}_{eta2}_{phi1}_{phi2}]", f"data_rdh[region_{eta1}_{eta2}_{phi1}_{phi2}]", [m], Import=data_histo)

                #Fit
                    model.fitTo(data_rdh)

                #Plot fits
                    cal = ROOT.TCanvas()
                    mframe = m.frame()
                    data_rdh.plotOn(mframe)
                    model.plotOn(mframe, LineColor="kGreen", LineWidth=2)
                    model.plotOn(mframe, Components = {signal}, LineColor="kBlue", LineWidth=2)
                    model.plotOn(mframe, Components = {bkg}, LineColor="kRed", LineWidth=2)
                    mframe.SetName(f"hist_region_{eta1}_{eta2}_{phi1}_{phi2}")
                    mframe.SetTitle(f"Fit region eta1={eta1} eta2={eta2} phi1={phi1} phi2={phi2}")
                    mframe.Draw("same")
                    cal.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fits/sim_fit_region_{eta1}_{eta2}_{phi1}_{phi2}.pdf")
                    cal.Delete()
                    
                # Plot eta, phi distribution of the means - pdg val
                    key_name = f"key_{eta1}_{phi1}"
                    m0_weight[key_name] = m0_weight[key_name] + m0.getValV()-91.1876
                    sigma_weight[key_name] =  sigma_weight[key_name] + sigma.getValV()
                    width_weight[key_name] = width_weight[key_name] + width.getValV()-2.4952
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
ca.SaveAs("m0_eta_phi.pdf")

sigma_eta_phi.SetStats(0)
sigma_eta_phi.Draw("COLZ")
ca.SaveAs("sigma_eta_phi.pdf")

width_eta_phi.SetStats(0)
width_eta_phi.Draw("COLZ")
ca.SaveAs("width_eta_phi.pdf")

infile.Close()

stop_time = time.time()
print("Time elapsed in seconds: ", stop_time-start_time)

# TO DO
# compute reduced chi squares

