import ROOT
from ROOT import RooRealVar, RooConstVar, RooFormulaVar, RooVoigtian, RooExponential
from ROOT import RooAddPdf, RooProdPdf, RooDataSet, RooPlot, RooArgList, RooCategory, RooSimultaneous
from ROOT import gROOT, TLegend, TArrayD, TH2F, gStyle
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
infile = ROOT.TFile.Open("data_histos_binning_6_6_6_6.root","READ")
#infile = ROOT.TFile.Open("data_histos_binning_6_6_3_3.root","READ")

# These must match the binning used to create the histograms
eta_1_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
phi_1_edges = np.array([-3,-2,-1,0,1,2,3])
eta_2_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
phi_2_edges = np.array([-3,-2,-1,0,1,2,3])

#eta_2_edges = np.array([-2.4,-0.8,0.8,2.4])
#phi_2_edges = np.array([-3,-1,1,3])

# Needed for eta, phi distribution of the mean - pdg val and others
m0_eta_phi = TH2F('m0_eta_phi', 'm0-pdg (fill bin method)', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
sigma_eta_phi = TH2F('sigma_eta_phi', 'sigma (fill bin method)', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
width_eta_phi = TH2F('width_eta_phi', 'width-pdg (fill bin method)', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

m0_eta_phi_weight = TH2F('m0_eta_phi_weight', 'm0-pdg (weight method)', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
sigma_eta_phi_weight = TH2F('sigma_eta_phi_weight', 'sigma (weight method)', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
width_eta_phi_weight = TH2F('width_eta_phi_weight', 'width-pdg (weight method)', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

fits_count = TH2F('fits_count', 'fits count', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

m0_weight = {}
sigma_weight = {}
width_weight = {}
m0_count = {}
sigma_count = {}
width_count = {}

# Performance histograms
fit_status = TH2F('fit_status', 'status', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
fit_cov_qual = TH2F('fit_cov_qual', 'covQual', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2 = TH2F('reduced_chi2', 'reduced Chi2', len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

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

for eta1 in range(1,len(eta_1_edges)):
    for phi1 in range(1,len(phi_1_edges)):
        # Initialise values for the eta, phi histograms bin method 
        key_name = f"key_{eta1}_{phi1}"
        m0_weight[key_name] = 0
        sigma_weight[key_name] =  0
        width_weight[key_name] = 0
        m0_count[key_name] = 0
        sigma_count[key_name] = 0
        width_count[key_name] = 0
        
        for eta2 in range(1, len(eta_2_edges)): 
            for phi2 in range(1, len(phi_2_edges)):
                region_name = f"region_{eta1}_{phi1}_{eta2}_{phi2}" 

                # Get data as RooDataHist
                data_histo = infile.Get(f"mll_fast_{eta1}_{phi1}_{eta2}_{phi2}")
                if data_histo: 
                    data_rdh = ROOT.RooDataHist(f"data_rdh[region_{eta1}_{phi1}_{eta2}_{phi2}]", f"data_rdh[region_{eta1}_{phi1}_{eta2}_{phi2}]", [m], Import=data_histo)

                #Fit
                    fit_result = model.fitTo(data_rdh, Save=True)

                #Plot  fits
                    cal = ROOT.TCanvas()
                    mframe = m.frame()
                    data_rdh.plotOn(mframe)
                    model.plotOn(mframe, LineColor="kGreen", LineWidth=2)
                    model.plotOn(mframe, Components = {signal}, LineColor="kBlue", LineWidth=2)
                    model.plotOn(mframe, Components = {bkg}, LineColor="kRed", LineWidth=2)
                    mframe.SetName(f"hist_region_{eta1}_{phi1}_{eta2}_{phi2}")
                    mframe.SetTitle(f"Fit region eta1 in [{eta_1_edges[eta1-1]},{eta_1_edges[eta1]}], phi1 in [{phi_1_edges[phi1-1]},{phi_1_edges[phi1]}], eta2 in [{eta_2_edges[eta2-1]},{eta_2_edges[eta2]}], phi2 in [{phi_2_edges[phi2-1]},{phi_2_edges[phi2]}]")
                    mframe.Draw("same")
                    cal.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fits/sim_fit_region_{eta1}_{phi1}_{eta2}_{phi2}.pdf")
                    cal.Delete()
                    
                # Collect eta, phi distribution of the means - pdg val and others for bin method
                    m0_weight[key_name] = m0_weight[key_name] + m0.getValV()-91.1876
                    sigma_weight[key_name] =  sigma_weight[key_name] + sigma.getValV()
                    width_weight[key_name] = width_weight[key_name] + width.getValV()-2.4952
                    m0_count[key_name] = m0_count[key_name]+1
                    sigma_count[key_name] = sigma_count[key_name]+1
                    width_count[key_name] = width_count[key_name]+1

                # Fill eta, phi histograms weight method
                    m0_eta_phi_weight.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, m0.getValV()-91.1876)
                    sigma_eta_phi_weight.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, sigma.getValV())
                    width_eta_phi_weight.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, width.getValV()-2.4952)

                    fits_count.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001)
                    
                # Fill performance histograms
                    fit_status.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, fit_result.status())
                    fit_cov_qual.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, fit_result.covQual())
                    reduced_chi2.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, mframe.chiSquare(5)) # no of floating parameters in bracket 
                    
        # Fill eta, phi histograms bin method
        if m0_count[key_name]: # this if will be important once histos with low stats are vetoed
            m0_weight[key_name] = m0_weight[key_name]/m0_count[key_name]
            m0_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, m0_weight[key_name])
            sigma_weight[key_name] = sigma_weight[key_name]/sigma_count[key_name]
            sigma_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, sigma_weight[key_name])
            width_weight[key_name] = width_weight[key_name]/width_count[key_name]
            width_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, width_weight[key_name])

# Needed to get averages from the width method
m0_eta_phi_weight.Divide(fits_count)
sigma_eta_phi_weight.Divide(fits_count)
width_eta_phi_weight.Divide(fits_count)

fit_status.Divide(fits_count)
fit_cov_qual.Divide(fits_count)
reduced_chi2.Divide(fits_count)

infile.Close()

# Plot eta, phi histograms
ca = ROOT.TCanvas()
gStyle.SetNumberContours(256)

m0_eta_phi.GetYaxis().SetTitle("#phi")
m0_eta_phi.GetXaxis().SetTitle("#eta")
m0_eta_phi.SetStats(0)
m0_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/m0_eta_phi.pdf")

m0_eta_phi_weight.GetYaxis().SetTitle("#phi")
m0_eta_phi_weight.GetXaxis().SetTitle("#eta")
m0_eta_phi_weight.SetStats(0)
m0_eta_phi_weight.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/m0_eta_phi_weight.pdf")

sigma_eta_phi.GetYaxis().SetTitle("#phi")
sigma_eta_phi.GetXaxis().SetTitle("#eta")
sigma_eta_phi.SetStats(0)
sigma_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/sigma_eta_phi.pdf")

sigma_eta_phi_weight.GetYaxis().SetTitle("#phi")
sigma_eta_phi_weight.GetXaxis().SetTitle("#eta")
sigma_eta_phi_weight.SetStats(0)
sigma_eta_phi_weight.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/sigma_eta_phi_weight.pdf")

width_eta_phi.GetYaxis().SetTitle("#phi")
width_eta_phi.GetXaxis().SetTitle("#eta")
width_eta_phi.SetStats(0)
width_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/width_eta_phi.pdf")

width_eta_phi_weight.GetYaxis().SetTitle("#phi")
width_eta_phi_weight.GetXaxis().SetTitle("#eta")
width_eta_phi_weight.SetStats(0)
width_eta_phi_weight.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/width_eta_phi_weight.pdf")

fits_count.GetYaxis().SetTitle("#phi")
fits_count.GetXaxis().SetTitle("#eta")
fits_count.SetStats(0)
fits_count.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/fits_count.pdf")

# Performance histograms

fit_status.GetYaxis().SetTitle("#phi")
fit_status.GetXaxis().SetTitle("#eta")
fit_status.SetStats(0)
fit_status.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/fit_status.pdf")

fit_cov_qual.GetYaxis().SetTitle("#phi")
fit_cov_qual.GetXaxis().SetTitle("#eta")
fit_cov_qual.SetStats(0)
fit_cov_qual.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/fit_cov_qual.pdf")

reduced_chi2.GetYaxis().SetTitle("#phi")
reduced_chi2.GetXaxis().SetTitle("#eta")
reduced_chi2.SetStats(0)
reduced_chi2.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos/reduced_chi2.pdf")

ca.Delete()

stop_time = time.time()
print("Time elapsed in seconds: ", stop_time-start_time)

# TO DO
# compute reduced chi squares

