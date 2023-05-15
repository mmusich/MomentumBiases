import ROOT
from ROOT import RooRealVar, RooConstVar, RooFormulaVar, RooVoigtian, RooExponential, RooFit
from ROOT import RooAddPdf, RooProdPdf, RooDataSet, RooPlot, RooArgList, RooCategory, RooSimultaneous
from ROOT import gROOT, TLegend, TArrayD, TH2D, gStyle
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
#infile = ROOT.TFile.Open("data_histos_binning_6_6_6_6.root","READ")
infile = ROOT.TFile.Open("data_histos_binning_6_6_6_6_nMuon_greater_equal_2_pt_cut_10_tight.root","READ")
#infile = ROOT.TFile.Open("data_histos_binning_6_6_3_3.root","READ")

#nametag = ''
nametag = '_nMuon_greater_equal_2_pt_cut_10_tight'

# These must match the binning used to create the histograms
eta_1_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
phi_1_edges = np.array([-3,-2,-1,0,1,2,3])
eta_2_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
phi_2_edges = np.array([-3,-2,-1,0,1,2,3])

#eta_2_edges = np.array([-2.4,-0.8,0.8,2.4])
#phi_2_edges = np.array([-3,-1,1,3])

# Eta, phi distribution of the mean - pdg val and others
m0_eta_phi = TH2D('m0_eta_phi', f"m0-pdg{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
sigma_eta_phi = TH2D('sigma_eta_phi', f"sigma{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
width_eta_phi = TH2D('width_eta_phi', f"width-pdg{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

fits_count = TH2D('fits_count', f"fits count{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

# Performance histograms
fit_status = TH2D('fit_status', f"status{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
fit_cov_qual = TH2D('fit_cov_qual', f"covQual{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2_average = TH2D('reduced_chi2_average', f"Average reduced Chi2{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2_median = TH2D('reduced_chi2_median', f"Median reduced Chi2{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2_variance = TH2D('reduced_chi2_variance', f"Variance reduced Chi2{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

red_chi2_array = np.empty(0)

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

                    red_chi2 = mframe.chiSquare(5) # no of floating parameters in bracket
                    model.paramOn(mframe, RooFit.Layout(0.7,0.9,0.9), RooFit.Label("chi2/ndf {0:.2f}".format(red_chi2)))
                    
                    model.plotOn(mframe, Components = {signal}, LineColor="kBlue", LineWidth=2)
                    model.plotOn(mframe, Components = {bkg}, LineColor="kRed", LineWidth=2)
                    
                    mframe.SetName(f"hist_region_{eta1}_{phi1}_{eta2}_{phi2}")
                    mframe.SetTitle(f"Fit region eta1 in [{eta_1_edges[eta1-1]},{eta_1_edges[eta1]}], phi1 in [{phi_1_edges[phi1-1]},{phi_1_edges[phi1]}], eta2 in [{eta_2_edges[eta2-1]},{eta_2_edges[eta2]}], phi2 in [{phi_2_edges[phi2-1]},{phi_2_edges[phi2]}]{nametag}")

                    mframe.Draw("same")
                    cal.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fits{nametag}/sim_fit_region_{eta1}_{phi1}_{eta2}_{phi2}.pdf")
                    cal.Delete()
                    
                # Fill eta, phi histograms
                    m0_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, m0.getValV()-91.1876)
                    sigma_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, sigma.getValV())
                    width_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, width.getValV()-2.4952)

                    fits_count.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001)
                    
                # Fill performance histograms
                    fit_status.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, fit_result.status())
                    fit_cov_qual.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, fit_result.covQual())
                    reduced_chi2_average.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, red_chi2)
                    red_chi2_array = np.append(red_chi2_array, red_chi2)

        # Fill the variance reduced chi2 histogram
        average = reduced_chi2_average.GetBinContent(eta1, phi1)/fits_count.GetBinContent(eta1, phi1)
        variance = 0
        
        for idx in range(0,len(red_chi2_array)):
            variance += pow(red_chi2_array[idx]-average, 2)
            
        variance = variance/len(red_chi2_array)
        reduced_chi2_variance.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, variance)
        
        # Fill the median reduced chi2 histogram
        red_chi2_array_sorted = np.sort(red_chi2_array, axis=None)
        if len(red_chi2_array_sorted)%2==1:
            median = red_chi2_array_sorted[len(red_chi2_array_sorted)//2]
        else:
            median = (red_chi2_array_sorted[len(red_chi2_array_sorted)//2-1] + red_chi2_array_sorted[len(red_chi2_array_sorted)//2])/2
        
        reduced_chi2_median.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, median)
        red_chi2_array = np.empty(0)
        
# Needed to get averages in the eta phi histograms
m0_eta_phi.Divide(fits_count)
sigma_eta_phi.Divide(fits_count)
width_eta_phi.Divide(fits_count)

fit_status.Divide(fits_count)
fit_cov_qual.Divide(fits_count)
reduced_chi2_average.Divide(fits_count)

infile.Close()

# Plot eta, phi histograms
ca = ROOT.TCanvas()
gStyle.SetNumberContours(256)

m0_eta_phi.GetYaxis().SetTitle("#phi")
m0_eta_phi.GetXaxis().SetTitle("#eta")
m0_eta_phi.SetStats(0)
m0_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/m0_eta_phi{nametag}.pdf")

sigma_eta_phi.GetYaxis().SetTitle("#phi")
sigma_eta_phi.GetXaxis().SetTitle("#eta")
sigma_eta_phi.SetStats(0)
sigma_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/sigma_eta_phi{nametag}.pdf")

width_eta_phi.GetYaxis().SetTitle("#phi")
width_eta_phi.GetXaxis().SetTitle("#eta")
width_eta_phi.SetStats(0)
width_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/width_eta_phi{nametag}.pdf")

fits_count.GetYaxis().SetTitle("#phi")
fits_count.GetXaxis().SetTitle("#eta")
fits_count.SetStats(0)
fits_count.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/fits_count{nametag}.pdf")

# Performance histograms

fit_status.SetMinimum(-0.0001)
fit_status.GetYaxis().SetTitle("#phi")
fit_status.GetXaxis().SetTitle("#eta")
fit_status.SetStats(0)
fit_status.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/fit_status{nametag}.pdf")

fit_cov_qual.GetYaxis().SetTitle("#phi")
fit_cov_qual.GetXaxis().SetTitle("#eta")
fit_cov_qual.SetStats(0)
fit_cov_qual.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/fit_cov_qual{nametag}.pdf")

reduced_chi2_average.GetYaxis().SetTitle("#phi")
reduced_chi2_average.GetXaxis().SetTitle("#eta")
reduced_chi2_average.SetStats(0)
reduced_chi2_average.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/reduced_chi2_average{nametag}.pdf")

reduced_chi2_median.GetYaxis().SetTitle("#phi")
reduced_chi2_median.GetXaxis().SetTitle("#eta")
reduced_chi2_median.SetStats(0)
reduced_chi2_median.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/reduced_chi2_median{nametag}.pdf")

reduced_chi2_variance.GetYaxis().SetTitle("#phi")
reduced_chi2_variance.GetXaxis().SetTitle("#eta")
reduced_chi2_variance.SetStats(0)
reduced_chi2_variance.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fit_histos{nametag}/reduced_chi2_variance{nametag}.pdf")

ca.Delete()

stop_time = time.time()
print("Time elapsed in seconds: ", stop_time-start_time)

# TO DO
# compute reduced chi squares

