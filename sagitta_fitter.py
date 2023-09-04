import ROOT
from ROOT import RooRealVar, RooConstVar, RooFormulaVar, RooVoigtian, RooExponential, RooFit, RooAbsPdf, RooCrystalBall, RooChebychev, RooGaussian
from ROOT import RooAddPdf, RooProdPdf, RooDataSet, RooPlot, RooArgList, RooCategory, RooSimultaneous, RooCmdArg, RooBreitWigner, RooCBShape, RooFFTConvPdf
from ROOT import gROOT, TLegend, TArrayD, TH2D, gStyle
import numpy as np
from array import array
import time
import sys
import math

# stop ROOT from displaying the plot when running
ROOT.gROOT.SetBatch(True)
# supress ROOT garbage collector for histograms
ROOT.TH1.AddDirectory(False)

start_time = time.time()

# Get data
#infile = ROOT.TFile.Open("data_histos_testfile.root","READ") #small sample for debugging 
#infile = ROOT.TFile.Open("mc_cvh_shorter_range_medium_mll.root","READ") #inclusive sample for debugging
#nametag = ''

#fit_pdf = 'Voigtian'
fit_pdf = 'BWCB' # Breit-Wigner convolved with Crystal-Ball
#fit_pdf = 'DSCB'

# Voigtian
if fit_pdf == 'Voigtian':

    # Get the command line arguments
    args = sys.argv
    # Check if the user provided the correct number of arguments
    if len(args) != 18:
        print("Usage: python sagitta_fitter.py <input_parameter1> <input_parameter2> ... <input_parameter18>")
        sys.exit()

    # Extract the fit parameters
    width_start = float(args[1].rstrip(','))
    width_low = float(args[2].rstrip(','))
    width_high = float(args[3].rstrip(','))

    sigma_start = float(args[4].rstrip(','))
    sigma_low = float(args[5].rstrip(','))
    sigma_high = float(args[6].rstrip(','))

    p1_start = float(args[7].rstrip(','))
    p1_low = float(args[8].rstrip(','))
    p1_high = float(args[9].rstrip(','))

    p2_start = float(args[10].rstrip(','))
    p2_low = float(args[11].rstrip(','))
    p2_high = float(args[12].rstrip(','))

    p3_start = float(args[13].rstrip(','))
    p3_low = float(args[14].rstrip(','))
    p3_high = float(args[15].rstrip(','))

    chi_label = args[16].rstrip(',') # will serve as fit attempt index

    filename = args[17].rstrip("'")

# BWCB
if fit_pdf == 'BWCB':

    # Get the command line arguments
    args = sys.argv
    # Check if the user provided the correct number of arguments
    if len(args) != 15:
        print("Usage: python sagitta_fitter.py <input_parameter1> <input_parameter2> ... <input_parameter15>")
        sys.exit()

    # Extract the fit parameters
    width_start = float(args[1].rstrip(','))
    width_low = float(args[2].rstrip(','))
    width_high = float(args[3].rstrip(','))

    sigma_start = float(args[4].rstrip(','))
    sigma_low = float(args[5].rstrip(','))
    sigma_high = float(args[6].rstrip(','))

    alpha_start = float(args[7].rstrip(','))
    alpha_low = float(args[8].rstrip(','))
    alpha_high = float(args[9].rstrip(','))

    n_start = float(args[10].rstrip(','))
    n_low = float(args[11].rstrip(','))
    n_high = float(args[12].rstrip(','))

    chi_label = args[13].rstrip(',') # will serve as fit attempt index

    filename = args[14].rstrip("'")

# DSCB
if fit_pdf == 'DSCB':

    # Get the command line arguments
    args = sys.argv
    # Check if the user provided the correct number of arguments
    if len(args) != 18:
        print("Usage: python sagitta_fitter.py <input_parameter1> <input_parameter2> ... <input_parameter18>")
        sys.exit()

    # Extract the fit parameters
    sigma_start = float(args[1].rstrip(','))
    sigma_low = float(args[2].rstrip(','))
    sigma_high = float(args[3].rstrip(','))

    alphaL_start = float(args[4].rstrip(','))
    alphaL_low = float(args[5].rstrip(','))
    alphaL_high = float(args[6].rstrip(','))

    nL_start = float(args[7].rstrip(','))
    nL_low = float(args[8].rstrip(','))
    nL_high = float(args[9].rstrip(','))

    alphaR_start = float(args[10].rstrip(','))
    alphaR_low = float(args[11].rstrip(','))
    alphaR_high = float(args[12].rstrip(','))

    nR_start = float(args[13].rstrip(','))
    nR_low = float(args[14].rstrip(','))
    nR_high = float(args[15].rstrip(','))

    chi_label = args[16].rstrip(',') # will serve as fit attempt index

    filename = args[17].rstrip("'")

# Open file
infile = ROOT.TFile.Open(filename,"READ")

filename = filename.replace('.root', '')
print(filename)
if filename[0] == 'm':
    filename = filename.replace('_histos', '')
    print(filename)
if filename[0] == 'd':
    filename = filename.replace('data_histos_', '')

nametag = f"_{filename}_{fit_pdf}"
    
    
# These must match the binning used to create the histograms
eta_1_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
#phi_1_edges = np.array([-math.pi,-2*math.pi/3,-math.pi/3,0,math.pi/3,2*math.pi/3,math.pi])
phi_1_edges = np.array([-3.14, -2.09, -1.05, 0, 1.05, 2.09, 3.14])
eta_2_edges = np.array([-2.4,-1.6,-0.8,0,0.8,1.6,2.4])
#phi_2_edges = np.array([-math.pi,-2*math.pi/3,-math.pi/3,0,math.pi/3,2*math.pi/3,math.pi])
phi_2_edges = np.array([-3.14, -2.09, -1.05, 0, 1.05, 2.09, 3.14])

#eta_2_edges = np.array([-2.4,-0.8,0.8,2.4])
#phi_2_edges = np.array([-3,-1,1,3])

# Eta, phi distribution of the mean - pdg val and others
#m0_eta_phi = TH2D('m0_eta_phi', f"m0-pdg{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
m0_eta_phi = TH2D('m0_eta_phi', "Fitted m_{#mu#mu} mean - PDG Z mass;Positive muon #eta;Positive muon #phi;m_{#mu#mu}-m_{Z} [GeV]", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

#sigma_eta_phi = TH2D('sigma_eta_phi', f"sigma{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
sigma_eta_phi = TH2D('sigma_eta_phi', "Fitted CB sigma ;Positive muon #eta;Positive muon #phi;CB sigma [GeV]", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

width_eta_phi = TH2D('width_eta_phi', f"width-pdg{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

fits_count = TH2D('fits_count', f"fits count{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
events_count = TH2D('events_count', f"events count{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
#events_failed_conv_count = TH2D('events_failed_conv_count', f"Events left after removing non-converged fits{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges)) 

events_failed_conv_count = TH2D('events_failed_conv_count', "Events from converged fits;Positive muon #eta;Positive muon #phi; Events", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

# Performance histograms
#fit_status = TH2D('fit_status', f"Percentage of non-converged fits{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
fit_status = TH2D('fit_status', "Percentage of non-converged fits; Positive muon #eta;Positive muon #phi; %", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

fit_cov_qual_1_3 = TH2D('fit_cov_qual_1_3', f"Percentage of fits with covQual 3{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
fit_cov_qual_m1_0 = TH2D('fit_cov_qual_m1_0', f"Percentage of fits with covQual -1 or 0{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2_average = TH2D('reduced_chi2_average', f"Average reduced Chi2{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
#reduced_chi2_median = TH2D('reduced_chi2_median', f"Median reduced Chi2{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2_median = TH2D('reduced_chi2_median', "Median #chi/ndf;Positive muon #eta;Positive muon #phi;Median #chi/ndf", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

#reduced_chi2_variance = TH2D('reduced_chi2_variance', f"Variance reduced Chi2{nametag}", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))
reduced_chi2_variance = TH2D('reduced_chi2_variance', "Variance #chi/ndf;Positive muon #eta;Positive muon #phi;Variance #chi/ndf", len(eta_1_edges)-1, array('f',eta_1_edges), len(phi_1_edges)-1, array('f',phi_1_edges))

red_chi2_array = np.empty(0)
    
# Fit variables
m = RooRealVar("m", "m", 75, 105)

# Define data
m0 = RooRealVar("mean", "mean", 91.1876, 75, 105)

# Define remaining variables Voigtian
if fit_pdf == 'Voigtian':
    width = RooRealVar("width", "width", width_start, width_low, width_high)
    sigma = RooRealVar("sigma", "sigma", sigma_start, sigma_low, sigma_high)

# Define remaining variables BWCB
if fit_pdf == 'BWCB':
    m0_cb = RooRealVar("mean_cb", "mean_cb", 0., -20., 20.)
    m0_cb.setConstant()
    #width = RooRealVar("width", "width", width_start, width_low, width_high)
    #sigma = RooRealVar("sigma", "sigma", sigma_start, sigma_low, sigma_high)
    #alpha = RooRealVar("alpha", "alpha", alpha_start, alpha_low, alpha_high)
    #n = RooRealVar("n", "n", n_start, n_low, n_high)

    width = RooRealVar("width", "width", 2.4952, 0., 10.)
    width.setConstant()
    sigma = RooRealVar("sigma", "sigma", 1.2, 0.7, 5.)
    alpha = RooRealVar("alpha", "alpha", 1.5, 0.05, 10.)
    n = RooRealVar("n", "n", 1, 0.01, 100.)

    
# Define remaining variables DSCB
if fit_pdf == 'DSCB':
    sigma = RooRealVar("sigma", "sigma", sigma_start, sigma_low, sigma_high) 

    alphaL = RooRealVar("alphaL", "alphaL", alphaL_start, alphaL_low, alphaL_high) 
    nL = RooRealVar("nL", "nL", nL_start, nL_low, nL_high) 
    alphaR = RooRealVar("alphaR", "alphaR", alphaR_start, alphaR_low, alphaR_high) 
    nR = RooRealVar("nR","nR", nR_start, nR_low, nR_high)

# Define exponential bkg variables
c = RooRealVar("c", "c", -1.,   -5.,   5.)

# Define Chebychev background variables 
#p1 = RooConstVar("p1", "p1", -0.35)
#p2 = RooConstVar("p2", "p2", -0.20)  
#p3 = RooConstVar("p3", "p3", -0.009)

#p1 = RooRealVar("p1", "p1", p1_start, p1_low, p1_high)
#p2 = RooRealVar("p2", "p2", p2_start, p2_low, p2_high)
#p3 = RooRealVar("p3", "p3", p3_start, p3_low, p3_high)

# Fraction of signal PDF in model PDF
frac_sig = RooRealVar("frac_sig", "frac_sig", 0.9, 0., 1.)

# Define PDFs
if fit_pdf == 'Voigtian':
    signal = RooVoigtian("signal", "signal", m, m0, width, sigma)

if fit_pdf == 'BWCB':
    bw = RooBreitWigner("bw", "bw", m, m0, width)
    cb = RooCBShape("cb", "cb", m, m0_cb, sigma, alpha, n)
    gaus = RooGaussian("gaus", "gauss", m, m0, sigma)
    signal = RooFFTConvPdf("signal","signal", m, bw, cb)
    
if fit_pdf == 'DSCB':
    signal = RooCrystalBall("signal", "signal", m, m0, sigma, alphaL, nL, alphaR, nR)

# CHANGE NDOF when changing bkg
bkg = RooExponential("bkg", "bkg", m, c)
#bkg = RooChebychev("bkg", "bkg", m, [p1, p2])

# Construct a signal and background PDF
model = RooAddPdf("model","s+b", signal, bkg, frac_sig)

open("getentries.txt", "w").close() # this deletes the previous content of the file

for eta1 in range(1,len(eta_1_edges)):
    for phi1 in range(1,len(phi_1_edges)):
        for eta2 in range(1, len(eta_2_edges)): 
            for phi2 in range(1, len(phi_2_edges)):
                region_name = f"region_{eta1}_{phi1}_{eta2}_{phi2}"
                
                # Get data as RooDataHist
                data_histo = infile.Get(f"mll_fast_{eta1}_{phi1}_{eta2}_{phi2}")
                if data_histo:
                    entries = data_histo.GetEntries()
                    if entries > 100:
                        f = open('getentries.txt', 'a')
                        f.write(f"Fit region eta1 in [{eta_1_edges[eta1-1]},{eta_1_edges[eta1]}], phi1 in [{phi_1_edges[phi1-1]},{phi_1_edges[phi1]}], eta2 in [{eta_2_edges[eta2-1]},{eta_2_edges[eta2]}], phi2 in [{phi_2_edges[phi2-1]},{phi_2_edges[phi2]}]{nametag}")
                        f.write(" {0:.0f}".format(entries))                   
                        f.write("\n")
                
                        data_rdh = ROOT.RooDataHist(f"data_rdh[region_{eta1}_{phi1}_{eta2}_{phi2}]", f"data_rdh[region_{eta1}_{phi1}_{eta2}_{phi2}]", [m], Import=data_histo)

                #Fit
                        fit_result = model.fitTo(data_rdh, Save=True)
                        #Chi fit method
                        #a = ROOT.RooLinkedList()
                        #a.Add(ROOT.RooFit.Save(True))
                        #fit_result = model.chi2FitTo(data_rdh, a)
                        
                #Plot  fits
                        cal = ROOT.TCanvas()
                        mframe = m.frame()

                        data_rdh.plotOn(mframe)
                        model.plotOn(mframe, LineColor="kGreen", LineWidth=2)

                        #CHANGE ndf for Voigtian/DSCB
                        if fit_pdf == 'Voigtian':
                            red_chi2 = mframe.chiSquare(5) # no of floating parameters in bracket
                        if fit_pdf == 'BWCB':
                            red_chi2 = mframe.chiSquare(6) # no of floating parameters in bracket
                        if fit_pdf == 'DSCB':
                            red_chi2 = mframe.chiSquare(8) # no of floating parameters in bracket
                            
                        string = "chi2/ndf {0:.2f}".format(red_chi2)
                        model.paramOn(mframe, RooFit.Layout(0.7,0.9,0.9), RooFit.Label(f"fit status {fit_result.status()}\n cov qual {fit_result.covQual()}\n {string}"))

                        model.plotOn(mframe, Components = {signal}, LineColor="kBlue", LineWidth=2)
                        model.plotOn(mframe, Components = {bkg}, LineColor="kRed", LineWidth=2)
                        
                        mframe.SetName(f"hist_region_{eta1}_{phi1}_{eta2}_{phi2}")
                        mframe.SetTitle(f"Fit region eta1 in [{eta_1_edges[eta1-1]},{eta_1_edges[eta1]}], phi1 in [{phi_1_edges[phi1-1]},{phi_1_edges[phi1]}], eta2 in [{eta_2_edges[eta2-1]},{eta_2_edges[eta2]}], phi2 in [{phi_2_edges[phi2-1]},{phi_2_edges[phi2]}]{nametag}")
                         
                        mframe.Draw("same")
                        cal.SaveAs(f"/home/users/alexe/workingarea/Sagitta/fits{nametag}/sim_fit_region_{eta1}_{phi1}_{eta2}_{phi2}.pdf")
                        cal.Delete()
                        
                # Fill eta, phi histograms
                        if fit_result.status() == 0:
                            m0_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, m0.getValV()-91.1876)
                            sigma_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, sigma.getValV())
                            if fit_pdf == 'Voigtian' or fit_pdf == 'BWCB':
                                width_eta_phi.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, width.getValV()-2.4952)
                        fits_count.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001)
                        events_count.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, entries)
                    
                # Fill performance histograms
                        if fit_result.status() > 0:
                            fit_status.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001)
                            events_failed_conv_count.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001, entries)
                        if fit_result.covQual() == -1 or fit_result.covQual() == 0:
                            fit_cov_qual_m1_0.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001)
                        if fit_result.covQual() == 3:
                            fit_cov_qual_1_3.Fill(eta_1_edges[eta1-1]+0.0001, phi_1_edges[phi1-1]+0.0001)
                            
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

# To track performace of fit attempts
average_red_chi2_median = 0
for eta1 in range(1,len(eta_1_edges)):
    for phi1 in range(1,len(phi_1_edges)):
        average_red_chi2_median += reduced_chi2_median.GetBinContent(eta1, phi1)
average_red_chi2_median = average_red_chi2_median/(len(eta_1_edges)*len(phi_1_edges))

average_fit_status = 0
for eta1 in range(1,len(eta_1_edges)):
    for phi1 in range(1,len(phi_1_edges)):
        average_fit_status += fit_status.GetBinContent(eta1, phi1)
average_fit_status = average_fit_status/(len(eta_1_edges)*len(phi_1_edges))

f = open('attempt_scores.txt', 'a')
f.write(chi_label)
f.write(" {0:.5f}".format(average_red_chi2_median))
f.write("-{0:.5f}".format(average_fit_status))
f.write("\n")

# Needed to get averages in the eta phi histograms
m0_eta_phi.Divide(fits_count)
sigma_eta_phi.Divide(fits_count)
if fit_pdf == 'Voigtian' or fit_pdf == 'BWCB':
    width_eta_phi.Divide(fits_count)

fit_status.Divide(fits_count)
fit_status.Scale(100)
fit_cov_qual_m1_0.Divide(fits_count)
fit_cov_qual_m1_0.Scale(100)
fit_cov_qual_1_3.Divide(fits_count)
fit_cov_qual_1_3.Scale(100)

events_failed_conv_count.Scale(-1)
events_failed_conv_count.Add(events_count)

reduced_chi2_average.Divide(fits_count)

infile.Close()

# Projections



# Plot eta, phi histograms
ca = ROOT.TCanvas()
ca.SetRightMargin(0.20)
gStyle.SetNumberContours(256)

m0_eta_phi.SetMinimum(-0.35)
m0_eta_phi.SetMaximum(0.05)
#m0_eta_phi.GetYaxis().SetTitle("#phi")
#m0_eta_phi.GetXaxis().SetTitle("#eta")
m0_eta_phi.GetZaxis().SetTitleOffset(2.)
m0_eta_phi.SetStats(0)
m0_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/m0_eta_phi{nametag}.pdf")

sigma_eta_phi.SetMinimum(0.85)
sigma_eta_phi.SetMaximum(1.55)
#sigma_eta_phi.GetYaxis().SetTitle("#phi")
#sigma_eta_phi.GetXaxis().SetTitle("#eta")
sigma_eta_phi.SetStats(0)
sigma_eta_phi.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/sigma_eta_phi{nametag}.pdf")

width_eta_phi.SetMinimum(-0.1)
width_eta_phi.SetMaximum(1.5)
width_eta_phi.GetYaxis().SetTitle("#phi")
width_eta_phi.GetXaxis().SetTitle("#eta")
width_eta_phi.SetStats(0)
if fit_pdf == 'Voigtian' or fit_pdf == 'BWCB':
    width_eta_phi.Draw("COLZ")
    ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/width_eta_phi{nametag}.pdf")

fits_count.GetYaxis().SetTitle("#phi")
fits_count.GetXaxis().SetTitle("#eta")
fits_count.SetStats(0)
fits_count.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/fits_count{nametag}.pdf")

if filename[0] == 'm':
    events_count.SetMaximum(1300000)
    events_failed_conv_count.SetMaximum(1300000)
else:
    events_count.SetMaximum(800000)
    events_failed_conv_count.SetMaximum(800000)

events_count.SetMinimum(-0.0001)
events_count.GetYaxis().SetTitle("#phi")
events_count.GetXaxis().SetTitle("#eta")
events_count.SetStats(0)
events_count.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/events_count{nametag}.pdf")

events_failed_conv_count.SetMinimum(-0.0001)
events_failed_conv_count.SetMaximum(810000)
#events_failed_conv_count.GetYaxis().SetTitle("#phi")
#events_failed_conv_count.GetXaxis().SetTitle("#eta")
events_failed_conv_count.SetTitleOffset(2.)
events_failed_conv_count.SetStats(0)
events_failed_conv_count.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/events_failed_conv_count{nametag}.pdf")

# Performance histograms

fit_status.SetMinimum(-0.0001)
#fit_status.SetMaximum(0.45)
#fit_status.GetYaxis().SetTitle("#phi")
#fit_status.GetXaxis().SetTitle("#eta")
fit_status.SetTitleOffset(2.)
fit_status.SetStats(0)
#fit_status.Draw("COLZ text")
fit_status.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/fit_status{nametag}.pdf")

fit_cov_qual_m1_0.GetYaxis().SetTitle("#phi")
fit_cov_qual_m1_0.GetXaxis().SetTitle("#eta")
fit_cov_qual_m1_0.SetStats(0)
fit_cov_qual_m1_0.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/fit_cov_qual_m1_0{nametag}.pdf")

fit_cov_qual_1_3.GetYaxis().SetTitle("#phi")
fit_cov_qual_1_3.GetXaxis().SetTitle("#eta")
fit_cov_qual_1_3.SetStats(0)
fit_cov_qual_1_3.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/fit_cov_qual_1_3{nametag}.pdf")

#reduced_chi2_average.SetMinimum(2)
#reduced_chi2_average.SetMaximum(8.5)
reduced_chi2_average.GetYaxis().SetTitle("#phi")
reduced_chi2_average.GetXaxis().SetTitle("#eta")
reduced_chi2_average.SetStats(0)
reduced_chi2_average.Draw("COLZ text")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/reduced_chi2_average{nametag}.pdf")

#reduced_chi2_median.SetMinimum(1.5)
#reduced_chi2_median.SetMaximum(4.5)
#reduced_chi2_median.GetYaxis().SetTitle("#phi")
#reduced_chi2_median.GetXaxis().SetTitle("#eta")
reduced_chi2_median.SetTitleOffset(2.)
reduced_chi2_median.SetStats(0)
reduced_chi2_median.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/reduced_chi2_median{nametag}.pdf")

#reduced_chi2_variance.SetMinimum(1)
#reduced_chi2_variance.SetMaximum(80)
#reduced_chi2_variance.GetYaxis().SetTitle("#phi")
#reduced_chi2_variance.GetXaxis().SetTitle("#eta")
reduced_chi2_variance.SetTitleOffset(2.)
reduced_chi2_variance.SetStats(0)
reduced_chi2_variance.Draw("COLZ")
ca.SaveAs(f"/home/users/alexe/workingarea/Sagitta/{fit_pdf}/fit_histos{nametag}/reduced_chi2_variance{nametag}.pdf")

ca.Delete()

stop_time = time.time()
print("Time elapsed in seconds: ", stop_time-start_time)

