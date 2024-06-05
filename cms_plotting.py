# ROOT imports 
import os, ROOT
import cmsstyle as CMS

# File reading 
f_constants_fitted = ROOT.TFile.Open('InOutputFiles/constants_fitted_for_plotting.root')
dummy_A = f_constants_fitted.Get("dummy_pars_A")
fitted_A = f_constants_fitted.Get("fitted_pars_A")

# Styling
CMS.SetExtraText("Simulation Private Work")
canv_name = 'A'
iPos = 11 
CMS.SetLumi("")
CMS.SetEnergy("13")
CMS.ResetAdditionalInfo()

# A #######################
max_y_coord = fitted_A.GetBinContent(fitted_A.GetMaximumBin())*2
min_y_coord = fitted_A.GetBinContent(fitted_A.GetMinimumBin())*1.5

# Plotting
canv = CMS.cmsCanvas(canv_name,-2.4,2.4,min_y_coord,max_y_coord,"#eta","Magnetic field correction",square=CMS.kSquare,extraSpace=0.07,iPos=iPos)
leg = CMS.cmsLeg(0.68, 0.78, 0.90, 0.90, textSize=0.035)

hdf = CMS.GetcmsCanvasHist(canv)
hdf.GetYaxis().SetMaxDigits(2)
hdf.GetYaxis().SetTitleOffset(1.4)

dummy_A.SetLineColor(ROOT.TColor.GetColor("#e42536"))
dummy_A.SetMarkerColor(ROOT.TColor.GetColor("#e42536"))
fitted_A.SetLineColor(ROOT.TColor.GetColor("#5790fc"))
fitted_A.SetMarkerColor(ROOT.TColor.GetColor("#5790fc"))

# Draw
dummy_A.Draw("same")
fitted_A.Draw("same")

leg.AddEntry(dummy_A, "input A","LP")
leg.AddEntry(fitted_A, "fitted A","LP")

# Output file
canv.SaveAs("InOutputFiles/constants_fitted_plots.pdf")
CMS.SaveCanvas(canv, "InOutputFiles/constants_fitted_plots.root", close=True)

