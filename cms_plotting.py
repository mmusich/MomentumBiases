# ROOT imports 
import os, ROOT
import cmsstyle as CMS

outfile = ROOT.TFile.Open("InOutputFiles/constants_fitted_plots.root","RECREATE")

# File reading 
f_constants_fitted = ROOT.TFile.Open('InOutputFiles/constants_fitted_for_plotting.root')

dict = {
    'A': 'Magnetic field correction',
    'e': 'Material correction',
    'M': 'Misalignment correction'
}

for key in dict:

    par_name = key
    par_title = dict[key]
    dummy_par = f_constants_fitted.Get(f"dummy_pars_{par_name}")
    fitted_par = f_constants_fitted.Get(f"fitted_pars_{par_name}")

    # Styling
    CMS.SetExtraText("Simulation Private Work")
    canv_name = par_name
    iPos = 11 
    CMS.SetLumi("")
    CMS.SetEnergy("13")
    CMS.ResetAdditionalInfo()

    max_y_coord = fitted_par.GetBinContent(fitted_par.GetMaximumBin())*2
    min_y_coord = fitted_par.GetBinContent(fitted_par.GetMinimumBin())*2
    
    # Plotting
    canv = CMS.cmsCanvas(canv_name,-2.4,2.4,min_y_coord,max_y_coord,"#eta",par_title,square=CMS.kSquare,extraSpace=0.07,iPos=iPos)
    leg = CMS.cmsLeg(0.68, 0.78, 0.90, 0.90, textSize=0.035)
    
    hdf = CMS.GetcmsCanvasHist(canv)
    hdf.GetYaxis().SetMaxDigits(2)
    hdf.GetYaxis().SetTitleOffset(1.4)
    
    dummy_par.SetLineColor(ROOT.TColor.GetColor("#e42536"))
    dummy_par.SetMarkerColor(ROOT.TColor.GetColor("#e42536"))
    dummy_par.SetMarkerStyle(6)
    fitted_par.SetLineColor(ROOT.TColor.GetColor("#5790fc"))
    fitted_par.SetMarkerColor(ROOT.TColor.GetColor("#5790fc"))
    fitted_par.SetMarkerStyle(6)
    
    # Draw
    dummy_par.Draw("same")
    fitted_par.Draw("same")
    
    if key == 'e':
        title = "#e"
        leg.AddEntry(dummy_par, title,"LP")
        leg.AddEntry(fitted_par, title,"LP")
    else:
        leg.AddEntry(dummy_par, f"input {par_name}","LP")
        leg.AddEntry(fitted_par, f"fitted {par_name}","LP")
        
    # Output file
    canv.SaveAs(f"InOutputFiles/PresentablePlots/constants_fitted_plots_{par_name}.pdf")
    #outfile.cd()
    #canv.Write()
    #outfile.Close()
    #CMS.SaveCanvas(canv, "InOutputFiles/constants_fitted_plots.root", close=False)

