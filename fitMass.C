#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <string>

// RooFit includes
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "TStopwatch.h" // for measuring time

bool debugMode_ = true;

/** A class that combines a value and it's associated uncertainty,
 *  or error, together. Provides a more explicit interface than
 *  a pair<double,double>. If you don't like the name, propose a better one!
 */

class Measurement1D {

public:
  // construct
  Measurement1D() : theValue(0.), theError(0.){};
  Measurement1D(const double& aValue) : theValue(aValue), theError(0.){};
  Measurement1D(const double& aValue, const double& aError) : theValue(aValue), theError(aError){};


  //destruct
  ~Measurement1D(){};

  double value() const { return theValue; }
  double error() const { return theError; }
  double significance() const {
    if (theError == 0)     
      return 0;
    else
      return theValue / theError;    
  }

private:
  double theValue;
  double theError;
};

// to store the data
namespace diMuonMassBias {

  struct fitOutputs {
  public:
    fitOutputs(const Measurement1D& bias, const Measurement1D& width) : m_bias(bias), m_width(width) {}

    // getters
    const Measurement1D getBias() { return m_bias; }
    const Measurement1D getWidth() { return m_width; }
    bool isInvalid() {
      return (m_bias.value() == 0.f && m_bias.error() == 0.f && m_width.value() == 0.f && m_width.error() == 0.f);
    }

  private:
    Measurement1D m_bias;
    Measurement1D m_width;
  };

  static constexpr int minimumHits = 10;
}

//-----------------------------------------------------------------------------------
diMuonMassBias::fitOutputs fitLineShape(TH1* hist, const bool& fitBackground_, const bool& useRooCBShape_)
//-----------------------------------------------------------------------------------
{
  if (hist->GetEntries() < diMuonMassBias::minimumHits) {
    std::cout << " Input histogram:" << hist->GetName() << " has not enough entries ("
	      << hist->GetEntries() << ") for a meaningful Voigtian fit!\n"
	      << "Skipping!";

    return diMuonMassBias::fitOutputs(Measurement1D(0., 0.), Measurement1D(0., 0.));
  }

  TCanvas* c1 = new TCanvas();
  if (debugMode_) {
    c1->Clear();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.10);
  }

  // silence messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  Double_t xmin = hist->GetXaxis()->GetXmin();
  Double_t xmax = hist->GetXaxis()->GetXmax();

  if (debugMode_) {
    std::cout << "fitting range: (" << xmin << "-" << xmax << ")" << std::endl;
  }

  RooRealVar InvMass("InvMass", "di-muon mass M(#mu^{+}#mu^{-}) [GeV]", xmin, xmax);
  std::unique_ptr<RooPlot> frame{InvMass.frame()};
  RooDataHist datahist("datahist", "datahist", InvMass, RooFit::Import(*hist));
  datahist.plotOn(frame.get());

  // parameters of the Voigtian
  RooRealVar mean("#mu", "mean",90.0, 60.0, 120.0); // (for Z)
  RooRealVar width("width", "width", 5.0,  0.0, 120.0); // (for Z)
  RooRealVar sigma("#sigma", "sigma", 5.0,  0.0, 120.0); // (for Z)
  RooVoigtian voigt("voigt", "voigt", InvMass, mean, width, sigma);

  // parameters of the Crystal-ball
  RooRealVar peakCB("peakCB", "peakCB", 90.0, 60.0, 120.0); // (for Z)
  RooRealVar sigmaCB("#sigma", "sigma", 5.0,  0.0, 120.0); // (for Z)
  RooRealVar alphaCB("#alpha", "alpha", 1., 0., 10.);
  RooRealVar nCB("n", "n", 1., 0., 100.);
  RooCBShape crystalball("crystalball", "crystalball", InvMass, peakCB, sigmaCB, alphaCB, nCB);

  // for the simple background fit
  RooRealVar lambda("#lambda", "slope", 0., -50., 50.);
  RooExponential expo("expo", "expo", InvMass, lambda);

  // define the signal and background fractions
  RooRealVar b("N_{b}", "Number of background events", 0, hist->GetEntries() / 10.);
  RooRealVar s("N_{s}", "Number of signal events", 0, hist->GetEntries());

  if (fitBackground_) {
    RooArgList listPdf;
    if (useRooCBShape_) {     
      // crystal-ball + exponential fit
      listPdf.add(crystalball);
      listPdf.add(expo);     
    } else {
      // voigtian + exponential fit
      listPdf.add(voigt);
      listPdf.add(expo);
    }
  
    RooAddPdf fullModel("fullModel", "Signal + Background Model", listPdf, RooArgList(s, b));
    fullModel.fitTo(datahist, RooFit::PrintLevel(-1));
    fullModel.plotOn(frame.get(), RooFit::LineColor(kRed));
    fullModel.plotOn(frame.get(), RooFit::Components(expo), RooFit::LineStyle(kDashed));  //Other option
    fullModel.paramOn(frame.get(), RooFit::Layout(0.65, 0.90, 0.90));
  } else {
    if (useRooCBShape_) {
      // use crystal-ball for a fit-only signal
      crystalball.fitTo(datahist, RooFit::PrintLevel(-1));
      crystalball.plotOn(frame.get(), RooFit::LineColor(kRed));  //this will show fit overlay on canvas
      crystalball.paramOn(frame.get(),
                          RooFit::Layout(0.65, 0.90, 0.90));  //this will display the fit parameters on canvas
    } else {
      // use voigtian for a fit-only signal
      voigt.fitTo(datahist, RooFit::PrintLevel(-1));
      voigt.plotOn(frame.get(), RooFit::LineColor(kRed));            //this will show fit overlay on canvas
      voigt.paramOn(frame.get(), RooFit::Layout(0.65, 0.90, 0.90));  //this will display the fit parameters on canvas
    }
  }

  // Redraw data on top and print / store everything
  datahist.plotOn(frame.get());
  frame->GetYaxis()->SetTitle("n. of events");
  TString histName = hist->GetName();
  frame->SetName("frame" + histName);
  frame->SetTitle(hist->GetTitle());
  frame->Draw();

  if (debugMode_) {
    c1->Print("fit_debug" + histName + ".pdf");
  }
  delete c1;

  float mass_mean = useRooCBShape_ ? peakCB.getVal() : mean.getVal();
  float mass_sigma = useRooCBShape_ ? sigmaCB.getVal() : sigma.getVal();

  float mass_mean_err = useRooCBShape_ ? peakCB.getError() : mean.getError();
  float mass_sigma_err = useRooCBShape_ ? sigmaCB.getError() : sigma.getError();

  Measurement1D resultM(mass_mean, mass_mean_err);
  Measurement1D resultW(mass_sigma, mass_sigma_err);

  return diMuonMassBias::fitOutputs(resultM, resultW);
}

// main function
int fitMass(TString input="data_histos_cuts.root"){

  // start a CPU timer
  TStopwatch timer;
  timer.Start();

  TFile *f = TFile::Open(input);
  TH1F* histo = nullptr;
  for(unsigned int eta1bin=1; eta1bin<=6; eta1bin++){
    for(unsigned int phi1bin=1; phi1bin<=6; phi1bin++){
      for(unsigned int eta2bin=1; eta2bin<=6; eta2bin++){
	for(unsigned int phi2bin=1; phi2bin<=6; phi2bin++){
	  TString address = Form("data_frame_%i_%i_%i_%i",eta1bin,phi1bin,eta2bin,phi2bin); 
	  histo = (TH1F*)f->Get(address);
	  if(histo){
	    std::cout << "histo: " << address << " has " << histo->GetEntries() << " entries!" << std::endl;
	    const auto result = fitLineShape(histo,true /* fit background */, false /* use Crystall-ball */);
	  } else {
	    return EXIT_FAILURE;
	  } // if histo is valid
	} // loop on second phi
      } // loop on second eta
    } // loop on first phi
  } // loop on first eta

  timer.Stop();
  
  // print some timing statistics
  double rtime = timer.RealTime();
  double ctime = timer.CpuTime();
  // timing printouts
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
    
  return 0;
} // end main


