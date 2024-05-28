#define _USE_MATH_DEFINES
 
#include "TMath.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <limits>

#include <sys/time.h>

using namespace ROOT::Minuit2;

class TheoryFunction {

public:
  
  TheoryFunction(const std::vector<double> &ini_par) : fpar(ini_par) {}
  ~TheoryFunction() {}

  std::vector<double> FPar() const {return fpar;}
  double operator()(double t) const;
  std::vector<double> Grad(double t) const;

private:
  std::vector<double> fpar;
};

// mycomm: function to build functional form of model
double TheoryFunction::operator()(double t) const {
  double asy1(fpar[2]);
  double asy2(fpar[4]);

  //  double expo(TMath::Exp(-TMath::Power(t*fpar[3],fpar[5])));
  double expo(TMath::Exp(-t*fpar[3]));

  double cosine(TMath::Cos(0.0174532925199432955*fpar[0]+6.28318530717958623*fpar[1]*t));

  double asymmetry;
  asymmetry = cosine*expo;
  asymmetry = asy1*asymmetry;
  asymmetry = asymmetry+asy2;

  // mycomm return functional form =  fpar[2]*cos(0.01*fpar[0] + 6.28*fpar[1]*t)*exp(-t*fpar[3]) + fpar[4]
  return asymmetry;
}

//mycomm: function to build gradient of functional form of model
std::vector<double> TheoryFunction::Grad(double t) const {
  double cosarg(0.0174532925199432955*fpar[0]+6.28318530717958623*fpar[1]*t);
  double exparg(t*fpar[3]);
  double cosine(TMath::Cos(cosarg));
  double sine(TMath::Sin(cosarg));
  //  double expo(TMath::Exp(-TMath::Power(exparg,fpar[5])));
  double expo(TMath::Exp(-exparg));

  std::vector<double> grad;

  grad.push_back(-0.0174532925199432955*fpar[2]*sine*expo); // mycomm: derivative wrt fpar[0]
  grad.push_back(-6.28318530717958623*t*fpar[2]*sine*expo); // mycomm: derivative wrt fpar[1]
  grad.push_back(cosine*expo); // mycomm: derivative wrt fpar[2]
  //  if(exparg > 0.0)
  //    grad.push_back(-fpar[2]*cosine*expo*fpar[5]/fpar[3]*TMath::Power(exparg,fpar[5]));
  //  else
  //    grad.push_back(0.0);
  grad.push_back(-t*fpar[2]*cosine*expo); // mycomm: derivative wrt fpar[3]

  grad.push_back(1.0); // mycomm: derivative wrt fpar[4]

  //  if(exparg > 0.0)
  //    grad.push_back(-fpar[2]*cosine*expo*TMath::Log(exparg)*TMath::Power(exparg,fpar[5]));
  //  else
  //    grad.push_back(0.0);
  grad.push_back(0.0); // mycomm: derivative wrt t, which is constant

  // mycomm return array with derivative of functional form wrt to each parameter, t is also a parameter that is not fitted, that's why you push_back 6 times and the last is 0
  return grad;
}

class TheoryFcn : public FCNGradientBase {

public:

  TheoryFcn(const std::vector<double> &meas, const std::vector<double> &err, double t0, double tstep, unsigned int nbins) : tMeasurements(meas),tErrors(err),tTzero(t0),tStep(tstep),tNbins(nbins),tErrorDef(1.) {}

  ~TheoryFcn() {}

  virtual double Up() const {return tErrorDef;}
  virtual double operator()(const std::vector<double>&) const;
  virtual std::vector<double> Gradient(const std::vector<double>& ) const;
  virtual bool CheckGradient() const {return true;}
  //mycomm: I have a feeling these functions have to stay here and with the names like this, for the minimizer to work with them

  std::vector<double> Measurements() const {return tMeasurements;} 
  void SetErrorDef(double def) {tErrorDef = def;} //mycomm: do I need this?

private:

  std::vector<double> tMeasurements;
  std::vector<double> tErrors;
  double tTzero;
  double tStep;
  unsigned int tNbins;
  double tErrorDef;
};

//mycomm: function to build chi2 
double TheoryFcn::operator()(const std::vector<double>& par) const {

  assert(par.size() == 8); //mycomm: why 8?

  double chi2(0.0);
  double diff(0.0);
  double x(0.0);

  TheoryFunction my_func(par);
  for(unsigned int n(0); n < tMeasurements.size() ; n++) {
    x=tTzero+(double)n*tStep;
    if(x>=0.1 && x<=10.0) {
      diff = tMeasurements[n]-(par[6]*TMath::Exp(-x/2.197019)*(1.0 + my_func(x)) + par[7]); // mycomm: par[6] seems to be normalisation signal and [7] for bkg

      chi2 += diff*diff/(tErrors[n]*tErrors[n]);
    }
  }

  return chi2;
}

//mycomm: function to build gradient of chi2 analytically 
std::vector<double> TheoryFcn::Gradient(const std::vector<double> &par ) const {

  assert(par.size() == 8);

  std::vector<double> grad(par.size(),0.0);

  std::vector<double> mygrad, theograd;
  double x(0.0), temp(0.0), myfunc(0.0);

  TheoryFunction my_func(par);
  for(unsigned int n(0); n < tMeasurements.size() ; n++) {
    x=tTzero+(double)n*tStep;

    if(x>=0.1 && x<=10.0) {
      myfunc = my_func(x);
      mygrad = my_func.Grad(x);
      for(unsigned int i(0); i<par.size()-2; i++){
        theograd.push_back(par[6]*TMath::Exp(-x/2.197019)*mygrad[i]); // mycomm: term in derivative of chi2 wrt functional form of model
      }
      theograd.push_back(TMath::Exp(-x/2.197019)*(1.0 + myfunc)); // df/dN // mycomm: term in derivative of chi2 wrt par[6]
      theograd.push_back(1.0); // df/dBg // mycomm: term in derivative of chi2 wrt par[7]

      temp=2.0*((par[6]*TMath::Exp(-x/2.197019)*(1.0 + myfunc) + par[7])-tMeasurements[n])/(tErrors[n]*tErrors[n]);  // mycomm: the minus compensates the lack of minus in the other terms of derivatives

      for(unsigned int i(0); i<par.size(); i++){
        grad[i] += temp*theograd[i]; // mycomm: size of Gradient is = number of measurements, why adding terms?
      }
      mygrad.clear();
      theograd.clear();
    }
  }

  return grad;
}

int main() {
  // dummy data
  double tstart=0.0, tstep=0.01;
  unsigned int noOfBins = 20;
  
  std::vector<double> data, error;

  for (unsigned int i=0; i<noOfBins; i++){
    data.push_back(i + 0.1);
    error.push_back(0.05);
  }


  TheoryFcn fFCN(data,error,tstart,tstep,noOfBins);

  // create parameters with initial starting values

  //mycomm: declare parameters here

  MnUserParameters upar;
  upar.Add("phase", 16.7545, 1.71896, 0.0, 100.0);
  upar.Add("freq", 0.37145, 0.02337, 0.0, 300.0);
  upar.Add("asymS", 0.2, 0.0632559, 0.0, 0.3);
  upar.Add("rateS", 0.6563, 0.07805, 0.0, 100.0);
  upar.Add("asymF", 0.0);
  upar.Add("betaS", 1.0);
  upar.Add("Norm_L", 351.856, 0.770441);
  upar.Add("BG_L", 35.4337, 0.104235);

  for (unsigned int i=0; i<upar.Params().size(); i++) {
    std::cout <<"par[" << i << "]: " << upar.Params()[i] << std::endl;
  }

  // create Migrad minimizer

  char a = 'a';
  char* a_ptr = &a;

  MnMigrad minimize(fFCN, upar, 1);

  // ... and Minimize

  unsigned int maxfcn(std::numeric_limits<unsigned int>::max());
  double tolerance(0.1);

  double t1(0.);
  // get start time
  struct timeval tv_start, tv_stop;
  gettimeofday(&tv_start, 0);

  FunctionMinimum min = minimize(maxfcn, tolerance);

  gettimeofday(&tv_stop, 0);
  t1 = (tv_stop.tv_sec - tv_start.tv_sec)*1000.0 + (tv_stop.tv_usec - tv_start.tv_usec)/1000.0;

  std::cout << "# Calculation time: " << t1 << " ms" << std::endl;

  std::cout << "CHI^2: " << min.Fval() << std::endl << std::endl;

  std::cout << "min is valid: " << min.IsValid() << std::endl;
  std::cout << "HesseFailed: " << min.HesseFailed() << std::endl;
  std::cout << "HasCovariance: " << min.HasCovariance() << std::endl;
  std::cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
  std::cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
  std::cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
  std::cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
  std::cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
  std::cout << "HasPosDefCovar : " << min.HasPosDefCovar() << std::endl;
  std::cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar() << std::endl;

  std::cout << min << "\n";

  for (unsigned int i(0); i<upar.Params().size(); i++) {
    std::cout <<"par[" << i << "]: " << min.UserState().Value(i) << std::endl;
  }

  return 0;
}


