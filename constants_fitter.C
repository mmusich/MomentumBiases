#define _USE_MATH_DEFINES
 
#include "TMath.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"

//#include <vector>
//#include <iostream>
//#include <cassert>
//#include <limits>

#include <sys/time.h>

using namespace ROOT::Minuit2;
//TODO use namespace std;

class TheoryFcn : public FCNGradientBase {

public:

  TheoryFcn(const std::vector<double> &meas, const std::vector<double> &err, const std::vector<std::string> &bin_labels, unsigned int nbins) : scaleSquared(meas),scaleSquaredError(err),binLabels(bin_labels),nBins(nbins), errorDef(1.0) {}
  //TODO do I need Nbins for minuit to accept this?
  ~TheoryFcn() {}

  virtual double Up() const {return errorDef;}
  virtual double operator()(const std::vector<double>&) const;
  virtual std::vector<double> Gradient(const std::vector<double>& ) const;
  virtual bool CheckGradient() const {return true;}
  //mycomm: figure out from documentatiosn which functions i want to overwrite 

  std::vector<int> getIndices(std::string bin_label) const; // Function to get bin index from label name
  double getK(const int pT_index) const;

  const std::vector<double> pT_binning {38.0, 39.3584, 40.4562, 41.2942, 42.9469, 43.0}; //{25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0}; //TODO am I allowed to use const here? //TODO change pT after debugging
  static constexpr double scaling_A = 0.001, scaling_e = 0.001 * 40.0, scaling_M = 0.001 / 40.0;
  static const int n_eta_bins = 2; //TODO automatize this from size of bin labels

  //TODO can I remove these?
  //std::vector<double> Measurements() const {return tMeasurements;} 
  //void SetErrorDef(double def) {tErrorDef = def;} //mycomm: do I need this?

private:

  std::vector<double> scaleSquared;
  std::vector<double> scaleSquaredError;
  std::vector<std::string> binLabels;
  unsigned int nBins; //TODO write method to get nBins, minuit shouldn't mind
  double errorDef;
};

//-----------------------------------------------
// Function to get bin index from label name

std::vector<int> TheoryFcn::getIndices(std::string bin_label) const{
  std::vector<int> indices;
  int pos = 0;
  while(pos < bin_label.size()){
    pos = bin_label.find("_");
    indices.push_back(stoi(bin_label.substr(0,pos)));
    bin_label.erase(0,pos+1); // 1 is the length of the delimiter, "_"
  }
  return indices;
}

//-----------------------------------------------
// Function to get average k from pT bin index

double TheoryFcn::getK(const int pT_index) const{
  // TODO what about errors on k?
  // TODO could add the pT histo per eta,pt,eta,pt to get more accurate average pT
  return 2.0 / (pT_binning[pT_index] + pT_binning[pT_index+1]); // k = 1 / avg_pT
}

//-----------------------------------------------
//mycomm: function to build chi2 

double TheoryFcn::operator()(const std::vector<double>& par) const {
  // par has size n_parameters = 3*number_eta_bins, idices from 0 to n-1 contain A, from n to 2n-1 epsilon, from 2n to 3n-1 M

  assert(par.size() == 3*n_eta_bins); 

  double chi2(0.0);
  double diff(0.0);

  double k_plus, k_minus, term_pos, term_neg;
  double my_func;
  int eta_pos_index, eta_neg_index;
  std::vector<int> bin_indices(4); // for 4D binning

  for(unsigned int n(0); n < scaleSquared.size() ; n++) {
    bin_indices = getIndices(binLabels[n]); //TODO I hope overwriting this vector is ok
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]); 
    k_minus = getK(bin_indices[3]);

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of e,M
    term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e*par[eta_pos_index + n_eta_bins]*k_plus +  scaling_M*par[eta_pos_index + 2*n_eta_bins]/k_plus);
    term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e*par[eta_neg_index + n_eta_bins]*k_minus - scaling_M*par[eta_neg_index + 2*n_eta_bins]/k_minus);
    
    my_func = term_pos*term_neg;
    
    diff = my_func - scaleSquared[n]; 
    chi2 += diff*diff/(scaleSquaredError[n]*scaleSquaredError[n]);
  }

  return chi2;
}

//-----------------------------------------------
//mycomm: function to build gradient of chi2 analytically 
std::vector<double> TheoryFcn::Gradient(const std::vector<double> &par ) const {

  assert(par.size() == 3*n_eta_bins);

  std::vector<double> grad(par.size(),0.0);

  double temp(0.0), local_func(0.0), local_grad(0.0);

  double k_plus, k_minus, term_pos, term_neg;
  int eta_pos_index, eta_neg_index;
  std::vector<int> bin_indices(4); // for 4D binning
  
  for(unsigned int n(0); n < scaleSquared.size() ; n++) { // loops over measurements
    // TODO I need everything from the chi2 function, maybe better to have operator() to return a tuple with the function and the grad vector
    bin_indices = getIndices(binLabels[n]); //TODO I hope overwriting this vector is ok
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]); 
    k_minus = getK(bin_indices[3]);

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of e,M
    term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e*par[eta_pos_index + n_eta_bins]*k_plus +  scaling_M*par[eta_pos_index + 2*n_eta_bins]/k_plus);
    term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e*par[eta_neg_index + n_eta_bins]*k_minus - scaling_M*par[eta_neg_index + 2*n_eta_bins]/k_minus);
    
    local_func = term_pos*term_neg;
    
    temp=2.0*(local_func - scaleSquared[n])/(scaleSquaredError[n]*scaleSquaredError[n]); 
    
    for (unsigned int i : {eta_pos_index, eta_pos_index + n_eta_bins, eta_pos_index + 2*n_eta_bins, eta_neg_index, eta_neg_index + n_eta_bins, eta_neg_index + 2*n_eta_bins}) { // loops over parameters for A,e,M model
      if (i == eta_pos_index) { // mycomm: derivative wrt A(+), which is fpar[eta_pos_index] for this bin label
	local_grad = scaling_A*term_neg;
      } else if (i == eta_pos_index + n_eta_bins) { // derivative wrt e(+)
      	local_grad = -scaling_e*k_plus*term_neg;
      } else if (i == eta_pos_index + 2*n_eta_bins) { // derivative wrt M(+)
	local_grad = scaling_M/k_plus*term_neg;
      } else if (i == eta_neg_index){ // derivative wrt A(-)
	local_grad = scaling_A*term_pos;
      } else if (i == eta_neg_index + n_eta_bins) { // derivative wrt e(-)
	local_grad = -scaling_e*k_minus*term_pos;
      } else if (i == eta_neg_index + 2*n_eta_bins){ // derivative wrt M(-)
	local_grad = -scaling_M/k_minus*term_pos;
      } else {
	std::cout<<"\n"<<"ERROR: indices in grad don't match"<<"\n"; //TOOD raise a proper error
      }
      
      grad[i] += temp*local_grad; // grad has size 3*n_eta_bins
    }
  }
  // mycomm: returns vector size nparam with derivatives of chi2 wrt to each parameter
  return grad;
}

//-----------------------------------------------
// Function to get parameter name (A0, e1, etc.) from index in parameters array

std::tuple<string,double> getParameterNameAndScaling(int index){

  int whole = index / TheoryFcn::n_eta_bins;
  int rest = index % TheoryFcn::n_eta_bins;
  double scaling;

  string name;
  if (whole == 0) {
    name = "A" + to_string(rest);
    scaling = TheoryFcn::scaling_A;
  }
  else if (whole == 1) {
    name = "e" + to_string(rest);
    scaling = TheoryFcn::scaling_e;
  }
  else if (whole == 2) {
    name = "M" + to_string(rest);
    scaling = TheoryFcn::scaling_M;
  }
  else {
    std::cout<<"\n"<<"ERROR counting parameters"<<"\n";
  }

  return make_tuple(name, scaling);
}

int main() {

  // dummy data
  unsigned int noOfBins = 12;
  unsigned int n_parameters = 6; //TODO refine after debugging 

  std::vector<double> data{0.993*0.993, 0.996*0.996, 0.992*0.992, 0.996*0.996, 0.998*0.998, 0.990*0.990, 0.992*0.992, 0.995*0.995, 0.994*0.994, 0.998*0.998, 0.996*0.996, 0.990*0.990};
  std::vector<double> error{2*0.993*0.0005, 2*0.996*0.0005, 2*0.992*0.0005, 2*0.996*0.0005, 2*0.998*0.0005, 2*0.990*0.0005, 2*0.992*0.0005, 2*0.995*0.0005, 2*0.994*0.0005, 2*0.998*0.0005, 2*0.996*0.0005, 2*0.990*0.0005};
  std::vector<std::string> labels{"0_0_1_1", "0_0_1_2", "0_0_1_3", "0_0_1_4", "0_1_1_0", "0_1_1_2", "0_1_1_3", "0_1_1_4", "0_2_1_0", "0_2_1_1", "0_2_1_3", "0_2_1_4"};

  TheoryFcn fFCN(data,error,labels,noOfBins);

  // create parameters with initial starting values
  std::vector<double> step(n_parameters, 0.1);
  std::vector<double> start(n_parameters, 0.0);

  MnUserParameters upar;
  for(unsigned int i=0 ; i<n_parameters; ++i){
    upar.Add(Form("param%d",i), start[i]); // TODO: bool ROOT::Minuit2::MnUserParameters::Add(const std::string & name, double val,double err,double low,double up )
    //TODO how do I set step? 
  }

  for (unsigned int i=0; i<upar.Params().size(); i++) {
    std::cout <<"par[" << i << "]: " << upar.Params()[i] << std::endl;
  }

  // create Migrad minimizer

  MnMigrad minimize(fFCN, upar, 1); //TODO strategy is 1, check others too

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

  std::cout << std::endl << "Fitted A,e,M parameters: " << std::endl;
  for (unsigned int i(0); i<upar.Params().size(); i++) {
    std::cout <<"par[" << i << "]: " << get<0>(getParameterNameAndScaling(i)) << " fitted to: "<< min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i)) << std::endl;
  }

  return 0;
}


