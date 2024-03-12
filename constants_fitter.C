// Stand alone code to fit for alignment parameters A,e,M from mass fits 

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"

#include <sys/time.h>

using namespace ROOT::Minuit2;
using namespace std;

class TheoryFcn : public FCNGradientBase {

public:

  TheoryFcn(const vector<double> &meas, const vector<double> &err, const vector<string> &bin_labels) : scaleSquared(meas), scaleSquaredError(err), binLabels(bin_labels), errorDef(1.0) {}
  ~TheoryFcn() {}

  virtual double Up() const {return errorDef;}
  virtual void SetErrorDef(double def) {errorDef = def;}

  virtual double operator()(const vector<double>&) const;
  virtual vector<double> Gradient(const vector<double>& ) const;
  virtual bool CheckGradient() const {return true;}
  // virtual std::vector< double > ROOT::Minuit2::FCNGradientBase::Hessian(const std::vector< double > & )const -> allows to do Hessian analytically? 

  vector<int> getIndices(string bin_label) const; 
  double getK(const int pT_index) const;

  const vector<double> pT_binning {25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0}; // pT binning goes here
  static constexpr double scaling_A = 0.01, scaling_e = 0.01 * 40.0, scaling_M = 0.01 / 40.0; // this scaling makes fitter parameters of order 1
  static const int n_eta_bins = 2; //TODO automatize this from size of bin labels and use the code that checks if we have enough data points for set of parameter

private:

  vector<double> scaleSquared;
  vector<double> scaleSquaredError;
  vector<string> binLabels;
  double errorDef;
};

//-----------------------------------------------
// Function to get bin index from label name

vector<int> TheoryFcn::getIndices(string bin_label) const{
  vector<int> indices;
  int pos = 0;
  for (int i=0; i<4-1;i++){ //4 for 4D binning
    pos = bin_label.find("_");
    indices.push_back(stoi(bin_label.substr(0,pos)));
    bin_label.erase(0,pos+1); // 1 is the length of the delimiter, "_"
  }
  indices.push_back(stoi(bin_label));
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
// Function to build chi2 

double TheoryFcn::operator()(const vector<double>& par) const {
  // par has size n_parameters = 3*number_eta_bins, idices from 0 to n-1 contain A, from n to 2n-1 epsilon, from 2n to 3n-1 M
  assert(par.size() == 3*n_eta_bins); 

  double chi2(0.0);
  double diff(0.0);

  double k_plus, k_minus, term_pos, term_neg;
  // TODO for debugging
  double k_plus_values[12] = {1.0/38.5, 1.0/45.2, 1.0/40.5, 1.0/44.2, 1.0/39.8, 1.0/37.4, 1.0/33.9, 1.0/46.2, 1.0/34.9, 1.0/45.9, 1.0/34.2, 1.0/38.4};
  double k_minus_values[12] = {1.0/35.9, 1.0/39.6, 1.0/42.0, 1.0/44.8, 1.0/41.0, 1.0/45.2, 1.0/39.1, 1.0/38.2, 1.0/35.0, 1.0/40.2, 1.0/46.9, 1.0/34.8};
  
  double my_func;
  int eta_pos_index, eta_neg_index;
  vector<int> bin_indices(4); // for 4D binning

  for(unsigned int n(0); n < scaleSquared.size() ; n++) {
    bin_indices = getIndices(binLabels[n]); //TODO I hope overwriting this vector is ok
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]); 
    k_minus = getK(bin_indices[3]);
    //TODO change after debugging
    //k_plus = k_plus_values[n];
    //k_minus = k_minus_values[n];

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
// Function to build gradient of chi2 analytically 

vector<double> TheoryFcn::Gradient(const vector<double> &par ) const {

  assert(par.size() == 3*n_eta_bins);

  vector<double> grad(par.size(),0.0);

  double temp(0.0), local_func(0.0), local_grad(0.0);

  double k_plus, k_minus, term_pos, term_neg;
  // TODO for debugging
  double k_plus_values[12] = {1.0/38.5, 1.0/45.2, 1.0/40.5, 1.0/44.2, 1.0/39.8, 1.0/37.4, 1.0/33.9, 1.0/46.2, 1.0/34.9, 1.0/45.9, 1.0/34.2, 1.0/38.4};
  double k_minus_values[12] = {1.0/35.9, 1.0/39.6, 1.0/42.0, 1.0/44.8, 1.0/41.0, 1.0/45.2, 1.0/39.1, 1.0/38.2, 1.0/35.0, 1.0/40.2, 1.0/46.9, 1.0/34.8};

  int eta_pos_index, eta_neg_index;
  vector<int> bin_indices(4); // for 4D binning
  
  for(unsigned int n(0); n < scaleSquared.size() ; n++) { // loops over measurements
    // TODO I need everything from the chi2 function, maybe better to have operator() to return a tuple with the function and the grad vector
    bin_indices = getIndices(binLabels[n]); //TODO I hope overwriting this vector is ok
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]); 
    k_minus = getK(bin_indices[3]);
    //TODO change after debugging
    //k_plus = k_plus_values[n];
    //k_minus = k_minus_values[n];

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of e,M
    term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e*par[eta_pos_index + n_eta_bins]*k_plus +  scaling_M*par[eta_pos_index + 2*n_eta_bins]/k_plus);
    term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e*par[eta_neg_index + n_eta_bins]*k_minus - scaling_M*par[eta_neg_index + 2*n_eta_bins]/k_minus);
    
    local_func = term_pos*term_neg;
    
    temp=2.0*(local_func - scaleSquared[n])/(scaleSquaredError[n]*scaleSquaredError[n]); 
    
    //TODO not calculate terms before checking if needed
    for (unsigned int i : {eta_pos_index, eta_pos_index + n_eta_bins, eta_pos_index + 2*n_eta_bins, eta_neg_index, eta_neg_index + n_eta_bins, eta_neg_index + 2*n_eta_bins}) { // loops over parameters
      if (i == eta_pos_index) { // derivative wrt A(+), which is fpar[eta_pos_index] 
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
	cout<<"\n"<<"ERROR: indices in grad don't match"<<"\n"; //TOOD raise a proper error
      }
      
      grad[i] += temp*local_grad; 
    }
  }
  return grad;
}

//-----------------------------------------------
// Function to get parameter name (A0, e1, etc.) and appropriate scaling from index in parameters array

tuple<string,double> getParameterNameAndScaling(int index){

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
    cout<<"\n"<<"ERROR counting parameters"<<"\n";
  }

  return make_tuple(name, scaling);
}

int main() {

  // dummy data
  unsigned int n_parameters = 6; //TODO refine after debugging 

  vector<double> data{0.990479, 0.989859, 0.98801, 0.984616, 0.992561, 0.986182, 0.986833, 0.986570, 0.993423, 0.988375, 0.987912, 0.988051}; // smeared values
  //vector<double> data{0.989709, 0.988419, 0.98794, 0.988019, 0.991833, 0.986387, 0.985908, 0.985987, 0.992988, 0.988824, 0.987057, 0.987135};
  /* params used to generate data for closure test
              A0 fitted to: -0.095517
      par[1]: A1 fitted to: -0.0950043
      par[2]: e0 fitted to: -2.20665
      par[3]: e1 fitted to: -1.21125
      par[4]: M0 fitted to: 0.00179058
      par[5]: M1 fitted to: -0.000555938
  */
  //vector<double> error{0.989709*0.001,0.988419*0.001, 0.98794*0.001, 0.988019*0.001, 0.991833*0.001, 0.986387*0.001, 0.985908*0.001, 0.985987*0.001, 0.992988*0.001, 0.988824*0.001, 0.987057*0.001, 0.987135*0.001 };
  vector<double> error{0.001, 0.001,0.001, 0.001,0.001, 0.001,0.001, 0.001,0.001, 0.001,0.001, 0.001};
  vector<string> labels{"0_0_1_1", "0_0_1_2", "0_0_1_3", "0_0_1_4", "0_1_1_0", "0_1_1_2", "0_1_1_3", "0_1_1_4", "0_2_1_0", "0_2_1_1", "0_2_1_3", "0_2_1_4"};

  TheoryFcn fFCN(data,error,labels);

  // create parameters with initial starting values
  vector<double> par_error(n_parameters, 1.0); //TODO is error on parameter before taking scaling into account
  vector<double> start(n_parameters, 0.0); // TODO if these are consts for every parameter, no need to make this a vector, but if we want to iterate fit is good

  MnUserParameters upar;
  for(unsigned int i=0 ; i<n_parameters; ++i){
    upar.Add(Form("param%d",i), start[i], par_error[i]);
    //TODO how do I set step? 
  }

  for (unsigned int i=0; i<upar.Params().size(); i++) {
    cout <<"par[" << i << "]: " << upar.Params()[i] << "\n";
  }

  // create Migrad minimizer

  MnMigrad minimize(fFCN, upar, 1); //TODO strategy is 1, check others too

  // ... and Minimize

  unsigned int maxfcn(1000000); //(numeric_limits<unsigned int>::max());
  double tolerance(0.001); //MIGRAD will stop iterating when edm (expected distance from minimum) will be: edm < tolerance * 10**-3

  double t1(0.);
  // get start time
  struct timeval tv_start, tv_stop;
  gettimeofday(&tv_start, 0);

  FunctionMinimum min = minimize(maxfcn, tolerance);  

  gettimeofday(&tv_stop, 0);
  t1 = (tv_stop.tv_sec - tv_start.tv_sec)*1000.0 + (tv_stop.tv_usec - tv_start.tv_usec)/1000.0;

  cout << "# Calculation time: " << t1 << " ms" << "\n";

  cout << "CHI^2: " << min.Fval() << "\n" << "\n";

  cout << "min is valid: " << min.IsValid() << std::endl;
  cout << "HesseFailed: " << min.HesseFailed() << std::endl;
  cout << "HasCovariance: " << min.HasCovariance() << std::endl;
  cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
  cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
  cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
  cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
  cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
  cout << "HasPosDefCovar : " << min.HasPosDefCovar() << std::endl;
  cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar() << std::endl;

  cout << min << "\n";

  cout << "\n" << "Fitted A [ ], e [GeV], M [GeV^-1] parameters: " << "\n";
  for (unsigned int i(0); i<upar.Params().size(); i++) {
    cout <<"par[" << i << "]: " << get<0>(getParameterNameAndScaling(i)) << " fitted to: "<< min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i)) << "\n";
  }

  return 0;
}


