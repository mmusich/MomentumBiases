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
  static constexpr double scaling_A = 0.001, scaling_e = 0.001 * 40.0, scaling_M = 0.001 / 40.0; // this scaling makes fitter parameters of order 1 TODO diff scaling for diff pt bin?
  static const int n_eta_bins = 2; //TODO automatize this from size of bin labels and use the code that checks if we have enough data points for set of parameter
  
  vector<string> binLabels;

private:

  vector<double> scaleSquared;
  vector<double> scaleSquaredError;
  //vector<string> binLabels; //TODO move to private after debugging
  double errorDef;
};

//-----------------------------------------------
class TheoryFcn2 : public TheoryFcn {

public:
  
  TheoryFcn2(const vector<double> &meas, const vector<double> &err, const vector<string> &bin_labels) : TheoryFcn(meas, err, bin_labels) {}
  ~TheoryFcn2() {}

  vector<vector<double>> DummyData(const int n_data_points, const double width) const; //TODO  vector<vector<double>> might be slow
  
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
  //double k_plus_values[12] = {1.0/38.5, 1.0/45.2, 1.0/40.5, 1.0/44.2, 1.0/39.8, 1.0/37.4, 1.0/33.9, 1.0/46.2, 1.0/34.9, 1.0/45.9, 1.0/34.2, 1.0/38.4};
  //double k_minus_values[12] = {1.0/35.9, 1.0/39.6, 1.0/42.0, 1.0/44.8, 1.0/41.0, 1.0/45.2, 1.0/39.1, 1.0/38.2, 1.0/35.0, 1.0/40.2, 1.0/46.9, 1.0/34.8};
  
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

  int eta_pos_index, eta_neg_index;
  vector<int> bin_indices(4); // for 4D binning
  
  for(unsigned int n(0); n < scaleSquared.size() ; n++) { // loops over measurements
    // TODO I need everything from the chi2 function, maybe better solution to this
    bin_indices = getIndices(binLabels[n]); // TODO is overwriting a vector like this safe?
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
    
    //TODO not calculate terms before checking if needed, ahm, not needed
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
// Function to generate dummy data for closure test

vector<vector<double>> TheoryFcn2::DummyData(const int n_data_points, const double width) const {
  // par has size n_parameters = 3*number_eta_bins, idices from 0 to n-1 contain A, from n to 2n-1 epsilon, from 2n to 3n-1 M
  
  vector<vector<double>> res(2);
  vector<double> test_data(n_data_points, 0.0), test_error(n_data_points, 0.0);
  double k_plus, k_minus, mean, width_rand, term_pos, term_neg;
  TRandom3 random = TRandom3(4357); 
  int eta_pos_index, eta_neg_index;

  vector<int> bin_indices(4); // for 4D binning
  vector<double> dummy_par_val{0.0003, 0.0002, 0.001, 0.002, -0.00001, -0.00002}; //TODO automatise

  /* // Parameter values used for dummy data
      par[0]: A0 0.0003
      par[1]: A1 0.0002
      par[2]: e0 0.001
      par[3]: e1 0.002
      par[4]: M0 -0.00001
      par[5]: M1 -0.00002
  */

  for(unsigned int n(0); n <n_data_points ; n++) {
    bin_indices = getIndices(binLabels[n]); //TODO I hope overwriting this vector is ok
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]);
    k_minus = getK(bin_indices[3]);
  
    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of e,M
    term_pos = (1. + dummy_par_val[eta_pos_index] - dummy_par_val[eta_pos_index + n_eta_bins]*k_plus +  dummy_par_val[eta_pos_index + 2*n_eta_bins]/k_plus);
    term_neg = (1. + dummy_par_val[eta_neg_index] - dummy_par_val[eta_neg_index + n_eta_bins]*k_minus - dummy_par_val[eta_neg_index + 2*n_eta_bins]/k_minus);

    mean = term_pos*term_neg;
    width_rand = random.Gaus(width, width/10.0);
    test_error[n] = width_rand;
    test_data[n] = random.Gaus(mean, width_rand);
    
  }

  res[0] = test_data;
  res[1] = test_error;

  return res;
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

//-----------------------------------------------
// Overload << operator

template <typename S>
ostream& operator<<(ostream& os,
                    const vector<S>& vector)
{
  for (auto element : vector) {
    os << element << " ";
  }
  return os;
}

//-----------------------------------------------

int main() {

  // dummy data
  unsigned int n_parameters = 6; //TODO refine after debugging 
  unsigned int n_data_points = 50; //TODO refine after debugging
  
  vector<double> empty_data(n_data_points, 0.0), empty_error(n_data_points, 0.0); // for closure test
  double error_start = 0.001;
  vector<string> labels{"0_0_1_0","0_0_1_1","0_0_1_2","0_0_1_3","0_0_1_4","0_1_1_0","0_1_1_1","0_1_1_2","0_1_1_3","0_1_1_4","0_2_1_0","0_2_1_1","0_2_1_2","0_2_1_3","0_2_1_4", "0_3_1_0","0_3_1_1","0_3_1_2","0_3_1_3","0_3_1_4", "0_4_1_0","0_4_1_1","0_4_1_2","0_4_1_3","0_4_1_4","1_0_0_0","1_0_0_1","1_0_0_2","1_0_0_3","1_0_0_4","1_1_0_0","1_1_0_1","1_1_0_2","1_1_0_3","1_1_0_4","1_2_0_0","1_2_0_1","1_2_0_2","1_2_0_3","1_2_0_4", "1_3_0_0","1_3_0_1","1_3_0_2","1_3_0_3","1_3_0_4", "1_4_0_0","1_4_0_1","1_4_0_2","1_4_0_3","1_4_0_4"}; //TODO automatise

  // Generate dummy data
  TheoryFcn2 f_dummy(empty_data, empty_error, labels);
  vector<vector<double>> dummy_call = f_dummy.DummyData(n_data_points, error_start);
  vector<double> data = dummy_call[0];
  vector<double> error = dummy_call[1];
  cout << data << "\n" << error << "\n";

  TheoryFcn fFCN(data,error,labels);

  // Create parameters with initial starting values
  vector<double> par_error(n_parameters, 1.0); // Will store error on parameter before taking scaling into account
  vector<double> start(n_parameters, 0.0); // TODO if these are consts for every parameter, no need to make this a vector, but if we want to iterate fit is good, does minuit not do it anyway?

  MnUserParameters upar;
  for(unsigned int i=0 ; i<n_parameters; ++i){
    upar.Add(Form("param%d",i), start[i], par_error[i]);
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
  // Get start time
  struct timeval tv_start, tv_stop;
  gettimeofday(&tv_start, 0);

  FunctionMinimum min = minimize(maxfcn, tolerance);  

  gettimeofday(&tv_stop, 0);
  t1 = (tv_stop.tv_sec - tv_start.tv_sec)*1000.0 + (tv_stop.tv_usec - tv_start.tv_usec)/1000.0;

  cout << "# Calculation time: " << t1 << " ms" << "\n";

  cout << "CHI^2: " << min.Fval() << " , chi^2/ndf: " << min.Fval()/(n_data_points - n_parameters) << "\n" << "\n";

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
    cout <<"par[" << i << "]: " << get<0>(getParameterNameAndScaling(i)) << " fitted to: "<< min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i)) << " +/- " << min.UserState().Error(i) * abs(get<1>(getParameterNameAndScaling(i))) << "\n";
  }

  return 0;
}


