// Simple minimizer for mass fits using interface provided by Math::Minimzer

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/IntegratorMultiDim.h"

std::unique_ptr<TFile> inputFile( TFile::Open("mass_fits_control_histos_smear_beta_val.root") );
std::unique_ptr<TH1D> scale(inputFile->Get<TH1D>("epsilon"));

const int number_4D_bins = 12; // TODO scale->GetEntries();
const int n_eta_bins = 2, n_pt_bins = 5;
const int n_parameters = n_eta_bins*3; // for A,e,M model
//pT binning must match the one used to make the histos
const double pt_binning[n_pt_bins+1] = {25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0}; //{38.0, 39.0, 40.0, 41.0, 42.0, 43.0};
const double scaling_A = 0.001, scaling_e = 0.001 * 40.0, scaling_M = 0.001 / 40.0; // this scaling makes fitter parameters of order 1

//-----------------------------------------------
// Function to get average k from pT bin index

double getK(int index){
  // TODO what about errors on k?
  // TODO could add the pT histo per eta,pt,eta,pt to get more accurate average pT
  return 2.0 / (pt_binning[index] + pt_binning[index+1]); // k = 1 / avg_pT
}

//-----------------------------------------------
// Function to get bin index from label name

vector<int> getIndices(string s){
  vector<int> indices;
  int pos = 0;
  for (int i=0; i<4-1;i++){ //4 for 4D binning
    pos = s.find("_");
    indices.push_back(stoi(s.substr(0,pos))); 
    s.erase(0,pos+1); // 1 is the length of the delimiter, "_"
  }
  indices.push_back(stoi(s));
  return indices;
}

//-----------------------------------------------
// Function to get parameter name (A0, e1, etc.) from index in parameters array

std::tuple<string,double> getParameterNameAndScaling(int index){
  
  int whole = index / n_eta_bins; 
  int rest = index % n_eta_bins;
  double scaling;
  
  string name;
  if (whole == 0) { 
    name = "A" + to_string(rest);
    scaling = scaling_A;
  }
  else if (whole == 1) {
    name = "e" + to_string(rest);
    scaling = scaling_e;
  }
  else if (whole == 2) {
    name = "M" + to_string(rest);
    scaling = scaling_M;
  }
  else { 
    std::cout<<"\n"<<"ERROR counting parameters"<<"\n";
  }
  
  return make_tuple(name, scaling); 
}


//-----------------------------------------------

double k_plus_values[number_4D_bins] = {0.025974, 0.025974, 0.025974, 0.025974, 0.0253165, 0.0253165, 0.0253165, 0.0253165, 0.0246914, 0.0246914, 0.0246914, 0.0246914};//{1.0/38.5, 1.0/45.2, 1.0/40.5, 1.0/44.2, 1.0/39.8, 1.0/37.4, 1.0/33.9, 1.0/46.2, 1.0/34.9, 1.0/45.9, 1.0/34.2, 1.0/38.4};
TVectorD k_plus_test(number_4D_bins, k_plus_values);

double k_minus_values[number_4D_bins] = {0.0253165, 0.0246914, 0.0240964, 0.0235294, 0.025974, 0.0246914, 0.0240964, 0.0235294, 0.025974, 0.0253165, 0.0240964, 0.0235294};//{1.0/35.9, 1.0/39.6, 1.0/42.0, 1.0/44.8, 1.0/41.0, 1.0/45.2, 1.0/39.1, 1.0/38.2, 1.0/35.0, 1.0/40.2, 1.0/46.9, 1.0/34.8};
TVectorD k_minus_test(number_4D_bins, k_minus_values);

//double scale_squared_values[number_4D_bins] = {0.993*0.993, 0.996*0.996, 0.992*0.992, 0.996*0.996, 0.998*0.998, 0.990*0.990, 0.992*0.992, 0.995*0.995, 0.994*0.994, 0.998*0.998, 0.996*0.996, 0.990*0.990}; 
//double scale_squared_values[number_4D_bins] = {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000};
double scale_squared_values[number_4D_bins] = {0.989709, 0.988419, 0.98794, 0.988019, 0.991833, 0.986387, 0.985908, 0.985987, 0.992988, 0.988824, 0.987057, 0.987135};
TVectorD scale_squared(number_4D_bins, scale_squared_values);

double scale_squared_error_values[number_4D_bins] = {2*0.993*0.0005, 2*0.996*0.0005, 2*0.992*0.0005, 2*0.996*0.0005, 2*0.998*0.0005, 2*0.990*0.0005, 2*0.992*0.0005, 2*0.995*0.0005, 2*0.994*0.0005, 2*0.998*0.0005, 2*0.996*0.0005, 2*0.990*0.0005}; 
//double scale_squared_error_values[number_4D_bins] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
TVectorD scale_squared_error(number_4D_bins, scale_squared_error_values); //here the scale is squared, not the error

string labels[number_4D_bins] = {"0_0_1_1", "0_0_1_2", "0_0_1_3", "0_0_1_4", "0_1_1_0", "0_1_1_2", "0_1_1_3", "0_1_1_4", "0_2_1_0", "0_2_1_1", "0_2_1_3", "0_2_1_4"};  

//-----------------------------------------------
// Define chi2 function to minimise

double chi2(const double *parameters){
  TVectorD params(n_parameters, parameters);
  // params is a 1D array of size n_parameters = 3*number_eta_bins, idices from 0 to n-1 contain A, from n to 2n-1 epsilon, from 2n to 3n-1 M

  double k_plus, k_minus, val_i, val_i_squared, val = 0.0;
  int eta_pos_index, pt_pos_index, eta_neg_index, pt_neg_index, index_A_plus, index_e_plus, index_M_plus, index_A_minus, index_e_minus, index_M_minus;
  vector<int> bin_indices(4); // for 4D binning
  
  for (int i=0; i<number_4D_bins; i++){
    val_i = 0.0;

    bin_indices = getIndices(labels[i]);
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]);
    k_minus = getK(bin_indices[3]);
    //k_plus = k_plus_test[i];
    //k_minus = k_minus_test[i];

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of e,M
    val_i = ( (1. + scaling_A*params[eta_pos_index] - scaling_e*params[eta_pos_index + n_eta_bins]*k_plus +  scaling_M*params[eta_pos_index + 2*n_eta_bins]/k_plus)*(1. + scaling_A*params[eta_neg_index] - scaling_e*params[eta_neg_index + n_eta_bins]*k_minus - scaling_M*params[eta_neg_index + 2*n_eta_bins]/k_minus) - scale_squared[i] ) / scale_squared_error[i];
    val_i_squared = val_i*val_i;
 
    val += val_i_squared;

    vector<double> fake_data(number_4D_bins);
    fake_data[i] = (1. -0.095517 - (-2.20665)*k_plus + 0.00179058/k_plus)*(1. -0.0950043 - (-1.21125)*k_minus - (-0.000555938)/k_minus);
    cout<<fake_data[i]<<" ";
    
    /*
    A0 fitted to: -0.095517
      par[1]: A1 fitted to: -0.0950043
      par[2]: e0 fitted to: -2.20665
      par[3]: e1 fitted to: -1.21125
      par[4]: M0 fitted to: 0.00179058
      par[5]: M1 -0.000555938
      */
  }
  
  cout<<"\n";
  return val;
}

//-----------------------------------------------
// Minimizer

void constants_fitter(const char * minName = "Minuit", const char *algoName = "Migrad"){

  // Get data
  
  //for(int i=0; i<number_4D_bins; i++){
  //scale_squared_values[i] = (1.+ (scale->GetBinContent(i+1)))*(1. + (scale->GetBinContent(i+1))); // TODO that 1.+ needs changed, change that we fit for beta in mass
  //scale_squared_error_values[i] = 2*abs(scale->GetBinContent(i+1))*(scale->GetBinError(i+1));
  //labels[i] = scale->GetXaxis()->GetLabels()->At(i)->GetName();
  //std::cout << "\n" << labels[i] << " scale_squared_values " << scale_squared_values[i] << " scale_squared_error_values " << scale_squared_error_values[i] << "\n";
  //}

  // Initiate the minimizer
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  minimum->SetMaxFunctionCalls(1000000); 
  minimum->SetTolerance(0.001);
  minimum->SetPrintLevel(1);

  // TODO analytical gradient to improve performance

  // Set up function and variables
  ROOT::Math::Functor f( &chi2, n_parameters);
  minimum->SetFunction(f);

  double step[n_parameters] = {};
  for(unsigned int i=0 ; i<n_parameters; ++i) step[i] += 0.01; // TODO what step to use? why +1?
  double start[n_parameters] = {}; // initial values of parameters to fit, all are set to 0.
  for(unsigned int i=0 ; i<n_parameters; ++i){
    minimum->SetVariable(i, Form("param%d",i), start[i], step[i]);
  }

  // Minimization starts
  minimum->Minimize();

  const double *parameters_fitted = minimum->X();
  
  std::cout << std::endl << "Fitted A,e,M parameters: " << std::endl;
  for(unsigned int i=0 ; i<n_parameters; ++i) std::cout << get<0>(getParameterNameAndScaling(i)) <<": " << parameters_fitted[i] * get<1>(getParameterNameAndScaling(i)) << ", ";
  std::cout << std::endl << "minimum chi2: " << minimum->MinValue() << "; minimum chi2/ndf: " << (minimum->MinValue()) / (number_4D_bins - n_parameters) << std::endl;

  return;
}
