#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/IntegratorMultiDim.h"

std::unique_ptr<TFile> inputFile( TFile::Open("mass_fits_control_histos_smear_beta_val.root") );
std::unique_ptr<TH1D> scale(inputFile->Get<TH1D>("epsilon"));

const int number_4D_bins = 6; // TODO scale->GetEntries();
const int n_eta_bins = 2, n_pt_bins = 5;
const int n_parameters = n_eta_bins*3; // for A,e,M model
//pT binning must match the one used to make the histos
const double pt_binning[n_pt_bins+1] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; // TODO after debugging, scaling {25.0, 33.3584, 38.4562, 42.2942, 45.9469, 55.0};

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
  while(pos < s.size()){
    pos = s.find("_");
    indices.push_back(stoi(s.substr(0,pos))); 
    s.erase(0,pos+1); // 1 is the length of the delimiter, "_"
  }
  return indices;
}

//-----------------------------------------------
// Function to get parameter name (A0, e1, etc.) from index in parameters array

string getParameterName(int index){
  
  int whole = index / n_eta_bins; 
  int rest = index % n_eta_bins;
  
  string name;
  if (whole == 0) { name = "A" + to_string(rest); }
  else if (whole == 1) { name = "e" + to_string(rest); }
  else if (whole == 2) { name = "M" + to_string(rest); }
  else { std::cout<<"\n"<<"ERROR counting parameters"<<"\n"; }
  
  return name; 
}


//-----------------------------------------------

float k_plus_values[number_4D_bins] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
TVectorF k_plus_test(number_4D_bins, k_plus_values);

float k_minus_values[number_4D_bins] = {1.1, 1.1, 1.1, 1.1, 1.1, 1.1};
TVectorF k_minus_test(number_4D_bins, k_minus_values);

double scale_squared_values[number_4D_bins] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
TVectorD scale_squared(number_4D_bins, scale_squared_values);

double scale_squared_error_values[number_4D_bins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; 
TVectorD scale_squared_error(number_4D_bins, scale_squared_error_values); //here the scale is squared, not the error

string labels[number_4D_bins] = {"0_0_1_0", "0_0_1_0", "0_0_1_0", "0_0_1_0", "0_0_1_0", "0_0_1_0"};  

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
    //k_plus = getK(bin_indices[1]);
    //k_minus = getK(bin_indices[3]);
    k_plus = k_plus_test[i];
    k_minus = k_minus_test[i];

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k)
    val_i = ( (1. + params[eta_pos_index] - params[eta_pos_index + n_eta_bins]*k_plus + params[eta_pos_index + 2*n_eta_bins]/k_plus)*(1. + params[eta_neg_index] - params[eta_neg_index + n_eta_bins]*k_minus - params[eta_neg_index + 2*n_eta_bins]/k_minus) - scale_squared[i] ) / scale_squared_error[i];
    val_i_squared = val_i*val_i;
 
    val += val_i_squared;
  }

  return val;
}

//-----------------------------------------------
// Minimizer

void constants_fitter(const char * minName = "Minuit", const char *algoName = "Migrad"){

  // Get data
  //int j = 0, k = 1; //for testing
  //for(int i=0; i<number_4D_bins; i++){
  //scale_squared_values[i] = (1.+ (scale->GetBinContent(i+1)))*(1. + (scale->GetBinContent(i+1))); // TODO that 1.+ needs changed
  // scale_squared_error_values[i] = 2*abs(scale->GetBinContent(i+1))*(scale->GetBinError(i+1));
    //for testing
    //scale_squared_error_values[i] = 1.;
    //labels[i] = scale->GetXaxis()->GetLabels()->At(i)->GetName();
    //for testing
    //labels[i] = "0_"+to_string(j)+"_0_"+to_string(k);
    //j++;
  // k++;
  //if (j == n_pt_bins || k == n_pt_bins){ k = 0; j = 1;}
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
  for(unsigned int i=0 ; i<n_parameters; ++i) step[i] += 0.01;
  double start[n_parameters] = {}; // initial values of parameters to fit, all are set to 0.
  for(unsigned int i=0 ; i<n_parameters; ++i){
    minimum->SetVariable(i, Form("param%d",i), start[i], step[i]);
  }

  // Minimization starts
  minimum->Minimize();

  const double *parameters_fitted = minimum->X();
  
  std::cout << std::endl << "Fitted A,e,M parameters: " << std::endl;
  for(unsigned int i=0 ; i<n_parameters; ++i) std::cout << getParameterName(i) <<": " << parameters_fitted[i] << ", ";

  std::cout << std::endl << "minimum chi2: " << minimum->MinValue() << std::endl;

  return;
}
