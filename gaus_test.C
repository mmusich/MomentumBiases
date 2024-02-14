// Standalone code to check gaussianity of pull distribution

#include "TFile.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TVectorT.h"
#include <string.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace ROOT;
using namespace ROOT::VecOps;

int main(){

  // define histogram with all pull values
  TH1D *pull_epsilon_control = new TH1D("pull_epsilon_control", " test: ((mean_yellow - mean_green) - epsilon) / error_epsilon ", 6, -3.0, 3.0);
  pull_epsilon_control->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control->GetYaxis()->SetTitle("Events");


  // initialise average pull
  double average_pull = 0.0;

  double mean = 0.0;
  double width = 1.0;
  double m_diff; 

  TRandom3* rand;
  TRandom3 rand_val;
  
  rand = &rand_val;
  // repeat toy for M times
  for(int M=1; M<=10; M++){
    //generate N random numbers between -5 and 5 from gaus(0, sigma) -> green
    rand_val = TRandom3(9361 + M);
    
    m_diff = rand->Gaus(mean, width);
    std::cout<< m_diff<<" ";

    // fill a histogram based on this

    // shift this using the weight trick -> yellow
    // make note of bias/sigma
    // fill a histogram based on this
    // fill variance matrix based on error

    // fill a vector green - yellow

    // calculate jacobian elements using green distribution (the generated numbers and their square)
    // fill a histogram based on this
    // fill a vector with the jacobian

    // define minimisation problem and solve for epsilon 
    // calculate pull as fitted - true / error_fitted, fill a histogram outside the for loop with this value
  }
   
  
  // calculate average pull
  // print & fit pull histogram, should be a normal gaussian 

  // repeat everything with a different N until pull histogram looks gaussian
  // additionally, try playing with the binning
}
