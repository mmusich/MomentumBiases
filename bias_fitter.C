// Standalone code to fit for biases in mass distribtions per eta, pT, eta, pT bins

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TVectorT.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"
#include <string.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <sys/time.h>

using namespace ROOT;
using namespace ROOT::VecOps;
using namespace ROOT::Minuit2;
using namespace std;
using ROOT::RDF::RNode;

using Eigen::MatrixXd;
using Eigen::VectorXd;

const int nbinseta = 24; 
const int nbinspt = 5;

//--------------------------------------------------------
// Functions for TTree
//--------------------------------------------------------

//dxy_significance
RVecF dxy_significance(RVecF Muon_dxy, RVecF Muon_dxyErr){
  return abs(Muon_dxy)/Muon_dxyErr;
}

// MuonisGood
RVecB MuonisGood(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVec<UChar_t> Muon_genPartFlav, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
    if (Muon_pt[i] > 10.0 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && Muon_genPartFlav[i]==1 && dxy_significance[i] < 4){
      muonisgood.push_back(1);
    } else {
      muonisgood.push_back(0);
    }
  }
  return muonisgood;
}

// MuonisGoodData
RVecB MuonisGoodData(RVecF Muon_pt, RVecF Muon_eta, RVecB Muon_isGlobal, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVecF dxy_significance){
  RVecB muonisgood;
  for(int i=0;i<Muon_pt.size();i++){
    if (Muon_pt[i] > 10.0 && abs(Muon_eta[i]) < 2.4 && Muon_isGlobal[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15 && dxy_significance[i] < 4){
      muonisgood.push_back(1);
    } else {
      muonisgood.push_back(0);
    }
  }
  return muonisgood;
}

//--------------------------------------------------------
// Functions for histograms
//--------------------------------------------------------

// Functions for label names
string stringify_name(int a, int b, int c, int d){
  string hyp = "_";
  return to_string(a) + hyp + to_string(b) + hyp + to_string(c) + hyp + to_string(d);
}

string stringify_title(int a, int b, int c, int d, vector<double> e, vector<double> p){
  string comma = ",";
  string txt1 = "#eta+ in [", txt2="] pT+ in [", txt3="] #eta- in [", txt4="] pT- in [", txt5 = "]";
  return txt1+to_string(e[a]).substr(0,4)+comma+to_string(e[a+1]).substr(0,4)+txt2+to_string(p[b]).substr(0,4)+comma+to_string(p[b+1]).substr(0,4)+txt3+to_string(e[c]).substr(0,4)+comma+to_string(e[c+1]).substr(0,4)+txt4+to_string(p[d]).substr(0,4)+comma+to_string(p[d+1]).substr(0,4)+txt5;
}

// Function for Gaussian fit and integral 
vector<double> fitHisto(TH1* histogram, int draw_option, int color, int nsigma){

  vector<double> fitresult;

  double mean = histogram->GetMean();
  double sigma = histogram->GetRMS();
  double mean_err, sigma_err, integral, chi2;

  TF1 *gaussianFunc = new TF1("gaussianFunc", "gaus", mean - 5 * sigma, mean + 5 * sigma);
  // first fit
  if(0 == histogram->Fit(gaussianFunc, "QNR")){
    mean = gaussianFunc->GetParameter(1);
    sigma = gaussianFunc->GetParameter(2);

    // second fit: few sigma of first fit around mean of first fit
    gaussianFunc->SetRange(mean - nsigma * sigma, mean + nsigma * sigma);
    if (0 == histogram->Fit(gaussianFunc, "Q0R")) { // don't draw yet
      if (histogram->GetFunction(gaussianFunc->GetName())) { 
	histogram->GetFunction(gaussianFunc->GetName())->SetLineColor(color);
        if (draw_option == 1){ histogram->GetFunction(gaussianFunc->GetName())->ResetBit(TF1::kNotDraw);} // draw fit
      }
      mean = gaussianFunc->GetParameter(1);
      sigma = gaussianFunc->GetParameter(2);
      mean_err = gaussianFunc->GetParError(1);
      sigma_err = gaussianFunc->GetParError(2);
      chi2 = gaussianFunc->GetChisquare();

      TF1 *gaussianNewFunc = new TF1("gaussianNewFunc", "gaus", mean - 5 * sigma, mean + 5 * sigma);
      gaussianNewFunc->SetParameter(0, 1.0/(sigma*pow(2*M_PI,0.5)));
      gaussianNewFunc->SetParameter(1, mean);
      gaussianNewFunc->SetParameter(2, sigma);

      // integral of fitted gaussian function in range of the histogram
      integral = gaussianNewFunc->Integral(75.0, 105.0);
      // TODO implement line below
      //integral = gaussianNewFunc->Integral(histogram->GetXaxis()->GetBinLowEdge(1), histogram->GetXaxis()->GetBinLowEdge(histogram->GetNbins()+1)); // TODO histogram->GetNbins() includes under/over flow??

    } else {
      mean = -90.0;
      sigma = -5.0;
      mean_err = 90.0;
      sigma_err = 5.0;
      integral = -100.0;
      chi2 = -100.0;
    } 
  } else {
    mean = -90.0;
    sigma = -5.0;
    mean_err = 90.0;
    sigma_err = 5.0;
    integral = -100.0;
    chi2 = -100.0;
  }
  
  fitresult.push_back(mean);
  fitresult.push_back(sigma);
  fitresult.push_back(mean_err);
  fitresult.push_back(sigma_err);
  fitresult.push_back(integral);
  fitresult.push_back(chi2);
  
  return fitresult;
}

//--------------------------------------------------------
// Class and functions for minimization
//--------------------------------------------------------

class TheoryFcn : public FCNGradientBase {
//class TheoryFcn : public FCNBase {

public:

  TheoryFcn(const vector<double> &meas, const vector<double> &err, const vector<string> &bin_labels, const TMatrixD &V_phys_to_int, const TMatrixD &V_int_to_phys) : scaleSquared(meas), scaleSquaredError(err), binLabels(bin_labels), VPhysToInt(V_phys_to_int), VIntToPhys(V_int_to_phys), errorDef(1.0) {}
  ~TheoryFcn() {}

  virtual double Up() const {return errorDef;}
  virtual void SetErrorDef(double def) {errorDef = def;}

  virtual double operator()(const vector<double>&) const;
  virtual vector<double> Gradient(const vector<double>& ) const;
  virtual bool CheckGradient() const {return true;} 
  // virtual std::vector< double > ROOT::Minuit2::FCNGradientBase::Hessian(const std::vector< double > & )const -> allows to do Hessian analytically? 

  vector<int> getIndices(string bin_label) const; 
  double getK(const int pT_index) const;
  
  const vector<double> pT_binning {25.0, 33.2011, 38.3067, 42.2411, 46.055, 55.0}; // pT binning goes here, this is for 2016 //TODO add fct SetPtbinning
  static constexpr double scaling_A = 0.001, scaling_e = 0.001 * 40.0, scaling_M = 0.001 / 40.0, scaling_e_prime = 0.001 / 0.01; // this scaling makes fitter parameters of order 1 

  vector<string> binLabels; // made public for the closure test
  TMatrixD VPhysToInt, VIntToPhys;
  
private:

  vector<double> scaleSquared;
  vector<double> scaleSquaredError;
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
// Function to build reduced chi2 

double TheoryFcn::operator()(const vector<double>& par) const {
  // par has size n_parameters = 3*number_eta_bins, idices from 0 to n-1 contain A, from n to 2n-1 epsilon, from 2n to 3n-1 M
  assert(par.size() == 3*nbinseta); 

  double chi2(0.0);
  double diff(0.0);

  double k_plus, k_minus, term_pos, term_neg;  
  double my_func;
  int eta_pos_index, eta_neg_index;
  vector<int> bin_indices(4); // for 4D binning
  
  double k_middle = (1.0/pT_binning[nbinspt] + 1.0/pT_binning[0])/2.0;

  vector<double> physical_M_pars(nbinseta);
  TMatrixD internal_M_matrix(nbinseta,1), physical_M_matrix(nbinseta,1);
  TArrayD internal_M_elements(nbinseta);
  
  for(int i=0; i<nbinseta; i++){
    internal_M_elements[i] = par[2*nbinseta + i];
  }
  
  internal_M_matrix.SetMatrixArray(internal_M_elements.GetArray());
  physical_M_matrix = TMatrixD(VIntToPhys,TMatrixD::kMult,internal_M_matrix);
 
  //std::cout<< "\n" << "physical M array"<< "\n";
  for(int i=0; i<nbinseta; i++){
    physical_M_pars[i] = TMatrixDColumn(physical_M_matrix,0)(i);
    //std::cout<< " "<< physical_M_pars[i];
  }
  //std::cout<<"\n";

  for(unsigned int n(0); n < scaleSquared.size() ; n++) {
    bin_indices = getIndices(binLabels[n]); 
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]); 
    k_minus = getK(bin_indices[3]);

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of A,e,M
    //term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e*par[eta_pos_index + nbinseta]*k_plus +  scaling_M*par[eta_pos_index + 2*nbinseta]/k_plus);
    //term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e*par[eta_neg_index + nbinseta]*k_minus - scaling_M*par[eta_neg_index + 2*nbinseta]/k_minus);
    
    // decorrelate A, e by shifting origin of k
    // (1 + A'(+) - e'(+)k' + M(+)/k)(1 + A'(-) - e'(-)k' - M(-)/k) and scaling of A,e,M
    //par[eta_pos_index] now has the meaning of A' rather than A
    //par[eta_pos_index + nbinseta] now has the meaning of e' rather than e
    //term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e_prime*par[eta_pos_index + nbinseta]*(k_plus - k_middle) + scaling_M*par[eta_pos_index + 2*nbinseta]/k_plus);
    //term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e_prime*par[eta_neg_index + nbinseta]*(k_minus - k_middle) - scaling_M*par[eta_neg_index + 2*nbinseta]/k_minus);

    // decorrelate A, e by shifting origin of k and decorrelate M by fitting a linear combination of physical parameters 
    term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e_prime*par[eta_pos_index + nbinseta]*(k_plus - k_middle) + scaling_M*physical_M_pars[eta_pos_index]/k_plus);
    term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e_prime*par[eta_neg_index + nbinseta]*(k_minus - k_middle) - scaling_M*physical_M_pars[eta_neg_index]/k_minus);

    my_func = term_pos*term_neg; 
    
    diff = my_func - scaleSquared[n]; 
    chi2 += diff*diff/(scaleSquaredError[n]*scaleSquaredError[n]);

  }
 
  int ndof = scaleSquared.size() - par.size();
 
  return (chi2/ndof - 1.0); // minimise reduced chi2 directly 

}

//-----------------------------------------------
// Function to build gradient of reduced chi2 analytically 

vector<double> TheoryFcn::Gradient(const vector<double> &par ) const {

  assert(par.size() == 3*nbinseta);

  vector<double> grad(par.size(),0.0);
  TArrayD grad_wrt_M_physical_array(nbinseta), grad_wrt_M_internal_array(nbinseta); // TODO and if not eta 24 bins?
  TMatrixD grad_wrt_M_physical_matrix(nbinseta,1), grad_wrt_M_internal_matrix(nbinseta,1);
  
  double temp(0.0), local_func(0.0), local_grad(0.0);

  double k_plus, k_minus, term_pos, term_neg;
  double k_middle = (1.0/pT_binning[nbinspt] + 1.0/pT_binning[0])/2.0;
  
  int eta_pos_index, eta_neg_index;
  vector<int> bin_indices(4); // for 4D binning
  
  int ndof = scaleSquared.size() - par.size();

  vector<double> physical_M_pars(nbinseta);
  TMatrixD internal_M_matrix(nbinseta,1), physical_M_matrix(nbinseta,1);
  TArrayD internal_M_elements(nbinseta);
  
  for(int i=0; i<nbinseta; i++){
    internal_M_elements[i] = par[2*nbinseta + i];
  }
    
  internal_M_matrix.SetMatrixArray(internal_M_elements.GetArray());
  physical_M_matrix = TMatrixD(VIntToPhys,TMatrixD::kMult,internal_M_matrix);
  
  for(int i=0; i<nbinseta; i++){
    physical_M_pars[i] = TMatrixDColumn(physical_M_matrix,0)(i);
  }
  
  for(unsigned int n(0); n < scaleSquared.size() ; n++) { // loops over measurements
    // TODO I need everything from the chi2 function, maybe better solution to this
    bin_indices = getIndices(binLabels[n]); 
    eta_pos_index = bin_indices[0];
    eta_neg_index = bin_indices[2];

    // get k value and TODO error from pT bin index
    k_plus = getK(bin_indices[1]); 
    k_minus = getK(bin_indices[3]);

    // (1 + A(+) - e(+)k + M(+)/k)(1 + A(-) - e(-)k - M(-)/k) and scaling of e,M
    //term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e*par[eta_pos_index + nbinseta]*k_plus +  scaling_M*par[eta_pos_index + 2*nbinseta]/k_plus);
    //term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e*par[eta_neg_index + nbinseta]*k_minus - scaling_M*par[eta_neg_index + 2*nbinseta]/k_minus);
    
    // decorrelate A, e by shifting origin of k
    // (1 + A'(+) - e'(+)k' + M(+)/k)(1 + A'(-) - e'(-)k' - M(-)/k) and scaling of A,e,M
    //par[eta_pos_index] now has the meaning of A' rather than A
    //par[eta_pos_index + nbinseta] now has the meaning of e' rather than e
    //term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e_prime*par[eta_pos_index + nbinseta]*(k_plus - k_middle) + scaling_M*par[eta_pos_index + 2*nbinseta]/k_plus);
    //term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e_prime*par[eta_neg_index + nbinseta]*(k_minus - k_middle) - scaling_M*par[eta_neg_index + 2*nbinseta]/k_minus);

    // decorrelate A, e by shifting origin of k and decorrelate M by fitting a linear combination of physical parameters
    term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e_prime*par[eta_pos_index + nbinseta]*(k_plus - k_middle) + scaling_M*physical_M_pars[eta_pos_index]/k_plus);
    term_neg = (1. + scaling_A*par[eta_neg_index] - scaling_e_prime*par[eta_neg_index + nbinseta]*(k_minus - k_middle) - scaling_M*physical_M_pars[eta_neg_index]/k_minus);

    local_func = term_pos*term_neg;
    
    temp=2.0*(local_func - scaleSquared[n])/(scaleSquaredError[n]*scaleSquaredError[n])/ndof; 
    
    // derivative wrt A'(+), which is fpar[eta_pos_index]
    local_grad = scaling_A*term_neg;
    grad[eta_pos_index] += temp*local_grad;
    
    // derivative wrt e'(+)
    local_grad = -scaling_e_prime*(k_plus-k_middle)*term_neg;
    grad[eta_pos_index + nbinseta] += temp*local_grad;

    // derivative wrt M(+), next transform to derivative wrt M_internal(+)
    local_grad = scaling_M/k_plus*term_neg;
    grad[eta_pos_index + 2*nbinseta] += temp*local_grad;
    
    // derivative wrt A'(-)
    local_grad = scaling_A*term_pos;
    grad[eta_neg_index] += temp*local_grad;

    // derivative wrt e'(-)
    local_grad = -scaling_e_prime*(k_minus-k_middle)*term_pos;
    grad[eta_neg_index + nbinseta] += temp*local_grad;

    // derivative wrt M(-), next transform to derivative wrt M_internal(-)
    local_grad = -scaling_M/k_minus*term_pos;
    grad[eta_neg_index + 2*nbinseta] += temp*local_grad;

  }

  // -----------------------------------------------
  // Transform M terms to derivatives wrt M_internal
  for(int i=0; i<nbinseta; i++){
    grad_wrt_M_physical_array[i] = grad[2*nbinseta + i];
  }
  
  grad_wrt_M_physical_matrix.SetMatrixArray(grad_wrt_M_physical_array.GetArray());
  grad_wrt_M_internal_matrix = TMatrixD(VPhysToInt,TMatrixD::kMult,grad_wrt_M_physical_matrix);

  for(int i=0; i<nbinseta; i++){
    grad[2*nbinseta + i] = TMatrixDColumn(grad_wrt_M_internal_matrix,0)(i);
  }
  // -----------------------------------------------
  
  return grad;
}

//-----------------------------------------------
// Function to get parameter name (A0, e1, etc.) and appropriate scaling from index in parameters array

tuple<string,double,string> getParameterNameAndScaling(int index){

  int whole = index / nbinseta;
  int rest = index % nbinseta;
  double scaling;

  string name, physical_name;
  if (whole == 0) {
    name = "A_prime" + to_string(rest);
    physical_name = "A" + to_string(rest);
    scaling = TheoryFcn::scaling_A;
  }
  else if (whole == 1) {
    name = "e_prime" + to_string(rest);
    physical_name = "e" + to_string(rest);
    scaling = TheoryFcn::scaling_e_prime;
  }
  else if (whole == 2) {
    name = "M_internal" + to_string(rest);
    physical_name = "M" + to_string(rest);
    scaling = TheoryFcn::scaling_M;
  }
  else {
    cout<<"\n"<<"ERROR counting parameters"<<"\n";
  }

  return make_tuple(name, scaling, physical_name);
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


//--------------------------------------------------------
// Main
//--------------------------------------------------------
int bias_fitter(){

  double t1(0.);
  // Get start time
  struct timeval tv_start, tv_stop;
  gettimeofday(&tv_start, 0);

  ROOT::EnableImplicitMT();

  // Choose validation/analysis mode
  //string mode_option("simulation"); //TODO pass as command line argument

  //TODO make it so that label of data sample is attached to output files

  //--------------------------------------------------------
  // Define dataframe
  
  TChain chain("Events");
  //chain.Add("/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_1.root");
  //chain.Add("/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0000/NanoV9MCPostVFP_1.root");
  //chain.Add("/scratch/wmass/y2017/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2017_TrackFitV722_NanoProdv3/NanoV9MC2017_1.root");
  
  string line;
  
  ifstream file("InOutputFiles/MCFilenames_2016.txt");
  if (file.is_open()) {
    while (getline(file, line)) {
      chain.Add(line.c_str());
    }
    file.close();
  }
  
  RDataFrame df(chain); // TODO write for analysis mode, will need 2 data frames, but honestly should be a different script that people use

  auto d0 = std::make_unique<RNode>(df);

  // Random numbers for pairs
  unsigned int nslots = d0->GetNSlots();
  std::cout<<"1st frame: nslots = "<<nslots<<"\n";
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(4357 + i*10) );
  }
  
  d0 = std::make_unique<RNode>(d0->Filter("HLT_IsoMu24 == 1"));
  d0 = std::make_unique<RNode>(d0->Filter("nMuon >= 2"));
  d0 = std::make_unique<RNode>(d0->Filter("PV_npvsGood >= 1"));
  
  //if (mode_option.compare("simulation") == 0) {
  d0 = std::make_unique<RNode>(d0->Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)"));
  d0 = std::make_unique<RNode>(d0->Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, Muon_genPartFlav, dxy_significance)"));
  //} else if (mode_option.compare("data") == 0){
    //d0 = std::make_unique<RNode>(d0->Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)"));
    //d0 = std::make_unique<RNode>(d0->Define("MuonisGood", "MuonisGoodData(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, dxy_significance)"));
    //}
  
  // Binning (pT done later)
  double pt_low = 25.0, pt_high = 55.0, eta_low = -2.4, eta_high = 2.4, mll_diff_low = -7.0,  mll_diff_high = 7.0, mll_low = 75.0, mll_high = 105.0;
  int nbinsmll_diff_over_gen=32, nbinsmll_diff=22, nbinsmll=32;
  vector<double> etabinranges, ptbinranges, mllbinranges;
  vector<double> mll_diffbinranges, mll_diff_over_genbinranges;
  vector<float> A_values(nbinseta), e_values(nbinseta), M_values(nbinseta);

  std::cout<<"\n etabinranges = [";
  for (int i=0; i<=nbinseta; i++){
    etabinranges.push_back(eta_low + i * (eta_high - eta_low)/nbinseta);
    std::cout<<etabinranges[i]<<", ";
  }
  std::cout<<"] \n";

  double myetaboundaries[etabinranges.size()];
  for (int i=0; i<etabinranges.size(); i++){
    myetaboundaries[i] = etabinranges[i];
  }

  for (int i=0; i<=nbinsmll_diff; i++) mll_diffbinranges.push_back(mll_diff_low + i*(mll_diff_high - mll_diff_low)/nbinsmll_diff);
  double mymll_diffboundaries[mll_diffbinranges.size()];
  for (int i=0; i<mll_diffbinranges.size(); i++){
    mymll_diffboundaries[i] = mll_diffbinranges[i];
  }

  
  for (int i=0; i<=nbinsmll; i++) mllbinranges.push_back(mll_low + i*(mll_high-mll_low)/nbinsmll);
  double mymllboundaries[mllbinranges.size()];
  for (int i=0; i<mllbinranges.size(); i++){
    mymllboundaries[i] = mllbinranges[i];
  }
  for (int i=0; i<=nbinsmll_diff_over_gen; i++) mll_diff_over_genbinranges.push_back(-0.1 + i*2.0*0.1/nbinsmll_diff_over_gen);

  int nbins = nbinseta*nbinseta*nbinspt*nbinspt;
  int nbins_plus_one = nbins+1;
  double mybinsboundaries[nbins_plus_one];
  for (int i=0; i<=nbins; i++){
    mybinsboundaries[i] = i;
  }

  float kmean_val = 0.5*( 1./25. + 1./55. ); // TODO don't input pt bounds by hand
  TRandom3* ran0 = new TRandom3(932);
  
  for (int i=0; i<nbinseta; i++){
    //std::cout<< "\n" << "i: " << i << "\n";
    // A from -0.0002 to 0.0005 and back
    A_values[i] = (-7.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 5.0)*0.0001;
    std::cout<< i << " A: "<<A_values[i] << " ";
    //A_values[i] = -0.001*( 1 + ((etabinranges[i]+etabinranges[i+1])/2.0/2.4)*((etabinranges[i]+etabinranges[i+1])/2.0/2.4) );
    //A_values[i] = 0.0004; 
    //std::cout<< "A: " << A_values[i] << ", ";
    // e from 0.01 to 0.001 and back
    e_values[i] = (9.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 1.0)*0.001;
    std::cout<< i << " e: "<<e_values[i] << " ";
    //e_values[i] = ran0->Uniform(-0.0001/kmean_val, 0.0001/kmean_val);
    //e_values[i] = 0.002; 
    //std::cout<< "e: " << e_values[i] << ", ";
    // M from 4*10^-5 to -2*10^-5 and back
    M_values[i] = (6.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) - 2.0)*0.00001;
    std::cout<< i << " M: "<<M_values[i] << "\n";
    //M_values[i] = ran0->Uniform(-0.001*kmean_val, 0.001*kmean_val);
    //M_values[i] = 0.00001; 
    //std::cout<< "M: " << M_values[i] << ", ";
  }

  auto GetEtaBin = [&](float eta)->int{
    for (int i=0; i<nbinseta; i++){
      if( etabinranges[i] <= eta && eta < etabinranges[i+1]) { return i; }
    }
    return -9;
  };

  //Pairs
  auto pairs = [&](unsigned int nslot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi)->std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,int,int>{

    // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 mll_smear_beta_val, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
    RVec<std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,int,int>> pairs; 
    std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,int,int> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco, firstPt_smear, secondPt_smear, mll_smear, firstPt_gen, secondPt_gen, mll_gen, firstPt_smear_beta_val, secondPt_smear_beta_val, smear_beta_weight;
    float smear_pt, mean, width, smeared_mean, smeared_curvature, smear_beta_weight_first_term, smear_beta_weight_second_term;
    //float A = 0.0004, e = 0.002, M = 0.00001; //for now constant per eta bin
    float A, e, M;
    int eta_bin_idx, first_gen_idx, second_gen_idx;

    for(int i=1;i<Muon_pt.size();i++){
      if(MuonisGood[i]){
	for(int j=0;j<i;j++){
	  if(MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	    TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen, firstSmearTrack, secondSmearTrack, motherSmear, firstSmearBiasTrack, secondSmearBiasTrack, motherSmearBias;
	    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	    mother = firstTrack + secondTrack;
	    mll_reco = mother.M();
	    firstPt_reco = Muon_pt[i];
	    secondPt_reco = Muon_pt[j];
	            
	    if(mll_low<mll_reco && mll_reco<mll_high){ //Cut in mll
	      //Gen match
	      bool firstGenMatched = false, secondGenMatched = false;
	      for (int k=0;k<GenPart_eta.size();k++){
		if (GenPart_status[k]==1 && abs(GenPart_pdgId[k])==13 && GenPart_pdgId[GenPart_genPartIdxMother[k]]==23){ // mu(-) has PDGID 13
		  if(pow(pow(GenPart_eta[k]-Muon_eta[i],2) + pow(TVector2::Phi_mpi_pi(GenPart_phi[k]-Muon_phi[i]),2),0.5)<0.3){
		    firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    firstPt_gen = GenPart_pt[k];
		    first_gen_idx = k;
		            
		    //smear 1st muon
		    mean = GenPart_pt[k]; //beta = 1 
		    width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		    firstPt_smear = rans[nslot]->Gaus(mean, width);
		    firstSmearTrack.SetPtEtaPhiM(firstPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    
		    //smear_beta_val, weight for beta != 1 use Eqn 22 in note AN2021_131_v5
		    eta_bin_idx = GetEtaBin(Muon_eta[i]);
		    if (eta_bin_idx < 0 || eta_bin_idx > 23){
		      std::cout<<"\n"<<"WARNING eta out of bounds";
		      pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
		      return pair_to_return ;
		    }
		    A = A_values[eta_bin_idx];
		    e = e_values[eta_bin_idx];
		    M = M_values[eta_bin_idx];
		    if (GenPart_pdgId[k]>0){ //mu(-)
		      if ( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) >= 0.0 ){
			smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) ,0.5) )/(-2.0)/e;
			if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
			  std::cout<<"\n"<<"WARNING WRONG CURVATURE";
			  pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
			  return pair_to_return ;
			}
		      } else {
			std::cout<<"\n"<<"WARNING negative delta";
			pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
			return pair_to_return ;
		      }
		      //smeared_curvature = 1.0/mean + M; // for M only
		    } else { //mu(+)
		      if ( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) >= 0.0 ){
			smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) ,0.5) )/(-2.0)/e;
			if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
			  std::cout<<"\n"<<"WARNING WRONG CURVATURE";
			  pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
			  return pair_to_return ;
			}
                      } else {
                        std::cout<<"\n"<<"WARNING negative delta";
                        pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
                        return pair_to_return ;
                      }
		      //smeared_curvature = 1.0/mean - M; // for M only
		    }
		    smeared_mean = 1.0/smeared_curvature;
		    firstPt_smear_beta_val = rans[nslot]->Gaus(smeared_mean, width);
                    firstSmearBiasTrack.SetPtEtaPhiM(firstPt_smear_beta_val, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_weight_first_term = TMath::Gaus(firstPt_smear, smeared_mean, width) / TMath::Gaus(firstPt_smear, mean, width);
		    		    
		    firstGenMatched = true;
		    if(secondGenMatched == true){break;}
		  } else if(pow(pow(GenPart_eta[k]-Muon_eta[j],2) + pow(TVector2::Phi_mpi_pi(GenPart_phi[k]-Muon_phi[j]),2),0.5)<0.3){
		    secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		    secondPt_gen = GenPart_pt[k];
		    second_gen_idx = k;
		    
		    //smear 2nd muon
		    mean = GenPart_pt[k]; //beta = 1
		    width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		    secondPt_smear = rans[nslot]->Gaus(mean, width);
		    secondSmearTrack.SetPtEtaPhiM(secondPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_val, weight for beta != 1
		    eta_bin_idx = GetEtaBin(Muon_eta[j]);
		    if (eta_bin_idx < 0 || eta_bin_idx > 23){
                      std::cout<<"\n"<<"WARNING eta out of bounds";
                      pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
                      return pair_to_return ;
                    }
                    A = A_values[eta_bin_idx];
                    e = e_values[eta_bin_idx];
                    M = M_values[eta_bin_idx];
		    if (GenPart_pdgId[k]>0){ //mu(-)
                      if ( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) >= 0.0 ){
                        smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(-1.0*M - 1.0/mean) ,0.5) )/(-2.0)/e;
                        if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
                          std::cout<<"\n"<<"WARNING WRONG CURVATURE";
                          pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
                          return pair_to_return ;
                        }
                      } else {
                        std::cout<<"\n"<<"WARNING negative delta";
                        pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
                        return pair_to_return ;
                      }
		      //smeared_curvature = 1.0/mean + M; // for M only
                    } else { //mu(+)
                      if ( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) >= 0.0 ){
                        smeared_curvature = ( -1.0*(1.0+A) + pow( (1.0+A)*(1.0+A) + 4.0*e*(M - 1.0/mean) ,0.5) )/(-2.0)/e;
                        if (smeared_curvature < 0.0 || std::isnan(smeared_curvature)){
                          std::cout<<"\n"<<"WARNING WRONG CURVATURE";
                          pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
                          return pair_to_return ;
                        }
                      } else {
                        std::cout<<"\n"<<"WARNING negative delta";
                        pair_to_return = make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
                        return pair_to_return ;
                      }
		      //smeared_curvature = 1.0/mean - M; // for M only
                    }
		    smeared_mean = 1.0/smeared_curvature;
		    secondPt_smear_beta_val = rans[nslot]->Gaus(smeared_mean, width);
                    secondSmearBiasTrack.SetPtEtaPhiM(secondPt_smear_beta_val, GenPart_eta[k], GenPart_phi[k], rest_mass);
		    //smear_beta_weight_second_term = TMath::Gaus(secondPt_smear, smeared_mean, width) / TMath::Gaus(secondPt_smear, mean, width); 
		    		    
		    secondGenMatched = true;
		    if(firstGenMatched == true){break;}
		  }
		}
	      }
	      if(firstGenMatched == false || secondGenMatched == false){
		continue;
	      }
	                  
	      motherSmear = firstSmearTrack + secondSmearTrack;
	      mll_smear = motherSmear.M();
	      //smear_beta_weight = smear_beta_weight_first_term * smear_beta_weight_second_term;

	      motherSmearBias = firstSmearBiasTrack + secondSmearBiasTrack;
              float mll_smear_beta_val = motherSmearBias.M();
	      
	      motherGen = firstGenTrack + secondGenTrack;
	      float mll_gen = motherGen.M();
	      //--------------------------------------------------------------------------------
              
	      // ATTENTION OVERWRITING smear
              // smear mass directly
              //mll_smear = rans[nslot]->Gaus(mll_gen, 0.02*mll_gen);
              //smear_beta_weight = TMath::Gaus(mll_smear, mll_gen*(1-0.001), 0.02*mll_gen) / TMath::Gaus(mll_smear, mll_gen, 0.02*mll_gen);
              
	      //--------------------------------------------------------------------------------
                    

	      if(Muon_charge[i]==1){
		// <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 mll_smear_beta_val, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val, 14 pos_gen_idx, 15 neg_gen_idx>
		temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco,mll_smear,firstPt_smear,secondPt_smear,mll_gen,firstPt_gen,secondPt_gen,mll_smear_beta_val,firstPt_smear_beta_val,secondPt_smear_beta_val,first_gen_idx, second_gen_idx);
	      } else {
		temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco,mll_smear,secondPt_smear,firstPt_smear,mll_gen,secondPt_gen,firstPt_gen,mll_smear_beta_val,secondPt_smear_beta_val,firstPt_smear_beta_val,second_gen_idx, first_gen_idx);
	      }
	      pairs.push_back(temp);
	    }
	  }
	}
      }
    }

    if(pairs.size()==1){
      pair_to_return=pairs.at(0);
    } else if(pairs.size()>1){
      float diff=100.0;
      int best=0;
      for(int i=0;i<pairs.size();i++){
	if(abs(get<2>(pairs.at(i))-91)<diff){
	  diff=(abs(get<2>(pairs.at(i))-91));
	  best=i;
	}
      }
      pair_to_return=pairs.at(best);
    } else {
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0);
    }
    return pair_to_return;
  };

  /*
  //PairsData
  auto pairsData = [&](RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecB MuonisGood, RVecF Muon_dxy, RVecF Muon_dz)->std::tuple<int,int,float,float,float>{
    
    // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco>
    RVec<std::tuple<int,int,float,float,float>> pairs; 
    std::tuple<int,int,float,float,float> temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco;
    
    for(int i=1;i<Muon_pt.size();i++){
      if(MuonisGood[i]){
	for(int j=0;j<i;j++){
	  if(MuonisGood[j] && Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){ //TODO better deltaR matching
	    TLorentzVector firstTrack, secondTrack, mother;
	    firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	    secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	    mother = firstTrack + secondTrack;
	    mll_reco = mother.M();
	    firstPt_reco = Muon_pt[i];
	    secondPt_reco = Muon_pt[j];
	            
	    if(mll_low<mll_reco && mll_reco<mll_high){ //Cut in mll
	      if(Muon_charge[i]==1){
		// <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco >
		temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco);
	      } else {
		temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco);
	      }
	      pairs.push_back(temp);
	    }
	  }
	}
      }
    }
    
    if(pairs.size()==1){
      pair_to_return=pairs.at(0);
    } else if(pairs.size()>1){
      float diff=100.0;
      int best=0;
      for(int i=0;i<pairs.size();i++){
	if(abs(get<2>(pairs.at(i))-91)<diff){
	  diff=(abs(get<2>(pairs.at(i))-91));
	  best=i;
	}
      }
      pair_to_return=pairs.at(best);
    } else {
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0);
    }
    return pair_to_return;
  };
  */
  
  //if (mode_option.compare("simulation") == 0) {
  d0 = std::make_unique<RNode>(d0->DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"}));
    //} else if (mode_option.compare("data") == 0){
    //d0 = std::make_unique<RNode>(d0->Define("pairs", pairsData, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz"}));
    //}
  
  d0 = std::make_unique<RNode>(d0->Define("mll_reco","return get<2>(pairs);")); 
  d0 = std::make_unique<RNode>(d0->Filter("mll_reco>1.0")); // this means only events with one mu pair are kept
    
  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
  
  d0 = std::make_unique<RNode>(d0->Define("posTrackPt","return get<3>(pairs);")); //TOOD make binning on smear not reco
  //dlast = std::make_unique<RNode>(dlast->Define("negTrackPt","return get<4>(pairs);"));
  //dlast = std::make_unique<RNode>(dlast->Define("posTrackEta","return Muon_eta[get<0>(pairs)];"));
  //dlast = std::make_unique<RNode>(dlast->Define("negTrackEta","return Muon_eta[get<1>(pairs)];"));
  ////dlast = std::make_unique<RNode>(dlast->Define("bin_index_reco", GetEtaPtEtaPtBin ,{"posTrackEta", "posTrackPt", "negTrackEta", "negTrackPt"}));

  ////if (mode_option.compare("simulation") == 0) { //NOTE ATTENTION can't yet run analysis mode with optimised pT binning
  d0 = std::make_unique<RNode>(d0->Define("weight", "std::copysign(1.0, genWeight)"));
  ////}
  
  //pT bin optimisation starts
  // Call 1st event loop
  std::cout<<"Call 1st event loop: pT binning optimization"<<"\n";
  auto pt_pos_uni = d0->Histo1D({"pt_pos_uni", "pt mu+", nbinspt*3, pt_low, pt_high},"posTrackPt","weight");
  
  // Get quartiles
  double xq[nbinspt+1], myptboundaries[nbinspt+1];
  for (int i=0;i<=nbinspt;i++) xq[i] = float(i)/nbinspt;
  pt_pos_uni->GetQuantiles(nbinspt+1,myptboundaries, xq);
  
  ptbinranges.clear();
  std::cout<<"ptbinranges = [";
  for (int i=0; i<=nbinspt; i++){
    ptbinranges.push_back(myptboundaries[i]);
    std::cout<<ptbinranges[i]<<", ";
  }
  std::cout<<"] \n";

  //---------------------------------------------------
  // Sort events in eta pt eta pt bins

  // GetEtaPtEtaPtBin
  auto GetEtaPtEtaPtBin = [&](float pos_eta, float pos_pt, float neg_eta, float neg_pt )->int{
    int idx = 0;
    for (int i=0; i<nbinseta; i++){
      for (int j=0; j<nbinspt; j++){
	for (int k=0; k<nbinseta; k++){
	  for (int l=0; l<nbinspt; l++){
	    if(etabinranges[i] <= pos_eta && pos_eta < etabinranges[i+1] && ptbinranges[j] <= pos_pt && pos_pt < ptbinranges[j+1] && etabinranges[k] <= neg_eta && neg_eta < etabinranges[k+1] && ptbinranges[l] <= neg_pt && neg_pt < ptbinranges[l+1] ) { return idx; } 
	    idx++;
	  }
	}
      }
    }
    return -9;
  };

  auto d1 = std::make_unique<RNode>(df);
  d1 = std::make_unique<RNode>(d1->Filter("event%2==0")); //MC like

  // Random numbers for pairs
  nslots = d1->GetNSlots();
  std::cout<<"2nd frame: nslots = "<<nslots<<"\n";
  for(unsigned int i = 0; i < nslots; i++){
    rans[i] = new TRandom3(8724 + i*10) ;
  }

  d1 = std::make_unique<RNode>(d1->Filter("HLT_IsoMu24 == 1"));
  d1 = std::make_unique<RNode>(d1->Filter("nMuon >= 2"));
  d1 = std::make_unique<RNode>(d1->Filter("PV_npvsGood >= 1"));
  
  d1 = std::make_unique<RNode>(d1->Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)"));
  d1 = std::make_unique<RNode>(d1->Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, Muon_genPartFlav, dxy_significance)"));
  
  d1 = std::make_unique<RNode>(d1->DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"}));
  
  d1 = std::make_unique<RNode>(d1->Define("mll_reco","return get<2>(pairs);"));
  d1 = std::make_unique<RNode>(d1->Filter("mll_reco>1.0")); // this means only events with one mu pair are kept
  
  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood
  
  // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
  
  d1 = std::make_unique<RNode>(d1->Define("posPtSmear","return get<6>(pairs);"));
  d1 = std::make_unique<RNode>(d1->Define("negPtSmear","return get<7>(pairs);"));
  d1 = std::make_unique<RNode>(d1->Define("posGenEta","return GenPart_eta[get<14>(pairs)];"));
  d1 = std::make_unique<RNode>(d1->Define("negGenEta","return GenPart_eta[get<15>(pairs)];"));
  d1 = std::make_unique<RNode>(d1->Define("bin_index_smear", GetEtaPtEtaPtBin ,{"posGenEta", "posPtSmear", "negGenEta", "negPtSmear"}));
  d1 = std::make_unique<RNode>(d1->Define("mll_smear","return get<5>(pairs);"));
  d1 = std::make_unique<RNode>(d1->Define("mll_gen","return get<8>(pairs);"));
  d1 = std::make_unique<RNode>(d1->Define("mll_diff_smear","return mll_smear - mll_gen;"));
  d1 = std::make_unique<RNode>(d1->Define("weight", "std::copysign(1.0, genWeight)"));
  
  // Take MC frame and fit diff smear per k bin, save results in a histo
  
  // MC diff
  // Call 2nd event loop
  std::cout<<"Call 2nd event loop: get MC resolution"<<"\n";
  auto mll_diff_smear_diffbin_4Dbin = d1->Histo2D({"mll_diff_smear_diffbin_4Dbin", "mll_diff_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"mll_diff_smear", "bin_index_smear", "weight");
  
  auto mDh_diff_mc_ptr = mll_diff_smear_diffbin_4Dbin;
  auto mDh_proj_diff_mc = mDh_diff_mc_ptr->ProjectionX("mDh_proj_diff_mc", 0, 1, "e");
  
  double evt_threshold = 300.0;
  int pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin;
  vector<double> fitresult;
  
  TH1D *mean_diff_by_idx = new TH1D("mean_diff_by_idx", "mean_diff_by_idx", nbins, 0, nbins);
  TH1D *sigma_diff_by_idx = new TH1D("sigma_smear_by_idx", "mean_diff_by_idx", nbins, 0, nbins);
  
  cout << "\n" << "First loop over 4D bins starts" << "\n";
  for (int multi_dim_bin=0; multi_dim_bin<nbins; multi_dim_bin++){
    
    pos_eta_bin = multi_dim_bin/(nbinspt*nbinseta*nbinspt);
    pos_pt_bin = (multi_dim_bin%(nbinspt*nbinseta*nbinspt))/(nbinseta*nbinspt);
    neg_eta_bin = ((multi_dim_bin%(nbinspt*nbinseta*nbinspt))%(nbinseta*nbinspt))/nbinspt;
    neg_pt_bin = ((multi_dim_bin%(nbinspt*nbinseta*nbinspt))%(nbinseta*nbinspt))%nbinspt;

    delete gROOT->FindObject("mDh_proj_diff_mc");
    mDh_proj_diff_mc = mDh_diff_mc_ptr->ProjectionX("mDh_proj_diff_mc", multi_dim_bin+1, multi_dim_bin+1, "e");
    if (mDh_proj_diff_mc->Integral(1,nbinsmll_diff) >= evt_threshold){

      // Fit diff_mc
      fitresult = fitHisto(mDh_proj_diff_mc, 0, 8, 5);
      if (fitresult[0] == -90.0 || fitresult[1] == -5.0){
        std::cout<<"WARNING bad diff_mc fit "<<" \n";
      }
      mean_diff_by_idx->SetBinContent(multi_dim_bin+1, fitresult[0]);
      mean_diff_by_idx->SetBinError(multi_dim_bin+1, fitresult[2]);
      sigma_diff_by_idx->SetBinContent(multi_dim_bin+1, fitresult[1]);
      sigma_diff_by_idx->SetBinError(multi_dim_bin+1, fitresult[3]);

    } else {
      mean_diff_by_idx->SetBinContent(multi_dim_bin+1, -10.0);
      mean_diff_by_idx->SetBinError(multi_dim_bin+1, 5.0);
      sigma_diff_by_idx->SetBinContent(multi_dim_bin+1, -10.0);
      sigma_diff_by_idx->SetBinError(multi_dim_bin+1, 5.0);
    }
  }

  // define jacobians with the info from that histo
  auto jacobian_weight_alpha_mass = [&](int multi_dim_bin, float m_reco, float m_gen, double MC_weight)->float{

    float mean_MC, sigma_MC;
    mean_MC = mean_diff_by_idx->GetBinContent(multi_dim_bin+1);
    sigma_MC = sigma_diff_by_idx->GetBinContent(multi_dim_bin+1);

    float jac_alpha = (m_reco-(mean_MC+m_gen))*(m_reco-(mean_MC+m_gen))/sigma_MC/sigma_MC-1.0;
    return jac_alpha*MC_weight;
  };

  auto jacobian_weight_beta_mass = [&](int multi_dim_bin, float m_reco, float m_gen, double MC_weight)->float{

    float mean_MC, sigma_MC;
    mean_MC = mean_diff_by_idx->GetBinContent(multi_dim_bin+1);
    sigma_MC = sigma_diff_by_idx->GetBinContent(multi_dim_bin+1);

    float jac_beta = (m_reco-(mean_MC+m_gen))*(mean_MC+m_gen)/sigma_MC/sigma_MC;
    return jac_beta*MC_weight;
  };

  auto jacobian_weight_alpha_diff = [&](int multi_dim_bin, float m_diff, double MC_weight)->float{

    float mean_MC, sigma_MC;
    mean_MC = mean_diff_by_idx->GetBinContent(multi_dim_bin+1);
    sigma_MC = sigma_diff_by_idx->GetBinContent(multi_dim_bin+1);

    float jac_alpha = (m_diff-mean_MC)*(m_diff-mean_MC)/sigma_MC/sigma_MC - 1.0;
    return jac_alpha*MC_weight;
  };

  auto jacobian_weight_epsilon_diff = [&](int multi_dim_bin, float m_diff, double MC_weight)->float{

    float mean_MC, sigma_MC;
    mean_MC = mean_diff_by_idx->GetBinContent(multi_dim_bin+1);
    sigma_MC = sigma_diff_by_idx->GetBinContent(multi_dim_bin+1);

    float jac_epsilon = (m_diff-mean_MC)/sigma_MC/sigma_MC;
    return jac_epsilon*MC_weight;
  };
  
  
  // Choose validation/analysis mode //TODO decide if we keep this in this script
  string mode_option("validation"), mc_name_root, data_name_root;

  // Input files
  if (mode_option.compare("analysis") == 0) {
    mc_name_root = "reco";
    data_name_root = "data2016";
  } else if (mode_option.compare("validation") == 0){
    mc_name_root = "smear";
    data_name_root = "smear_beta_val";
  }

  auto d2 = std::make_unique<RNode>(df);
  
  // Random numbers for pairs
  nslots = d2->GetNSlots();
  std::cout<<"3rd frame: nslots = "<<nslots<<"\n";
  for(unsigned int i = 0; i < nslots; i++){
    rans[i] = new TRandom3(17724 + i*10) ;
  }

  d2 = std::make_unique<RNode>(d2->Filter("HLT_IsoMu24 == 1"));
  d2 = std::make_unique<RNode>(d2->Filter("nMuon >= 2"));
  d2 = std::make_unique<RNode>(d2->Filter("PV_npvsGood >= 1"));

  d2 = std::make_unique<RNode>(d2->Define("dxy_significance","dxy_significance(Muon_dxy, Muon_dxyErr)"));
  d2 = std::make_unique<RNode>(d2->Define("MuonisGood", "MuonisGood(Muon_pt, Muon_eta, Muon_isGlobal, Muon_mediumId, Muon_pfRelIso04_all, Muon_genPartFlav, dxy_significance)"));

  d2 = std::make_unique<RNode>(d2->DefineSlot("pairs", pairs, {"Muon_pt", "Muon_charge", "Muon_eta", "Muon_phi", "MuonisGood", "Muon_dxy", "Muon_dz", "GenPart_status", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi"}));

  // This below works because actually we kept only one pair per event
  // Accessing properties through the idices of pairs ensures the muons passed MuonisGood

  // <0 pos_muon_index, 1 neg_muon_index, 2 mll_reco, 3 posPt_reco, 4 negPt_reco , 5 mll_smear, 6 posPt_smear, 7 negPt_smear, 8 mll_gen, 9 posPt_gen, 10 negPt_gen, 11 smear_beta_weight, 12 posPt_smear_beta_val, 13 negPt_smear_beta_val>
  
  d2 = std::make_unique<RNode>(d2->Define("mll_reco","return get<2>(pairs);"));
  d2 = std::make_unique<RNode>(d2->Filter("mll_reco>1.0")); // this means only events with one mu pair are kept
  d2 = std::make_unique<RNode>(d2->Define("mll_gen","return get<8>(pairs);"));
  d2 = std::make_unique<RNode>(d2->Define("mll_diff_reco","return mll_reco - mll_gen;"));
  d2 = std::make_unique<RNode>(d2->Define("posGenEta","return GenPart_eta[get<14>(pairs)];"));
  d2 = std::make_unique<RNode>(d2->Define("negGenEta","return GenPart_eta[get<15>(pairs)];"));
  d2 = std::make_unique<RNode>(d2->Define("weight", "std::copysign(1.0, genWeight)"));
  
  // define a frame that plays the role of data -> odd events, MC -> even events
  auto d_mc = std::make_unique<RNode>(df);
  auto d_sim_data = std::make_unique<RNode>(df);
  d_sim_data = std::make_unique<RNode>(d2->Filter("event%2==1"));
    
  d_sim_data = std::make_unique<RNode>(d_sim_data->Define("posPtSmearBetaVal","return get<12>(pairs);"));
  d_sim_data = std::make_unique<RNode>(d_sim_data->Define("negPtSmearBetaVal","return get<13>(pairs);"));
  d_sim_data = std::make_unique<RNode>(d_sim_data->Define("bin_index_smear_beta_val", GetEtaPtEtaPtBin ,{"posGenEta", "posPtSmearBetaVal", "negGenEta", "negPtSmearBetaVal"}));
  d_sim_data = std::make_unique<RNode>(d_sim_data->Define("mll_smear_beta_val","return get<11>(pairs);"));
  d_sim_data = std::make_unique<RNode>(d_sim_data->Define("mll_diff_smear_beta_val","return mll_smear_beta_val - mll_gen;"));

  // 2D histos of data mll_diff, mass vs bin_index_(smear_beta_val)
  // Call 3rd event loop
  std::cout<<"Call 3rd event loop: data histograms"<<"\n";

  // Dummy data mll_diff
  auto mll_diff_smear_beta_val_diffbin_4Dbin = d_sim_data->Histo2D({"mll_diff_smear_beta_val_diffbin_4Dbin", "mll_diff_smear_beta_val_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"mll_diff_smear_beta_val", "bin_index_smear_beta_val", "weight");
  // Dummy data mll
  auto mll_smear_beta_val_mllbin_4Dbin = d_sim_data->Histo2D({"mll_smear_beta_val_mllbin_4Dbin", "mll_smear_beta_val_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"mll_smear_beta_val", "bin_index_smear_beta_val", "weight");
  
  d_mc = std::make_unique<RNode>(d2->Filter("event%2==0"));

  d_mc = std::make_unique<RNode>(d_mc->Define("posPtSmear","return get<6>(pairs);"));
  d_mc = std::make_unique<RNode>(d_mc->Define("negPtSmear","return get<7>(pairs);"));
  d_mc = std::make_unique<RNode>(d_mc->Define("posGenPhi","return GenPart_phi[get<14>(pairs)];"));
  d_mc = std::make_unique<RNode>(d_mc->Define("negGenPhi","return GenPart_phi[get<15>(pairs)];"));

  // Call 4th event loop -> save MC TTree
  std::cout<<"Call 4th event loop: save MC TTree"<<"\n";
  // Take snapshot of MC frame with weight, posPtSmear, negPtSmear, posGenEta, negGenEta, mll_gen to correct the pT and recalculate mll_smear and jacs and 4D bin
  d_mc->Snapshot("Events", "InOutputFiles/snapshot_mc.root", {"posPtSmear", "negPtSmear", "posGenEta", "negGenEta", "posGenPhi", "negGenPhi", "mll_gen", "weight"});

  vector<double> A_from_iter = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  vector<double> e_from_iter = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  vector<double> M_from_iter = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Funct correct MC pT with A, e, M -> output of previous fit iteration
  auto CorrectPt = [&](float pos_pt, float neg_pt, float pos_eta, float neg_eta, float pos_phi, float neg_phi){
    TLorentzVector posTrack, negTrack, mother;
    float mll, rest_mass = 0.105658; // muMass = 0.105658 GeV
    vector<float> result(3);
    int pos_eta_bin = GetEtaBin(pos_eta);
    int neg_eta_bin = GetEtaBin(neg_eta);

    cout<<"A_from_iter: "<<A_from_iter<<"\n";
    cout<<"e_from_iter: "<<e_from_iter<<"\n";
    cout<<"M_from_iter: "<<M_from_iter<<"\n";
    
    float pos_pt_corr = pos_pt/( 1. + A_from_iter[pos_eta_bin] - e_from_iter[pos_eta_bin]/pos_pt + M_from_iter[pos_eta_bin]*pos_pt );
    float neg_pt_corr = neg_pt/( 1. + A_from_iter[neg_eta_bin] - e_from_iter[neg_eta_bin]/neg_pt - M_from_iter[neg_eta_bin]*neg_pt );

    result[0] = pos_pt_corr;
    result[1] = neg_pt_corr;
    posTrack.SetPtEtaPhiM(pos_pt_corr, pos_eta, pos_phi, rest_mass);
    negTrack.SetPtEtaPhiM(neg_pt_corr, neg_eta, neg_phi, rest_mass);
    mother = posTrack + negTrack;
    mll = mother.M();

    result[2] = mll;
    return result;
  };

  unique_ptr<TFile> f_control_const( TFile::Open("InOutputFiles/constants_fitted.root", "RECREATE") );
  
  // Start loop over iterations of correcting MC with: for iter < last_iter, Redefine() d_mc_corr frame
  std::cout<<"Start loop over iterations"<<"\n";

  int total_iter = 2;
  for(int iter=0; iter<total_iter; iter++){
    std::cout<<"Iter "<<iter<<"\n";
    
  // define corrected MC frame with the jacobian branches
  RDataFrame df_mc_corr("Events", "InOutputFiles/snapshot_mc.root");
  auto d_mc_corr = std::make_unique<RNode>(df_mc_corr);

  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("corrections", CorrectPt ,{"posPtSmear", "negPtSmear", "posGenEta", "negGenEta", "posGenPhi", "negGenPhi"}));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("posPtSmearCorr","return corrections[0];"));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("negPtSmearCorr","return corrections[1];"));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("bin_index_smear", GetEtaPtEtaPtBin ,{"posGenEta", "posPtSmearCorr", "negGenEta", "negPtSmearCorr"}));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("mll_smear","return corrections[2];"));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("mll_diff_smear","return mll_smear - mll_gen;"));
  
  // Smear jacobians
  // Additive k bin dependent PDF
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("jacobian_weight_mll_alpha_smear", jacobian_weight_alpha_mass, {"bin_index_smear", "mll_smear","mll_gen","weight"}));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("jacobian_weight_mll_beta_smear", jacobian_weight_beta_mass, {"bin_index_smear", "mll_smear","mll_gen","weight"}));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("jacobian_weight_mll_diff_alpha_smear", jacobian_weight_alpha_diff, {"bin_index_smear", "mll_diff_smear","weight"}));
  d_mc_corr = std::make_unique<RNode>(d_mc_corr->Define("jacobian_weight_mll_diff_epsilon_smear", jacobian_weight_epsilon_diff, {"bin_index_smear", "mll_diff_smear","weight"}));

  // ONE PDF for all k BIN
  // TODO define it too if you want
  //
  
  // Reco
  // TODO define it too
  //
  
  // 2D histos of mll_diff, mass, jacobian terms vs bin_index_(smear/smear_beta_val)
  // Call MC TTree event loop
  std::cout<<"Call MC TTree event loop: MC histograms"<<"\n";
  
  // binned in mll_diff
  // MC
  mll_diff_smear_diffbin_4Dbin = d_mc_corr->Histo2D({"mll_diff_smear_diffbin_4Dbin", "mll_diff_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries},"mll_diff_smear", "bin_index_smear", "weight");
  // Jac, mll_diff weighted by jac
  auto jac_weight_mll_diff_alpha_smear_diffbin_4Dbin = d_mc_corr->Histo2D({"jac_weight_mll_diff_alpha_smear_diffbin_4Dbin", "jac_weight_mll_diff_alpha_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries}, "mll_diff_smear", "bin_index_smear", "jacobian_weight_mll_diff_alpha_smear");
  auto jac_weight_mll_diff_epsilon_smear_diffbin_4Dbin = d_mc_corr->Histo2D({"jac_weight_mll_diff_epsilon_smear_diffbin_4Dbin", "jac_weight_mll_diff_epsilon_smear_diffbin_4Dbin", nbinsmll_diff, mymll_diffboundaries, nbins, mybinsboundaries}, "mll_diff_smear", "bin_index_smear", "jacobian_weight_mll_diff_epsilon_smear");
    
  // binned in mll
  // MC
  auto mll_smear_mllbin_4Dbin = d_mc_corr->Histo2D({"mll_smear_mllbin_4Dbin", "mll_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"mll_smear", "bin_index_smear", "weight");
  // Jac, mll weighted by jac
  auto jac_weight_mll_alpha_smear_mllbin_4Dbin = d_mc_corr->Histo2D({"jac_weight_mll_alpha_smear_mllbin_4Dbin", "jac_weight_mll_alpha_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"mll_smear","bin_index_smear","jacobian_weight_mll_alpha_smear");
  auto jac_weight_mll_beta_smear_mllbin_4Dbin = d_mc_corr->Histo2D({"jac_weight_mll_beta_smear_mllbin_4Dbin", "jac_weight_mll_beta_smear_mllbin_4Dbin", nbinsmll, mymllboundaries, nbins, mybinsboundaries},"mll_smear", "bin_index_smear", "jacobian_weight_mll_beta_smear");
    
  std::unique_ptr<TFile> fz( TFile::Open("InOutputFiles/test.root", "RECREATE") );
  mll_diff_smear_diffbin_4Dbin->Write();
  mll_smear_mllbin_4Dbin->Write();
  jac_weight_mll_alpha_smear_mllbin_4Dbin->Write();
  jac_weight_mll_beta_smear_mllbin_4Dbin->Write();
  fz->Close();

  /*
  //Control histograms
  auto mll_smear_hist = d_mc->Histo1D({"mll_smear", "mll inclusive all bins", 20, mll_low, mll_high},"mll_smear","weight");
  auto mll_diff_smear_hist = d_mc->Histo1D({"mll_diff_smear", "mll_diff inclusive all bins", 20, mll_diff_low, mll_diff_high},"mll_diff_smear","weight");
  auto mll_smear_beta_val_hist = d_sim_data->Histo1D({"mll_smear_beta_val", "mll inclusive all bins #beta != 1", 20, mll_low, mll_high},"mll_smear_beta_val","weight");
  auto pt_smear_hist = d_mc->Histo1D({"pt_smear", "pt smear beta = 1", 15, pt_low, pt_high},"posPtSmear","weight");
  auto pt_smear_beta_val_hist = d_sim_data->Histo1D({"pt_smear_beta_val", "pt smear #beta != 1", 15, pt_low, pt_high},"posPtSmearBetaVal","weight");

  std::unique_ptr<TFile> f2( TFile::Open("InOutputFiles/control_histo.root", "RECREATE") );
  mll_smear_hist->Write();
  mll_diff_smear_hist->Write();
  mll_smear_beta_val_hist->Write();
  pt_smear_hist->Write();
  pt_smear_beta_val_hist->Write();
  f2->Close();

  return 0;
  */
  
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Mass fitting part
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Variables declaration //TODO can move a lot of variable and hist declaration before loop over pT corrections

  //TODO comment explain where the variables are needed
  double total_nevents=0.0, nevents=0.0, all_histos_count=0.0, remaining_nevents=0.0, empty_histos_count=0.0, hfrac=-1.0, efrac=-1.0; // to check number of passing events
  double value_alpha, error_alpha, value_beta, error_beta, value_epsilon, error_epsilon; // for jac calculation
  double mean_mc = 0.0, error_mean_mc = 0.0, sigma_mc=0.0, error_sigma_mc = 0.0, integral_mc, chi2_mc;
  double mean_diff_data = 1.0, error_mean_diff_data = 1.0, sigma_diff_data = 1.0, error_sigma_diff_data = 1.0, integral_diff_data = 1.0; //TODO just initialise with 0, are they called in analysis mode at all?
  double fit_beta_error, value_corrected, error_corrected, mean_corrected_diff, error_mean_corrected_diff, sigma_corrected_diff, error_sigma_corrected_diff;
  double max_hist_mll_diff, max_hist_mll;
  int filled_bins_mll, filled_bins_mll_diff, position_to_fill;
  string leg_entry, name, title;

  // MC mll_diff already defined above
  
  auto mDh_mll_mc_ptr = mll_smear_mllbin_4Dbin;
  auto mDh_proj_mll_mc = mDh_mll_mc_ptr->ProjectionX("mDh_proj_mll_mc", 0, 1, "e");

  auto mDh_jac_alpha_mll_ptr = jac_weight_mll_alpha_smear_mllbin_4Dbin;
  auto mDh_proj_jac_alpha_mll = mDh_jac_alpha_mll_ptr->ProjectionX("mDh_proj_jac_alpha_mll", 0, 1, "e");

  auto mDh_jac_beta_mll_ptr = jac_weight_mll_beta_smear_mllbin_4Dbin;
  auto mDh_proj_jac_beta_mll = mDh_jac_beta_mll_ptr->ProjectionX("mDh_proj_jac_beta_mll", 0, 1, "e");

  auto mDh_mll_data_ptr = mll_smear_beta_val_mllbin_4Dbin;
  auto mDh_proj_mll_data = mDh_mll_data_ptr->ProjectionX("mDh_proj_mll_data", 0, 1, "e");

  auto mDh_jac_alpha_diff_ptr = jac_weight_mll_diff_alpha_smear_diffbin_4Dbin;
  auto mDh_proj_jac_alpha_diff = mDh_jac_alpha_diff_ptr->ProjectionX("mDh_proj_jac_alpha_diff", 0, 1, "e");

  auto mDh_jac_epsilon_diff_ptr = jac_weight_mll_diff_epsilon_smear_diffbin_4Dbin;
  auto mDh_proj_jac_epsilon_diff = mDh_jac_epsilon_diff_ptr->ProjectionX("mDh_proj_jac_epsilon_diff", 0, 1, "e");

  auto mDh_diff_data_ptr = mll_diff_smear_beta_val_diffbin_4Dbin;
  auto mDh_proj_diff_data = mDh_diff_data_ptr->ProjectionX("mDh_proj_diff_data", 0, 1, "e");

  // Save quantities relevant to control histograms in a TTree
  std::unique_ptr<TFile> f_control_tree( TFile::Open("InOutputFiles/control_tree.root", "RECREATE") );
  //auto control_tree = std::make_unique<TTree>("control_tree", "control_tree");
  TTree *control_tree = new TTree("control_tree", "control_tree");
  double fitted_beta = 0.0, fitted_beta_error = 0.0, fitted_alpha = 0.0, fitted_alpha_error = 0.0, fitted_nu = 0.0, fitted_nu_error = 0.0; 
  double* nevents_ptr = &nevents;
  double* fitted_beta_ptr = &fitted_beta;
  double* fitted_beta_error_ptr = &fitted_beta_error;
  double* fitted_alpha_ptr = &fitted_alpha;
  double* fitted_alpha_error_ptr = &fitted_alpha_error;
  double* fitted_nu_ptr = &fitted_nu;
  double* fitted_nu_error_ptr = &fitted_nu_error;
  double* mean_mc_ptr = &mean_mc;
  double* error_mean_mc_ptr = &error_mean_mc;
  double* sigma_mc_ptr = &sigma_mc;
  double* error_sigma_mc_ptr = &error_sigma_mc;
  string* name_ptr = &name;

  control_tree->Branch("nevents", nevents_ptr);
  control_tree->Branch("beta", fitted_beta_ptr);
  control_tree->Branch("error_beta", fitted_beta_error_ptr);
  control_tree->Branch("alpha", fitted_alpha_ptr);
  control_tree->Branch("error_alpha", fitted_alpha_error_ptr);
  control_tree->Branch("nu", fitted_nu_ptr);
  control_tree->Branch("error_nu", fitted_nu_error_ptr);
  control_tree->Branch("mean_mc", mean_mc_ptr);
  control_tree->Branch("error_mean_mc", error_mean_mc_ptr);
  control_tree->Branch("sigma_mc", sigma_mc_ptr);
  control_tree->Branch("error_sigma_mc", error_sigma_mc_ptr);
  control_tree->Branch("label", name_ptr);
    
  // Histograms for mll_diff distribution properties

  title = mc_name_root + "-gen mll mean";
  TH1D *mean_diff = new TH1D("mean_diff", title.c_str(), 3, 0, 3);
  mean_diff->SetCanExtend(TH1::kAllAxes);
  //mean_diff->SetMarkerStyle(kPlus);
  mean_diff->GetXaxis()->SetTitle("Bin number");
  title = mc_name_root + "-gen mll mean [GeV]";
  mean_diff->GetYaxis()->SetTitle(title.c_str());

  title = mc_name_root + "-gen mll sigma";
  TH1D *sigma_diff = new TH1D("sigma_smear", title.c_str(), 3, 0, 3);
  sigma_diff->SetCanExtend(TH1::kAllAxes);
  sigma_diff->GetXaxis()->SetTitle("Bin number");
  title = mc_name_root + "-gen mll sigma [GeV]";
  sigma_diff->GetYaxis()->SetTitle(title.c_str());
  
  TH1D *occupancy = new TH1D("bin_occupancy", "Events in mll_diff", 3, 0, 3);
  occupancy->SetCanExtend(TH1::kAllAxes);
  occupancy->SetMarkerStyle(kPlus);
  occupancy->SetMarkerColor(kBlue);
  occupancy->SetLineColor(kBlue);
  occupancy->GetXaxis()->SetTitle("Bin number");
  occupancy->GetYaxis()->SetTitle("Events");

  // Histograms for mll distribution properties

  TH1D *gaus_integral = new TH1D("gaus_integral", "mll integral 75.5-105.0", 3, 0, 3);
  gaus_integral->SetCanExtend(TH1::kAllAxes);
  gaus_integral->SetMarkerStyle(kPlus);
  gaus_integral->SetMarkerColor(kBlue);
  gaus_integral->SetLineColor(kBlue);
  gaus_integral->GetXaxis()->SetTitle("Bin number");
  gaus_integral->GetYaxis()->SetTitle("integral [evts]");

  // Histograms for jacobians inclusive in eta, pt
  
  TH1D *jac_alpha_inclusive = new TH1D("jacobian_alpha_inclusive", "jacobian alpha inclusive in eta, pt", nbinsmll, mll_low, mll_high);
  jac_alpha_inclusive->GetXaxis()->SetTitle("mll");
  jac_alpha_inclusive->GetYaxis()->SetTitle("jacobian alpha");

  TH1D *jac_beta_inclusive = new TH1D("jacobian_beta_inclusive", "jacobian beta inclusive in eta, pt", nbinsmll, mll_low, mll_high);
  jac_beta_inclusive->GetXaxis()->SetTitle("mll");
  jac_beta_inclusive->GetYaxis()->SetTitle("jacobian beta");

  // Histograms for pull of fitted variables per k bin

  TH1D *alpha = new TH1D("alpha", "#Delta#alpha", 3, 0, 3);
  alpha->SetCanExtend(TH1::kAllAxes);
  alpha->SetMarkerStyle(kPlus);
  alpha->SetMarkerColor(kBlue);
  alpha->GetXaxis()->SetTitle("Bin number");
  alpha->GetYaxis()->SetTitle("#alpha");

  TH1D *beta = new TH1D("beta", "#Delta#beta", 3, 0, 3);
  beta->SetCanExtend(TH1::kAllAxes);
  beta->SetMarkerStyle(kPlus);
  beta->SetMarkerColor(kBlue);
  beta->GetXaxis()->SetTitle("Bin number");
  beta->GetYaxis()->SetTitle("#Delta#beta");

  TH1D *nu = new TH1D("nu", "#Delta#nu", 3, 0, 3);
  nu->SetCanExtend(TH1::kAllAxes);
  nu->SetMarkerStyle(kPlus);
  nu->SetMarkerColor(kBlue);
  nu->GetXaxis()->SetTitle("Bin number");
  nu->GetYaxis()->SetTitle("#nu");

  TH1D *h_scales = new TH1D("h_scales", "scales", nbins, 0, nbins);
  h_scales->GetXaxis()->SetTitle("Bin idx");
  h_scales->GetYaxis()->SetTitle("scale");

  TH1D *h_masks = new TH1D("h_masks", "masks", nbins, 0, nbins);
  h_masks->GetXaxis()->SetTitle("Bin idx");
  h_masks->GetYaxis()->SetTitle("mask");

  TH1D *alpha_control=nullptr, *epsilon_control=nullptr, *nu_control=nullptr, *epsilon_test1=nullptr;
  TH2D *epsilon_test2=nullptr;
  
  alpha_control = new TH1D("alpha_control", "alpha from mll_diff", 3, 0, 3);
  alpha_control->SetCanExtend(TH1::kAllAxes);
  alpha_control->SetMarkerStyle(kPlus);
  alpha_control->SetMarkerColor(kBlue);
  alpha_control->GetXaxis()->SetTitle("Bin number");
  alpha_control->GetYaxis()->SetTitle("#alpha");

  epsilon_control = new TH1D("epsilon_control", "epsilon from mll_diff", 3, 0, 3);
  epsilon_control->SetCanExtend(TH1::kAllAxes);
  epsilon_control->SetMarkerStyle(kPlus);
  epsilon_control->SetMarkerColor(kBlack);
  epsilon_control->GetXaxis()->SetTitle("Bin number");
  epsilon_control->GetYaxis()->SetTitle("#epsilon");

  nu_control = new TH1D("nu_control", "nu from mll_diff", 3, 0, 3);
  nu_control->SetCanExtend(TH1::kAllAxes);
  nu_control->SetMarkerStyle(kPlus);
  nu_control->SetMarkerColor(kBlue);
  nu_control->GetXaxis()->SetTitle("Bin number");
  nu_control->GetYaxis()->SetTitle("#nu");
  
  epsilon_test1 = new TH1D("epsilon_test1", " (mean_yellow - mean_green) / sigma_green", 3, 0, 3);
  epsilon_test1->SetCanExtend(TH1::kAllAxes);
  epsilon_test1->SetMarkerStyle(kPlus);
  epsilon_test1->SetMarkerColor(kBlue);
  epsilon_test1->GetXaxis()->SetTitle("Bin number");
  epsilon_test1->GetYaxis()->SetTitle("Pull");
  
  epsilon_test2 = new TH2D("epsilon_test2", "epsilon_test2", 8, 0.1, 0.4, 8, -3.0, 3.0);
  epsilon_test2->GetXaxis()->SetTitle("abs((mean_yellow - mean_green) / sigma_green)");
  epsilon_test2->GetYaxis()->SetTitle("((mean_yellow - mean_green) - epsilon) / error_epsilon");

  // Histograms for pull distributions

  double limit_p = 3.0, limit_m = -3.0;
  int bins = (limit_p - limit_m)/0.1;

  TH1D *pull_alpha_control=nullptr, *pull_nu_control=nullptr, *pull_epsilon_control=nullptr, *pull_epsilon_control1=nullptr, *pull_epsilon_control2=nullptr;
  
  pull_alpha_control = new TH1D("pull_alpha_control", " ((sigma_yellow - sigma_green) - sigma_green*Delta_alpha) / error_alpha ", bins, limit_m, limit_p);
  pull_alpha_control->GetXaxis()->SetTitle("Pull");
  pull_alpha_control->GetYaxis()->SetTitle("Events");
  
  pull_nu_control = new TH1D("pull_nu_control", " ((integral_yellow - integral_green) - integral_green*Delta_nu) / error_nu ", bins, limit_m, limit_p);
  pull_nu_control->GetXaxis()->SetTitle("Pull");
  pull_nu_control->GetYaxis()->SetTitle("Events");
  
  pull_epsilon_control = new TH1D("pull_epsilon_control", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze none) ", bins, limit_m, limit_p);
  pull_epsilon_control->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control->GetYaxis()->SetTitle("Events");
  
  pull_epsilon_control1 = new TH1D("pull_epsilon_control1", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze nu) ", bins, limit_m, limit_p);
  pull_epsilon_control1->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control1->GetYaxis()->SetTitle("Events");
  
  pull_epsilon_control2 = new TH1D("pull_epsilon_control2", " ((mean_yellow - mean_green) - epsilon) / error_epsilon (freeze nu, alpha) ", bins, limit_m, limit_p);
  pull_epsilon_control2->GetXaxis()->SetTitle("Pull");
  pull_epsilon_control2->GetYaxis()->SetTitle("Events");

  //--------------------------------------------------------------------------------------
  // Files to write results
  title = "InOutputFiles/mass_fits_control_histos_" + data_name_root + ".root";
  std::unique_ptr<TFile> f_control( TFile::Open(title.c_str(), "RECREATE") ); // histos inclusive in k bins
  title = "InOutputFiles/mass_fits_" + data_name_root + ".root";
  std::unique_ptr<TFile> f_fits( TFile::Open(title.c_str(), "RECREATE") ); // histos per k bin
  ofstream f_pass_reg("InOutputFiles/passed_regions.txt"); // to check if there are enough k bins to constrain all sagitta parameters TODO get from python
  
  // Canvas
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  c1->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  auto leg2 = new TLegend(0.58, 0.68, 0.90, 0.90);

  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.015);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.015);
  leg2->SetFillColor(10);
  leg2->SetNColumns(1);
  leg2->SetHeader("");

  //--------------------------------------------------------------------------------------
  // Prepare counting events in mll_diff

  total_nevents = 0.0;
  all_histos_count = 0.0;
  empty_histos_count = 0.0;
  hfrac = -1.0;  
  efrac = -1.0;
  remaining_nevents = 0.0;

  //--------------------------------------------------------------------------------------  
  // Loop over 4D bins
  
  cout << "\n" << "Loop over 4D bins starts" << "\n";
  for (int multi_dim_bin=0; multi_dim_bin<nbins; multi_dim_bin++){ 

    all_histos_count++;

    pos_eta_bin = multi_dim_bin/(nbinspt*nbinseta*nbinspt);
    pos_pt_bin = (multi_dim_bin%(nbinspt*nbinseta*nbinspt))/(nbinseta*nbinspt);
    neg_eta_bin = ((multi_dim_bin%(nbinspt*nbinseta*nbinspt))%(nbinseta*nbinspt))/nbinspt;
    neg_pt_bin = ((multi_dim_bin%(nbinspt*nbinseta*nbinspt))%(nbinseta*nbinspt))%nbinspt;

    h_scales->SetBinContent(multi_dim_bin, 0.0);
    h_scales->SetBinError(multi_dim_bin, 10.0);
    
    h_masks->SetBinContent(multi_dim_bin, 0);
    
    // Require enough stats in diff_mc for sigma_MC fit
    delete gROOT->FindObject("mDh_proj_diff_mc");
    mDh_proj_diff_mc = mDh_diff_mc_ptr->ProjectionX("mDh_proj_diff_mc", multi_dim_bin+1, multi_dim_bin+1, "e");
    nevents = mDh_proj_diff_mc->Integral(1,nbinsmll_diff); 
    total_nevents += nevents;
    
    if (nevents < evt_threshold){ // reject low stats 
      empty_histos_count++;

    } else {

      // Require most of mll_mc Gaussian to be in fit window
      // TODO and require same from mll_data
      delete gROOT->FindObject( "mDh_proj_mll_mc" );
      mDh_proj_mll_mc = mDh_mll_mc_ptr->ProjectionX("mDh_proj_mll_mc", multi_dim_bin+1, multi_dim_bin+1, "e");
      mDh_proj_mll_mc->GetXaxis()->SetTitle("mll [GeV]");
      mDh_proj_mll_mc->GetYaxis()->SetTitle("Events");
      fitresult = fitHisto(mDh_proj_mll_mc, 0, 8, 5);
      
      if (fitresult[4] < 0.75){ // reject small gaus integral 
	empty_histos_count++;

      } else {

	name = stringify_name(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin);
	//std::cout<< name <<"\n";
	//f_pass_reg << name << "\n"; do this after you check there are at least 3 points in mass jacobian
	                  
	remaining_nevents += nevents;      
	occupancy->SetBinError(occupancy->Fill(name.c_str(), nevents), 100);
	gaus_integral->SetBinError(gaus_integral->Fill(name.c_str(), fitresult[4]), 0.01);

	//--------------------------------------------------------------------------
	// Diff histograms
	//--------------------------------------------------------------------------
	
	//--------------------------------------------------------------------------
	// diff_mc histogram
	
	// Already projected diff_mc
	mDh_proj_diff_mc->GetXaxis()->SetTitle("mll_diff [GeV]");
	mDh_proj_diff_mc->GetYaxis()->SetTitle("Events");
	mDh_proj_diff_mc->SetLineColor(kGreen);
	// Fit diff_mc
	fitresult = fitHisto(mDh_proj_diff_mc, 1, 8, 5);
	if (fitresult[0] == -90.0 || fitresult[1] == -5.0){
	  std::cout<<"WARNING bad diff_mc fit "<< name.c_str() <<" \n";
	} 
	mean_diff->SetBinError(mean_diff->Fill(name.c_str(), fitresult[0]), fitresult[2]);
	sigma_diff->SetBinError(sigma_diff->Fill(name.c_str(), fitresult[1]), fitresult[3]);
	
	// Save for nu, beta/epsilon, alpha fit //TODO specify is for diff?
	mean_mc = fitresult[0];
	if (mean_diff_by_idx->GetBinContent(multi_dim_bin+1) != mean_mc) cout << " mean_diff_by_idx->GetBinContent(multi_dim_bin+1) != mean_mc "<< name.c_str() << mean_diff_by_idx->GetBinContent(multi_dim_bin+1) <<" "<< mean_mc <<"\n";
	error_mean_mc = fitresult[2];
	sigma_mc = fitresult[1];
	error_sigma_mc = fitresult[3];
	chi2_mc = fitresult[5];
	
	//-----------------------------------------------
	// Prepare to draw mll_diff panel
	
	c1->cd(1); 
	max_hist_mll_diff = -1.0;      
	
	if(max_hist_mll_diff < mDh_proj_diff_mc->GetBinContent(mDh_proj_diff_mc->GetMaximumBin())){
	  max_hist_mll_diff = mDh_proj_diff_mc->GetBinContent(mDh_proj_diff_mc->GetMaximumBin());
	}
	
	//--------------------------------------------------------------------------
	// diff_data histogram 
	
	if (mode_option.compare("validation") == 0){
	  // Project diff_data histogram
	  delete gROOT->FindObject("mDh_proj_diff_data");
	  mDh_proj_diff_data = mDh_diff_data_ptr->ProjectionX("mDh_proj_diff_data", multi_dim_bin+1, multi_dim_bin+1, "e");
	  mDh_proj_diff_data->GetXaxis()->SetTitle("mll_diff [GeV]");
	  mDh_proj_diff_data->GetYaxis()->SetTitle("Events");
	  mDh_proj_diff_data->SetLineColor(kOrange+7);
	  // Fit diff_data
	  fitresult = fitHisto(mDh_proj_diff_data, 1, 95, 5);
		      
	  // Save these for pull histograms 
	  mean_diff_data = fitresult[0]; 
	  error_mean_diff_data = fitresult[2];
	  sigma_diff_data = fitresult[1];
	  error_sigma_diff_data = fitresult[3];
		
	  if(max_hist_mll_diff < mDh_proj_diff_data->GetBinContent(mDh_proj_diff_data->GetMaximumBin())){
	    max_hist_mll_diff = mDh_proj_diff_data->GetBinContent(mDh_proj_diff_data->GetMaximumBin());
	  }
	} 

	//--------------------------------------------------------------------------
	// Start drawing mll_diff 
	                  
	max_hist_mll_diff = max_hist_mll_diff * 1.3;
	mDh_proj_diff_mc->SetMaximum(max_hist_mll_diff);
	mDh_proj_diff_mc->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	mDh_proj_diff_mc->Draw("HIST");
	if (mode_option.compare("validation") == 0) mDh_proj_diff_data->Draw("HIST SAME");
	
	// Legend
	leg1->Clear();
	leg_entry = "Region " + name + ": " + nevents + " events";
	leg1->SetHeader(leg_entry.c_str(),"C");
                    
	leg1->AddEntry(mDh_proj_diff_mc, "smeared-gen, #beta_pT=1, #alpha=1", "l");
	if (mode_option.compare("validation") == 0) leg1->AddEntry(mDh_proj_diff_data, "smeared-gen, #beta_pT!=1, #alpha=1", "l");
	leg_entry = "green fit: #mu=" + to_string(mean_mc).substr(0, 5) + ", #sigma=" + to_string(sigma_mc).substr(0, 5);
	leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	
	//--------------------------------------------------------------------------
	// Mll histograms
	//--------------------------------------------------------------------------
	
	//--------------------------------------------------------------------------
	// Prepare to draw mll panel
	
	c1->cd();
	c1->cd(2); 
	max_hist_mll = -1.0;
	
	//--------------------------------------------------------------------------
	// mll mc histogram

	// Already projected and fitted mll mc
	mDh_proj_mll_mc->SetTitle(("mll "+ stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	mDh_proj_mll_mc->GetXaxis()->SetTitle("mll [GeV]");
	mDh_proj_mll_mc->GetYaxis()->SetTitle("Events");
	mDh_proj_mll_mc->SetLineColor(kGreen);
	                  
	if(max_hist_mll < mDh_proj_mll_mc->GetBinContent(mDh_proj_mll_mc->GetMaximumBin())){
	  max_hist_mll = mDh_proj_mll_mc->GetBinContent(mDh_proj_mll_mc->GetMaximumBin());
	}

	//--------------------------------------------------------------------------
	// mll data
	                  
	// Project mll data
	delete gROOT->FindObject("mDh_proj_mll_data");
	mDh_proj_mll_data = mDh_mll_data_ptr->ProjectionX("mDh_proj_mll_data", multi_dim_bin+1, multi_dim_bin+1, "e");
	mDh_proj_mll_data->GetXaxis()->SetTitle("mll [GeV]");
	mDh_proj_mll_data->GetYaxis()->SetTitle("Events");
	mDh_proj_mll_data->SetLineColor(kOrange+7);
	// Fit mll data
	//fitresult = fitHisto(mDh_proj_mll_data, 0, 5, 5);
              
	if(max_hist_mll < mDh_proj_mll_data->GetBinContent(mDh_proj_mll_data->GetMaximumBin())){
	  max_hist_mll = mDh_proj_mll_data->GetBinContent(mDh_proj_mll_data->GetMaximumBin());
	}

	//--------------------------------------------------------------------------
	// Start drawing mll
	      
	max_hist_mll = max_hist_mll * 1.3;

	mDh_proj_mll_mc->SetMaximum(max_hist_mll);
	mDh_proj_mll_mc->Draw("HIST");
	mDh_proj_mll_data->Draw("HIST SAME");

	// Legend 
	leg2->Clear();
	leg_entry = "Region " + name + ": " + nevents + " events";
	leg2->SetHeader(leg_entry.c_str(),"C"); 
	
	leg2->AddEntry(mDh_proj_mll_mc, "smeared, #beta=1", "l"); //TODO automatise labels
	leg2->AddEntry(mDh_proj_mll_data, "smeared, #beta!=1", "l"); 
        
	c1->cd();
	
	//--------------------------------------------------------------------------
	// Start mass fit
	//--------------------------------------------------------------------------
	
	//--------------------------------------------------------------------------
	// Fill vectors and variance for minimisation
	
	// Find mll bins with defined jacobian
	
	filled_bins_mll=0;
	vector<int> good_indices_mll;
	for(int i=1; i<=nbinsmll; i++){
	  //asking jacobian to be average over at least 20 events
	  if (mDh_proj_mll_mc->GetBinContent(i) >= 20 && mDh_proj_mll_data->GetBinContent(i) > 0){ //TODO refine the other criterium
	    good_indices_mll.push_back(i);
	    filled_bins_mll++;
	  }
	}
	if (filled_bins_mll < 5 ){
	  //std::cout<< "WARNING not enough points in mll to fit nu, beta, alpha in "<< name.c_str() << "\n";
	  remaining_nevents -= nevents;
	  empty_histos_count++;
	  continue;
	} 
	f_pass_reg << name << "\n";
	      
	// Declare mass vectors, jacobians and variance
	      
	VectorXd h_data_minus_mc_mll_vector(filled_bins_mll);
	Eigen::MatrixXd V_inv_sqrt(filled_bins_mll, filled_bins_mll), J(filled_bins_mll, 3); //J.col(0) is nu, (1) is beta, (2) is alpha

	V_inv_sqrt = MatrixXd::Zero(filled_bins_mll, filled_bins_mll);
	position_to_fill=0;
	for(int i=1; i<=nbinsmll; i++){
	  if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
	    h_data_minus_mc_mll_vector(position_to_fill) = mDh_proj_mll_data->GetBinContent(i) - mDh_proj_mll_mc->GetBinContent(i);
	    J(position_to_fill,0) = mDh_proj_mll_mc->GetBinContent(i); //J.col(0) is nu
	    V_inv_sqrt(position_to_fill,position_to_fill) = 1.0 / mDh_proj_mll_data->GetBinError(i); // 1/sqrt(data_stat**2) 
	    position_to_fill++;
	  }
	}
	if (position_to_fill != filled_bins_mll){ std::cout<<"problem counting vector size \n"; }
	
	//--------------------------------------------------------------------------
	// Define jacobian histograms
	      
	// Jacobian histogram alpha mass
	delete gROOT->FindObject("jacobian_alpha");
	TH1D *jac_alpha = new TH1D("jacobian_alpha", "jacobian alpha", nbinsmll, mll_low, mll_high);
	jac_alpha->GetXaxis()->SetTitle("mll [GeV]");
	jac_alpha->GetYaxis()->SetTitle("jacobian [GeV]");
	      
	// Jacobian histogram beta mass
	delete gROOT->FindObject("jacobian_beta");
	TH1D *jac_beta = new TH1D("jacobian_beta", "jacobian beta", nbinsmll, mll_low, mll_high);
	jac_beta->GetXaxis()->SetTitle("mll [GeV]");
	jac_beta->GetYaxis()->SetTitle("jacobian [GeV]");
	
	//--------------------------------------------------------------------------
	// Compute jacobians mass
		
	delete gROOT->FindObject("mDh_proj_jac_alpha_mll");
	mDh_proj_jac_alpha_mll = mDh_jac_alpha_mll_ptr->ProjectionX("mDh_proj_jac_alpha_mll", multi_dim_bin+1, multi_dim_bin+1, "e");
	
	delete gROOT->FindObject("mDh_proj_jac_beta_mll");
        mDh_proj_jac_beta_mll = mDh_jac_beta_mll_ptr->ProjectionX("mDh_proj_jac_beta_mll", multi_dim_bin+1, multi_dim_bin+1, "e");
	
	position_to_fill=0;
	for(int i=1; i<=nbinsmll; i++){
	  if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
	    
	    // alpha mass
	    value_alpha = mDh_proj_jac_alpha_mll->GetBinContent(i);
	    error_alpha = mDh_proj_jac_alpha_mll->GetBinError(i);
	    J(position_to_fill,2) = value_alpha; //J.col(2) is alpha
	    
	    jac_alpha->SetBinContent(i, value_alpha);
	    jac_alpha->SetBinError(i, error_alpha);
	    jac_alpha_inclusive->Fill(mllbinranges[i-1], value_alpha);
	    
	    // beta mass
	    value_beta = mDh_proj_jac_beta_mll->GetBinContent(i);
	    error_beta = mDh_proj_jac_beta_mll->GetBinError(i);
	    J(position_to_fill,1) = value_beta; //J.col(1) is beta
	    
	    jac_beta->SetBinContent(i, value_beta);
	    jac_beta->SetBinError(i, error_beta);
	    jac_beta_inclusive->Fill(mllbinranges[i-1], value_beta);
	    
	    position_to_fill++;
	    
	  }
	}
				
	if (position_to_fill != filled_bins_mll){ std::cout<<"PROBLEM: counting jac mass size \n"; }

	//--------------------------------------------------------------------------
	// Write jacobian histograms
	
	f_fits->WriteObject(jac_alpha, ("jac_alpha" + name).c_str());
	f_fits->WriteObject(jac_beta, ("jac_beta" + name).c_str());
	
	//--------------------------------------------------------------------------
	// Solve for nu, beta, alpha mass
	//--------------------------------------------------------------------------
	
	// Solve for nu, beta, alpha mass
	Eigen::MatrixXd A = V_inv_sqrt*J;
	Eigen::MatrixXd b = V_inv_sqrt*h_data_minus_mc_mll_vector;
	// ATTENTION: n_b_a_vector contains nu, beta, alpha
	Eigen::VectorXd n_b_a_vector = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	
	// Error on nu, beta, alpha mass
	Eigen::MatrixXd V_nu = (J.col(0).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	Eigen::MatrixXd V_beta = (J.col(1).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	fit_beta_error = pow(V_beta(0,0),0.5); // save for closure test
	Eigen::MatrixXd V_alpha = (J.col(2).transpose()*(V_inv_sqrt*V_inv_sqrt)*J.col(2)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	
	// Write nu, beta, alpha mass
	fitted_nu = n_b_a_vector(0);
	fitted_nu_error = pow(V_nu(0,0),0.5);
	fitted_beta = n_b_a_vector(1);
	fitted_beta_error = pow(V_beta(0,0),0.5);
	fitted_alpha = n_b_a_vector(2);
	fitted_alpha_error = pow(V_alpha(0,0),0.5);
	
	nu->SetBinError(nu->Fill(name.c_str(), fitted_nu), fitted_nu_error);
	beta->SetBinError(beta->Fill(name.c_str(), fitted_beta), fitted_beta_error);
	alpha->SetBinError(alpha->Fill(name.c_str(), fitted_alpha), fitted_alpha_error);

	h_scales->SetBinContent(multi_dim_bin, 1.+ fitted_beta);
	h_scales->SetBinError(multi_dim_bin, fitted_beta_error);
	h_masks->SetBinContent(multi_dim_bin, 1);
	
	//--------------------------------------------------------------------------
	// Add fitted parameters to mass fit panel
	
	c1->cd(2); //TODO  #Delta?#beta=
	leg_entry = "Fit #Delta#beta=" + to_string(fitted_beta).substr(0, 6) + "#pm" + to_string(fitted_beta_error).substr(0, 6);
	leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	leg_entry = "Fit #Delta#nu=" + to_string(fitted_nu).substr(0, 6) + "#pm" + to_string(fitted_nu_error).substr(0, 6);
	leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	leg_entry = "Fit #Delta#alpha=" + to_string(fitted_alpha).substr(0, 6) + "#pm" + to_string(fitted_alpha_error).substr(0, 6);
	leg2->AddEntry((TObject*)0, leg_entry.c_str(), "");
	
	//--------------------------------------------------------------------------
	// Apply fitted nu, beta, alpha correction mass
        
	// Define corrected_mll histogram
	delete gROOT->FindObject("corrected_mll");
	TH1D *corrected_mll = new TH1D("corrected_mll", "corrected_mll", nbinsmll, mll_low, mll_high);
	corrected_mll->SetTitle(("mll " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	corrected_mll->GetXaxis()->SetTitle("mll [GeV]");
	corrected_mll->GetYaxis()->SetTitle("Events");
	
	position_to_fill = 0;
	for(int i=1; i<=nbinsmll; i++){
	  if ( find(good_indices_mll.begin(), good_indices_mll.end(), i) != good_indices_mll.end() ){
	    value_corrected = mDh_proj_mll_mc->GetBinContent(i) + fitted_beta*J(position_to_fill,1) + fitted_nu*J(position_to_fill,0) + fitted_alpha*J(position_to_fill,2);
	    error_corrected = pow(
				  pow(mDh_proj_mll_mc->GetBinError(i),2) +
				  pow(fitted_beta*jac_beta->GetBinError(i),2) +
				  pow(J(position_to_fill,1),2)*V_beta(0,0) +
				  pow(fitted_nu*mDh_proj_mll_mc->GetBinError(i),2) +
				  pow(J(position_to_fill,0),2)*V_nu(0,0) +
				  pow(fitted_alpha*jac_alpha->GetBinError(i),2) +
				  pow(J(position_to_fill,2),2)*V_alpha(0,0)
			      ,0.5);

	    corrected_mll->SetBinContent(i, value_corrected);
	    corrected_mll->SetBinError(i, error_corrected);

	    position_to_fill++;
	  }
	}
	if (position_to_fill != filled_bins_mll){ std::cout<<"PROBLEM: counting corrected mll bins \n"; }
	      
	//--------------------------------------------------------------------------
	// Draw mll
	      
	c1->cd();
	c1->cd(2);
	corrected_mll->Draw("HIST SAME");
	      
	// Legend
	leg2->AddEntry(corrected_mll, "corrected mll", "l");
	leg2->Draw("");
	      
	if (mode_option.compare("validation") == 0){
	  
	  //--------------------------------------------------------------------------
	  // Start diff fit
	  //--------------------------------------------------------------------------
	  
	  //--------------------------------------------------------------------------
	  // Find mll_diff bins with defined jacobian
		
	  filled_bins_mll_diff=0;
	  vector<int> good_indices_mll_diff;
	  for(int i=1; i<=nbinsmll_diff; i++){
	    //asking jacobian to be average over at least 20 events
	    if( mDh_proj_diff_mc->GetBinContent(i) >= 20 && mDh_proj_diff_data->GetBinContent(i) > 0 ){ // TODO refine the other criterium
	      good_indices_mll_diff.push_back(i);
	      filled_bins_mll_diff++;
	    }
	  }
	  if (filled_bins_mll_diff < 3 ){ //TODO refine, maybe I want to carry on with mass fit
	    std::cout<<"WARNING not enough points in mll_diff to fit nu, epsilon, alpha in "<< name.c_str() << "\n";
	    remaining_nevents -= nevents;
	    empty_histos_count++;
	    continue;
	  } 
		
	  // Declare mll_diff vectors, jacobians and variance
	  VectorXd h_data_minus_mc_diff_vector(filled_bins_mll_diff);
	  Eigen::MatrixXd V_inv_sqrt_control(filled_bins_mll_diff, filled_bins_mll_diff), J_control(filled_bins_mll_diff, 3); //J_control.col(0) is nu, (1) is epsilon, (2) is alpha 
		
	  V_inv_sqrt_control = MatrixXd::Zero(filled_bins_mll_diff, filled_bins_mll_diff);
	  position_to_fill = 0;
	  integral_diff_data = 0.0;
	  integral_mc = 0.0;
	  
	  for(int i=1; i<=nbinsmll_diff; i++){
	    if( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
	      h_data_minus_mc_diff_vector(position_to_fill) = mDh_proj_diff_data->GetBinContent(i) - mDh_proj_diff_mc->GetBinContent(i);
	      J_control(position_to_fill,0) = mDh_proj_diff_mc->GetBinContent(i); //J_control.col(0) is nu
	      V_inv_sqrt_control(position_to_fill,position_to_fill) = 1.0 / mDh_proj_diff_data->GetBinError(i); // 1/sqrt(data_stat**2)
	      integral_diff_data += mDh_proj_diff_data->GetBinContent(i);
	      integral_mc += mDh_proj_diff_mc->GetBinContent(i);
	      
	      position_to_fill++;
	    }
	  }
	  if (position_to_fill != filled_bins_mll_diff){ std::cout<<"problem counting vector size \n"; }
	  
	  //--------------------------------------------------------------------------
	  // Define jacobian histograms
		
	  // Jacobian histogram alpha diff 
	  delete gROOT->FindObject("jacobian_alpha_control");
	  TH1D *jac_alpha_control = new TH1D("jacobian_alpha_control", "jacobian alpha", nbinsmll_diff, mll_diff_low, mll_diff_high);
	  jac_alpha_control->GetXaxis()->SetTitle("mll_diff [GeV]");
	  jac_alpha_control->GetYaxis()->SetTitle("jacobian [GeV]");
	  
	  // Jacobian histogram epsilon diff
	  delete gROOT->FindObject("jacobian_epsilon_control");
	  TH1D *jac_epsilon_control = new TH1D("jacobian_epsilon_control", "jacobian epsilon", nbinsmll_diff, mll_diff_low, mll_diff_high);
	  jac_epsilon_control->GetXaxis()->SetTitle("mll_diff [GeV]");
	  jac_epsilon_control->GetYaxis()->SetTitle("jacobian [GeV]");
	  
	  //--------------------------------------------------------------------------
	  // Compute jacobians diff 
		
	  delete gROOT->FindObject("mDh_proj_jac_alpha_diff");
          mDh_proj_jac_alpha_diff = mDh_jac_alpha_diff_ptr->ProjectionX("mDh_proj_jac_alpha_diff", multi_dim_bin+1, multi_dim_bin+1, "e");
	  
	  delete gROOT->FindObject("mDh_proj_jac_epsilon_diff");
	  mDh_proj_jac_epsilon_diff = mDh_jac_epsilon_diff_ptr->ProjectionX("mDh_proj_jac_epsilon_diff", multi_dim_bin+1, multi_dim_bin+1, "e");
	  	  
	  position_to_fill = 0;
	  for(int i=1; i<=nbinsmll_diff; i++){
	    if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
	      // alpha diff
	      value_alpha = mDh_proj_jac_alpha_diff->GetBinContent(i);
	      error_alpha = mDh_proj_jac_alpha_diff->GetBinError(i);
	      
	      J_control(position_to_fill,2) = value_alpha; //J_control.col(2) is for alpha
	      jac_alpha_control->SetBinContent(i, value_alpha);
	      jac_alpha_control->SetBinError(i, error_alpha);
		    
	      // epsilon diff
	      value_epsilon = mDh_proj_jac_epsilon_diff->GetBinContent(i);
	      error_epsilon = mDh_proj_jac_epsilon_diff->GetBinError(i);
		    
	      J_control(position_to_fill,1) = value_epsilon; //J_control.col(1) is for epsilon
	      jac_epsilon_control->SetBinContent(i, value_epsilon);
	      jac_epsilon_control->SetBinError(i, error_epsilon);
	      
	      position_to_fill++;
	      
	    } 
	  }
	  
	  if (position_to_fill != filled_bins_mll_diff){ std::cout<<"PROBLEM: counting jac diff size \n"; }
		
	  //--------------------------------------------------------------------------
	  // Write jacobian histograms
		
	  f_fits->WriteObject(jac_alpha_control, ("jac_alpha_control" + name).c_str());
	  f_fits->WriteObject(jac_epsilon_control, ("jac_epsilon_control" + name).c_str());
	  	  
	  //--------------------------------------------------------------------------
		
	  // Solve for nu, epsilon, alpha diff 
	  Eigen::MatrixXd A_control = V_inv_sqrt_control*J_control;
	  Eigen::MatrixXd b_control = V_inv_sqrt_control*h_data_minus_mc_diff_vector;
	  // ATTENTION: n_e_a_vector_control contains nu, epsilon, alpha in this order
	  Eigen::VectorXd n_e_a_vector_control = A_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_control);
	  
	  // Solve for epsilon, alpha diff only
	  Eigen::MatrixXd A1_control = V_inv_sqrt_control*J_control.rightCols(2); 
	  // ATTENTION: e_a_vector_control contains epsilon, alpha
	  Eigen::VectorXd e_a_vector_control = A1_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_control);
	  
	  // Solve for epsilon diff only
	  Eigen::MatrixXd A2_control = V_inv_sqrt_control*J_control.col(1);
	  // ATTENTION: e_vector_control contains epsilon
	  Eigen::VectorXd e_vector_control = A2_control.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_control);
	  
	  // Error on epsilon diff 
	  Eigen::MatrixXd V_epsilon_control(1,1); 
	  V_epsilon_control = (J_control.col(1).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(1)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	  epsilon_control->SetBinError(epsilon_control->Fill(name.c_str(), n_e_a_vector_control(1)), pow(V_epsilon_control(0,0),0.5));
	  // Error on nu diff
	  Eigen::MatrixXd V_nu_control(1,1); 
	  V_nu_control = (J_control.col(0).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(0)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	  nu_control->SetBinError(nu_control->Fill(name.c_str(), n_e_a_vector_control(0)), pow(V_nu_control(0,0),0.5)); 
	  // Error on alpha diff 
	  Eigen::MatrixXd V_alpha_control(1,1);
	  V_alpha_control = (J_control.col(2).transpose()*(V_inv_sqrt_control*V_inv_sqrt_control)*J_control.col(2)).completeOrthogonalDecomposition().solve(MatrixXd::Identity(1,1));
	  alpha_control->SetBinError(alpha_control->Fill(name.c_str(), n_e_a_vector_control(2)), pow(V_alpha_control(0,0),0.5));
	  
	  //--------------------------------------------------------------------------
	  // Add fitted parameters to diff fit panel 
	  
	  c1->cd();
	  c1->cd(1);
	  
	  leg_entry = "Fit #varepsilon=" + to_string(n_e_a_vector_control(1)).substr(0, 6) + "#pm" + to_string(pow(V_epsilon_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(1) is epsilon
	  leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	  leg_entry = "Fit #Delta#nu=" + to_string(n_e_a_vector_control(0)).substr(0, 6) + "#pm" + to_string(pow(V_nu_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(0) is nu 
	  leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	  leg_entry = "Fit #Delta#alpha=" + to_string(n_e_a_vector_control(2)).substr(0, 6) + "#pm" + to_string(pow(V_alpha_control(0,0),0.5)).substr(0, 6); // n_e_a_vector_control(2) is alpha
	  leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	  
	  //--------------------------------------------------------------------------
	  // Define corrected_diff histogram
	  delete gROOT->FindObject("corrected_diff");
	  TH1D *corrected_diff = new TH1D("corrected_diff", "corrected_diff", nbinsmll_diff, mll_diff_low, mll_diff_high);
	  corrected_diff->SetTitle(("mll diff " + stringify_title(pos_eta_bin, pos_pt_bin, neg_eta_bin, neg_pt_bin, etabinranges, ptbinranges)).c_str());
	  corrected_diff->GetXaxis()->SetTitle("mll_diff [GeV]");
	  corrected_diff->GetYaxis()->SetTitle("Events");
		
	  position_to_fill = 0;
	  for(int i=1; i<=nbinsmll_diff; i++){
	    if ( find(good_indices_mll_diff.begin(), good_indices_mll_diff.end(), i) != good_indices_mll_diff.end() ){
	      value_corrected = mDh_proj_diff_mc->GetBinContent(i) + n_e_a_vector_control(1)*J_control(position_to_fill,1) + n_e_a_vector_control(0)*J_control(position_to_fill,0) + n_e_a_vector_control(2)*J_control(position_to_fill,2);
	      error_corrected = pow(
				    pow(mDh_proj_diff_mc->GetBinError(i),2) +
				    pow(n_e_a_vector_control(1)*jac_epsilon_control->GetBinError(i),2) +
				    pow(J_control(position_to_fill,1),2)*V_epsilon_control(0,0) +
				    pow(n_e_a_vector_control(0)*mDh_proj_diff_mc->GetBinError(i),2) + 
				    pow(J_control(position_to_fill,0),2)*V_nu_control(0,0) +
				    pow(n_e_a_vector_control(2)*jac_alpha_control->GetBinError(i),2) + 
				    pow(J_control(position_to_fill,2),2)*V_alpha_control(0,0)
				,0.5);
	      
	      corrected_diff->SetBinContent(i, value_corrected);
	      corrected_diff->SetBinError(i, error_corrected);
	      
	      position_to_fill++;
	    }
	  }
	  if (position_to_fill != filled_bins_mll_diff){ std::cout<<"PROBLEM: counting corrected_diff bins \n"; }
	  
	  // Draw corrected_diff in diff fit panel
	  c1->cd();
	  c1->cd(1);
	  
	  fitresult = fitHisto(corrected_diff, 1, 1, 5);
	  
	  // Save for pull distribution epsilon_control
	  mean_corrected_diff = fitresult[0];
	  sigma_corrected_diff = fitresult[1];
	  error_mean_corrected_diff = fitresult[2];
	  error_sigma_corrected_diff = fitresult[3];
	  corrected_diff->SetLineColor(kBlack);
	  corrected_diff->Draw("HIST SAME");
	  
	  // Draw pull distributions
	  
	  // nu diff
	  pull_nu_control->Fill( (integral_diff_data - integral_mc - integral_mc*n_e_a_vector_control(0)) / integral_mc / pow(V_nu_control(0,0),0.5) );
	  
	  // epsilon diff
	  pull_epsilon_control->Fill( (mean_diff_data - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) );
	  pull_epsilon_control1->Fill( (mean_diff_data - mean_mc - e_a_vector_control(0))/pow(V_epsilon_control(0,0),0.5) );
	  pull_epsilon_control2->Fill( (mean_diff_data - mean_mc - e_vector_control(0))/pow(V_epsilon_control(0,0),0.5) );
	  
	  // epsilon test
	  epsilon_test1->SetBinError(epsilon_test1->Fill(name.c_str(), (mean_diff_data - mean_mc) / sigma_mc), 0.001);
	  if ( abs((mean_diff_data - mean_mc) / sigma_mc) > 0.05 ){
	    epsilon_test2->Fill(abs((mean_diff_data - mean_mc) / sigma_mc), (mean_diff_data - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5));
	  }
	  
	  // alpha diff
	  pull_alpha_control->Fill( (sigma_diff_data - sigma_mc - sigma_mc*n_e_a_vector_control(2))/ sigma_mc / pow(V_alpha_control(0,0),0.5) );
	  
	  //--------------------------------------------------------------------------
	  // Add to diff panel
	  leg1->AddEntry(corrected_diff, "corrected diff", "l");
	  leg_entry = "Pull #varepsilon=" + to_string( (mean_diff_data - mean_mc - n_e_a_vector_control(1))/pow(V_epsilon_control(0,0),0.5) ).substr(0, 6);
	  leg1->AddEntry((TObject*)0, leg_entry.c_str(), "");
	}
	
	// Finish drawing diff panel
	c1->cd();
	c1->cd(1);
	leg1->Draw("");
	
	c1->cd();
	// Write plots per k bin
	f_fits->WriteObject(c1, name.c_str());
	
	// Fill control TTree
	control_tree->Fill();
	//std::cout<<"Filled tree in bin "<< name.c_str() <<" \n";
	
      }
      
    }
  }

  f_pass_reg.close();
  
  hfrac = empty_histos_count / all_histos_count;
  efrac = remaining_nevents / total_nevents;
  
  std::cout <<"For nbinspt="<<nbinspt<<" there are " << total_nevents << " events in the inclusive projection \n";
  std::cout<<"CHECKPOINT: "<<stringify_name(nbinseta, nbinspt, nbinseta, nbinspt)<<" histos  empty/all="<< hfrac <<"; events remaining/all "<< efrac <<"\n";
  
  f_control_tree->WriteObject(control_tree, "control_tree");
  
  //--------------------------------------------------------------------------
  // Write remaining histograms 
  
  gaus_integral->SetStats(0);
  gaus_integral->LabelsDeflate();
  f_control->WriteObject(gaus_integral, "gaus_integral");
  
  occupancy->SetStats(0);
  occupancy->LabelsDeflate();
  f_control->WriteObject(occupancy, "bin_occupancy");

  mean_diff->SetStats(0);
  mean_diff->LabelsDeflate();
  f_control->WriteObject(mean_diff, "mean_diff");

  sigma_diff->SetStats(0);
  sigma_diff->LabelsDeflate();
  f_control->WriteObject(sigma_diff, "sigma_diff");

  mean_diff_by_idx->SetStats(0);
  f_control->WriteObject(mean_diff_by_idx, "mean_diff_by_idx");

  sigma_diff_by_idx->SetStats(0);
  f_control->WriteObject(sigma_diff_by_idx, "sigma_diff_by_idx");
  
  f_control->WriteObject(alpha, "alpha");
  f_control->WriteObject(jac_alpha_inclusive, "jacobian_alpha_inclusive");

  f_control->WriteObject(beta, "beta");
  f_control->WriteObject(jac_beta_inclusive, "jacobian_beta_inclusive");

  f_control->WriteObject(nu, "nu");

  f_control->WriteObject(h_scales, "h_scales");
  f_control->WriteObject(h_masks, "h_masks");
  
  // Validation plots
  if (mode_option.compare("validation") == 0){
    f_control->WriteObject(epsilon_control, "epsilon_control");
    fitresult = fitHisto(pull_epsilon_control, 1, 1, 5);
    f_control->WriteObject(pull_epsilon_control, "pull_epsilon_control");
    
    f_control->WriteObject(epsilon_test1, "epsilon_test1");
    
    f_control->WriteObject(alpha_control, "alpha_control");
    fitresult = fitHisto(pull_alpha_control, 1, 1, 5);
    f_control->WriteObject(pull_alpha_control, "pull_alpha_control");
    
    f_control->WriteObject(nu_control, "nu_control");
    fitresult = fitHisto(pull_nu_control, 1, 1, 5);
    f_control->WriteObject(pull_nu_control, "pull_nu_control");

    //--------------------------------------------------------------------
    // Superimposed histograms

    // pull_epsilon_many

    TCanvas *c4 = new TCanvas("c4","c4",800,600);
    auto leg4 = new TLegend(0.68, 0.78, 0.90, 0.90);
    
    leg4->SetFillStyle(0);
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.025);
    leg4->SetFillColor(10);
    leg4->SetNColumns(1);
    leg4->SetHeader("");
    
    pull_epsilon_control->SetMaximum( pull_epsilon_control->GetBinContent(pull_epsilon_control->GetMaximumBin())*1.3 );
    pull_epsilon_control->SetStats(0);
    fitresult = fitHisto(pull_epsilon_control, 1, 1, 5);
    leg_entry = "freeze none, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
    leg4->AddEntry(pull_epsilon_control, leg_entry.c_str(), "l");
    pull_epsilon_control->SetLineColor(kBlack);
    pull_epsilon_control->Draw("E");
    /*
      fitresult = fitHisto(pull_epsilon_control1, 1, 2, 5);
      leg_entry = "freeze nu, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
      leg4->AddEntry(pull_epsilon_control1, leg_entry.c_str(), "l");
      pull_epsilon_control1->SetLineColor(kRed);
      pull_epsilon_control1->Draw("SAME");
      
      fitresult = fitHisto(pull_epsilon_control2, 1, 3, 5);
      leg_entry = "freeze nu, alpha, #mu = " + to_string(fitresult[0]).substr(0, 6) + ", #sigma= " + to_string(fitresult[1]).substr(0, 6);
      leg4->AddEntry(pull_epsilon_control2, leg_entry.c_str(), "l");
      pull_epsilon_control2->SetLineColor(kGreen);
      pull_epsilon_control2->Draw("SAME");
      
    */
    leg4->Draw("");
    
    f_control->WriteObject(c4, "pull_epsilon_many");
    
    //-----------------------------------------------------------
    // Other histograms
    
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    epsilon_test2->SetStats(0);
    epsilon_test2->Draw("COLZ");
    f_control->WriteObject(c3, "epsilon_test2");
  }
  
  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Dataframe control histograms
  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
      //Save tree for debugging
      //std::unique_ptr<TFile> f1( TFile::Open("InOutputFiles/snapshot_output.root", "RECREATE") );
      //dlast->Snapshot("Events", "snapshot_output.root", {"GenPart_status", "GenPart_pt", "posPtSmearBetaVal", "negPtSmearBetaVal", "Muon_charge", "GenPart_pdgId", "GenPart_genPartIdxMother"});
      //f1->Close();
        
 /*
      std::unique_ptr<TFile> f3( TFile::Open("InOutputFiles/control_bin_histo.root", "RECREATE") );
      f3->WriteObject(pt_pos_uni.GetPtr(), "pt_pos_uni");
      
      //TH1 in pT+ with variable bin size -> should be uniform
      auto pt_pos = d_mc->Histo1D({"pt_pos", "pt mu+", nbinspt, myptboundaries},"posTrackPt","weight");
      f3->WriteObject(pt_pos.GetPtr(), "pt_pos");
      
      auto pt_eta_pos = d_mc->Histo2D({"pt_eta_pos", "pt eta mu+", nbinseta, myetaboundaries, nbinspt, myptboundaries},"posTrackEta", "posTrackPt", "weight");
      f3->WriteObject(pt_eta_pos.GetPtr(), "pt_eta_pos");
      
      // small uniform pt binning
      vector<double> smallptunibinning;
      for (int i=0; i<=nbinspt*10; i++){
      smallptunibinning.push_back(pt_low + i * (pt_high - pt_low)/nbinspt/10);
      }
      
      double mysmallptunibinning[smallptunibinning.size()];
      for (int i=0; i<smallptunibinning.size(); i++){
      mysmallptunibinning[i] = smallptunibinning[i];
      }
      
      auto pt_pos_bin = d_mc->Histo2D({"pt_pos_bin", "pt mu+ distr in bin", nbinspt*10, mysmallptunibinning, nbinspt, myptboundaries}, "posTrackPt", "posTrackPt", "weight");
      f3->WriteObject(pt_pos_bin.GetPtr(), "pt_pos_bin");
      
      auto pt_neg_bin = d_mc->Histo2D({"pt_neg_bin", "pt mu- distr in bin", nbinspt*10, mysmallptunibinning, nbinspt, myptboundaries}, "negTrackPt", "negTrackPt", "weight");
      f3->WriteObject(pt_neg_bin.GetPtr(), "pt_neg_bin");
      
      f3->Close();
  */
  /*
      std::unique_ptr<TFile> f4( TFile::Open("InOutputFiles/multiD_histo_smear_beta_val_easy.root", "RECREATE") );
      // mll_diff_smear_beta_val_easy 
      auto mDh_diff_smear_beta_val_easy = d_sim_data->HistoND<float, float, float, float, float, double>({"multi_data_histo_diff_smear_beta_val_easy", "multi_data_histo_diff_smear_beta_val_easy", 5, {nbinseta, nbinspt, nbinseta, nbinspt, nbinsmll_diff}, {etabinranges, ptbinranges, etabinranges, ptbinranges, mll_diffbinranges}}, {"posTrackEta","posPtSmear","negTrackEta","negPtSmear","mll_diff_smear","weight"});
      f4->WriteObject(mDh_diff_smear_beta_val_easy.GetPtr(), "multi_data_histo_diff_smear_beta_val_easy");
      f4->Close();
  */
  
  //} closes if (mode_option.compare("simulation") == 0) {


  //---------------------------------------------------------------------------------------------------------------------------------------------------
  // Fitting constants part
  //---------------------------------------------------------------------------------------------------------------------------------------------------

  // Set verbosity
  int verbosity = 2;
  ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);

  //-------------------------------------------------
  // Internal to physical parameter transformation

  TMatrixD V_internal_to_physical(nbinseta, nbinseta);
  TArrayD V_internal_to_physical_elements(nbinseta*nbinseta);
  
  //identity e.g. M not decorrelated
  for(int j=0; j<nbinseta; j++){
    for(int i=0; i<nbinseta; i++){
      if (j == i){
        V_internal_to_physical_elements[j*nbinseta+i] = 1.0;
      } else {
        V_internal_to_physical_elements[j*nbinseta+i] = 0.0;
      }
    }
  }
  
  V_internal_to_physical.SetMatrixArray(V_internal_to_physical_elements.GetArray());
  
  //std::cout<< "\n" << "VIntToPhys"<< "\n";
  //V_internal_to_physical.Print();

  //-------------------------------------------------
  // Physical to internal parameter transformation
  
  TMatrixD V_physical_to_internal(nbinseta, nbinseta);
  TArrayD V_physical_to_internal_elements(nbinseta*nbinseta);
  
  //identity e.g. M not decorrelated
  for(int j=0; j<nbinseta; j++){
    for(int i=0; i<nbinseta; i++){
      if (j == i){
	V_physical_to_internal_elements[j*nbinseta+i] = 1.0;
      } else {
	V_physical_to_internal_elements[j*nbinseta+i] = 0.0;
      }
    }
  }
  
  V_physical_to_internal.SetMatrixArray(V_physical_to_internal_elements.GetArray());
  //std::cout<< "\n" << "VPhysToInt"<< "\n";
  //V_physical_to_internal.Print();
  
  unsigned int n_data_points, n_parameters;
  vector<string> labels;
  vector<double> scale_squared_values, scale_squared_error_values;

  //-----------------------------------------------
  // Get data
    
  int n_data_points_initial = beta->GetEntries();
  n_data_points = 0;

  std::cout << "\n" << n_data_points_initial <<" entries initially in beta histogram";
  
  for(int i=0; i<n_data_points_initial; i++){
    if( occupancy->GetBinContent(i+1) > 300.0 && gaus_integral->GetBinContent(i+1) > 0.7 ){
      scale_squared_values.push_back((1.0 + beta->GetBinContent(i+1))*(1.0 + beta->GetBinContent(i+1)));
      scale_squared_error_values.push_back(2*abs(1.0 + beta->GetBinContent(i+1))*(beta->GetBinError(i+1)));
      labels.push_back(beta->GetXaxis()->GetLabels()->At(i)->GetName());
      
      n_data_points += 1;
    }
  }
  
  for(int i=0; i<10; i++){
    std::cout << "\n" << labels[i] << ": beta = " << beta->GetBinContent(i+1) << "; scale_squared =  " << scale_squared_values[i] << "; scale_squared_error = " << scale_squared_error_values[i] << "\n";
  }
  
  std::cout << "\n" << n_data_points << "entries remain from beta histogram";
  std::cout << "\n" << "scale_squared_values.size() = " << scale_squared_values.size() << "; scale_squared_error_values.size() = " << scale_squared_error_values.size() << "labels.size() = " << labels.size();

  n_parameters = 3*nbinseta; // 3 for A,e,M model
  
  //-------------------------------------------------
  // Use data

  // TheoryFcn(const vector<double> &meas, const vector<double> &err, const vector<string> &bin_labels, const TMatrixD &V_phys_to_int, const TMatrixD &V_int_to_phys)
  TheoryFcn fFCN(scale_squared_values, scale_squared_error_values, labels, V_physical_to_internal, V_internal_to_physical);
  fFCN.SetErrorDef(1.0/(n_data_points - n_parameters)); // new error definition when minimising reduced chi2

  // Create parameters with initial starting values

  double start=0.0, par_error=0.01; // will store error on parameter before taking scaling into account
  MnUserParameters upar;

  for (int i=0; i<n_parameters/3; i++){ 
    upar.Add(Form("param%d",i), start, par_error); // A_internal
    //upar.Add(Form("param%d",i), start); // A_internal const
  }
  for (int i=n_parameters/3; i<2*n_parameters/3; i++){ // e_internal
    upar.Add(Form("param%d",i), start, par_error);
    //upar.Add(Form("param%d",i), start); // e_internal set constant to 0
  }
  for (int i=2*n_parameters/3; i<n_parameters; i++){ // M_internal
    upar.Add(Form("param%d",i), start, par_error); //upar.Add(Form("param%d",i), double((i-2*nbinseta)/10.0), par_error);
  }

  // create Migrad minimizer

  MnMigrad minimize(fFCN, upar, 1); //TODO strategy is 1, check others too

  // ... and Minimize

  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  
  double tolerance(0.001); //MIGRAD will stop iterating when edm : 0.002 * tolerance * UPERROR
  //fFCN.SetErrorDef(1.0/(n_data_points - n_parameters)); 
  //3829.158
    
  FunctionMinimum min = minimize(maxfcn, tolerance);
    
  // save results
  
  double k_middle = (1.0/55.0 + 1.0/25.0)/2.0; //TODO change to not input by hand
  vector<double> A_e_M_values(n_parameters), A_e_M_errors(n_parameters);
  double cov_prime;
  
  TMatrixD internal_M_matrix(nbinseta,1), internal_M_err_matrix(nbinseta,1), physical_M_matrix(nbinseta,1), physical_M_err_matrix(nbinseta,1);
  TArrayD internal_M_data(nbinseta), internal_M_error_data(nbinseta);
  vector<double> physical_M_parameters(nbinseta), physical_M_errors(nbinseta);
  
  for(int i=0; i<nbinseta; i++){
    internal_M_data[i] = min.UserState().Value(2*nbinseta + i);
    internal_M_error_data[i] = min.UserState().Error(2*nbinseta + i);
  }
    
  internal_M_matrix.SetMatrixArray(internal_M_data.GetArray());
  physical_M_matrix = TMatrixD(V_internal_to_physical,TMatrixD::kMult,internal_M_matrix);
  
  internal_M_err_matrix.SetMatrixArray(internal_M_error_data.GetArray());
  physical_M_err_matrix = TMatrixD(V_internal_to_physical,TMatrixD::kMult,internal_M_err_matrix);
  
  for(int i=0; i<nbinseta; i++){
    physical_M_parameters[i] = TMatrixDColumn(physical_M_matrix,0)(i);
    physical_M_errors[i] =TMatrixDColumn(physical_M_err_matrix,0)(i);
  }

  for (int i=0; i<n_parameters; i++){ 
    int whole = i / nbinseta;

    //term_pos = (1. + scaling_A*par[eta_pos_index] - scaling_e_prime*par[eta_pos_index + nbinseta]*(k_plus - k_middle) + scaling_M*physical_M_pars[eta_pos_index]/k_plus);
    if (whole == 0) { // A = scaling_A_prime*par[A_prime] + scaling_e_prime*par[e_prime]*k_middle
      cov_prime = min.UserState().Covariance().Data()[(i+nbinseta)*(i+nbinseta+1)/2.0 + i];
      //corr_hist->SetBinContent(bin, covariance[(i_temp-1)*(i_temp)/2+(j_temp-1)] / min.UserState().Error(i-1) / min.UserState().Error(j-1));
      A_e_M_values[i] = min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i)) + min.UserState().Value(i+nbinseta) * get<1>(getParameterNameAndScaling(i+nbinseta)) * k_middle;
      A_e_M_errors[i] = pow(pow(get<1>(getParameterNameAndScaling(i)),2)*pow(min.UserState().Error(i),2)+pow(get<1>(getParameterNameAndScaling(i+nbinseta))*k_middle,2)*pow(min.UserState().Error(i+nbinseta),2)+2.0*get<1>(getParameterNameAndScaling(i))*get<1>(getParameterNameAndScaling(i+nbinseta))*k_middle*cov_prime, 0.5); // with cov[A',e']
    }
    else if (whole == 1) { // e = scaling_e_prime*par[e_prime]
      A_e_M_values[i] = min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i));
      A_e_M_errors[i] = min.UserState().Error(i) * abs(get<1>(getParameterNameAndScaling(i)));
    }
    else if (whole == 2) { // M
      //cout << i<<" " << physical_M_parameters[i - 2*nbinseta] * get<1>(getParameterNameAndScaling(i)) << "\n";
      A_e_M_values[i] = physical_M_parameters[i - 2*nbinseta] * get<1>(getParameterNameAndScaling(i)); //min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i));
      //cout << i<<" " << A_e_M_values[i] << "\n";
      A_e_M_errors[i] = physical_M_errors[i - 2*nbinseta] * abs(get<1>(getParameterNameAndScaling(i))); //min.UserState().Error(i) *get<1>(getParameterNameAndScaling(i));
    }
  }

  cout << "\n" << "chi^2/ndf: " << min.Fval() + 1. << "\n";
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
  
  cout << "\tHesse..." << endl;
  MnHesse hesse(1);
  hesse(fFCN, min);
  
  vector<double> hessian = min.UserState().Hessian().Data();
  vector<double> covariance = min.UserState().Covariance().Data();
  
  TH2D *hessian_hist = new TH2D("hessian_hist", "Hessian", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);
  TH2D *covariance_hist = new TH2D("covariance_hist", "Covariance", n_parameters, 0, n_parameters, 0, n_parameters);
  TH2D *corr_hist = new TH2D("corr_hist", "Correlation", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);

  TH2D *corr_physical_hist = new TH2D("corr_physical_hist", "Correlation physical parameters", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);

  int bin, i_temp, j_temp, bin_x, bin_y;
  for (int i=1; i<=n_parameters; i++){ // counts x axis bins left to right
    hessian_hist->GetXaxis()->SetBinLabel(i,get<0>(getParameterNameAndScaling(i-1)).c_str());
    hessian_hist->GetYaxis()->SetBinLabel(i,get<0>(getParameterNameAndScaling(n_parameters-i)).c_str());
    covariance_hist->GetXaxis()->SetBinLabel(i,get<0>(getParameterNameAndScaling(i-1)).c_str());
    covariance_hist->GetYaxis()->SetBinLabel(i,get<0>(getParameterNameAndScaling(n_parameters-i)).c_str());
    corr_hist->GetXaxis()->SetBinLabel(i,get<0>(getParameterNameAndScaling(i-1)).c_str());
    corr_hist->GetYaxis()->SetBinLabel(i,get<0>(getParameterNameAndScaling(n_parameters-i)).c_str());
    for (int j=1; j<=i; j++){ // counts y axis bins, n_parameters+1-j counts bins up to down
      //if ((i<=nbinseta || i>nbinseta*2) && (j<=nbinseta || j>nbinseta*2)){ //TODO remove if and else block when fitting 3 pars
      bin = hessian_hist->GetBin(i, n_parameters+1-j);
      i_temp = i; // i_temp = (i > nbinseta*2) ? i - nbinseta : i; //TODO when fitting 3 pars i_temp can be just i
      j_temp = j; // j_temp = (j > nbinseta*2) ? j - nbinseta : j; //TODO when fitting 3 pars j_temp can be just j
      hessian_hist->SetBinContent(bin, hessian[(i_temp-1)*(i_temp)/2+(j_temp-1)]); 
      covariance_hist->SetBinContent(bin, covariance[(i_temp-1)*(i_temp)/2+(j_temp-1)]); 
      corr_hist->SetBinContent(bin, covariance[(i_temp-1)*(i_temp)/2+(j_temp-1)] / min.UserState().Error(i-1) / min.UserState().Error(j-1)); // TODO when fitting 2 pars is i i_temp, j j_temp?
      
      bin = hessian_hist->GetBin(j, n_parameters+1-i); 
      
      hessian_hist->SetBinContent(bin, hessian[(i_temp-1)*(i_temp)/2+(j_temp-1)]); 
      covariance_hist->SetBinContent(bin, covariance[(i_temp-1)*(i_temp)/2+(j_temp-1)]); 
      corr_hist->SetBinContent(bin, covariance[(i_temp-1)*(i_temp)/2+(j_temp-1)] / min.UserState().Error(i-1) / min.UserState().Error(j-1)); // TODO when fitting 2 pars is i i_temp, j j_temp?
      
      //} else {
      //bin = hessian_hist->GetBin(i, n_parameters+1-j); 
      //hessian_hist->SetBinContent(bin, 0.0);
      //covariance_hist->SetBinContent(bin, 0.0);
      //corr_hist->SetBinContent(bin, 0.0);
      //bin = hessian_hist->GetBin(j, n_parameters+1-i); 
      //hessian_hist->SetBinContent(bin, 0.0);
      //covariance_hist->SetBinContent(bin, 0.0);
      //corr_hist->SetBinContent(bin, 0.0);
      // }
    }
  }
  
  double hessian_ar[n_parameters*n_parameters]; 
  double covariance_ar[n_parameters*n_parameters];
  double corr_ar[n_parameters*n_parameters];

  for (int i=1; i<=n_parameters; i++){
    for (int j=1; j<=n_parameters; j++){
      // define arrays collumn wise
      hessian_ar[(i-1)*n_parameters+(j-1)] = hessian_hist->GetBinContent(hessian_hist->GetBin(i, n_parameters+1-j)); 
      covariance_ar[(i-1)*n_parameters+(j-1)] = covariance_hist->GetBinContent(covariance_hist->GetBin(i, n_parameters+1-j)); 
      corr_ar[(i-1)*n_parameters+(j-1)] = corr_hist->GetBinContent(corr_hist->GetBin(i, n_parameters+1-j));
      //TODO scaling???
    }
  }
  
  //cout << "Hessian: " << hessian <<"\n"<< "Hessian: " <<"\n";
  
  TMatrixTSym<double> Hessian_matrix(n_parameters, hessian_ar, "F"); //need option F to unroll collumn wise //TODO pass by ptr
  
  for(int i=0; i<n_parameters; i++){ 
    for(int j=0; j<n_parameters; j++){
      //cout << Hessian_matrix(i,j) << " ";
    }
    //cout << "\n";
  }
  //cout << "\n";
  //cout << "Covariance: " << covariance << "\n" << "Covariance: " << "\n";
 
  TMatrixTSym<double> Covariance_matrix(n_parameters, covariance_ar, "F"); //need option F to unroll collumn wise //TODO pass by ptr

  for(int i=0; i<n_parameters; i++){
    for(int j=0; j<n_parameters; j++){
      //cout << Covariance_matrix(i,j) << " ";
    }
    //cout << "\n";
  }

  TMatrixD U_internal_to_physical(n_parameters, n_parameters); // TODO what if n_parameters != 3*netabins?
  TArrayD U_internal_to_physical_elements(n_parameters*n_parameters);

  for(int j=0; j<nbinseta*3; j++){ // counts rows 
    if (j<nbinseta){ // A rows
      for(int i=0; i<nbinseta*3; i++){ 
	if(i == j){
	  U_internal_to_physical_elements[j*nbinseta*3+i] = get<1>(getParameterNameAndScaling(i));  // A collumns non empty
	} else if (i >= nbinseta && i == j + nbinseta && i < nbinseta*2){
	  U_internal_to_physical_elements[j*nbinseta*3+i] = get<1>(getParameterNameAndScaling(i))*k_middle; // e collumns non empty
	} else {
	  U_internal_to_physical_elements[j*nbinseta*3+i] = 0.0; // M collumns and the rest
	}
      }
    } else if(j<2*nbinseta){ // e rows
      for(int i=0; i<nbinseta*3; i++){
        if(i == j){
          U_internal_to_physical_elements[j*nbinseta*3+i] = get<1>(getParameterNameAndScaling(i));  // e collumns non empty
        } else {
          U_internal_to_physical_elements[j*nbinseta*3+i] = 0.0; // the rest
        }
      }
    } else if(j<3*nbinseta){ // M rows
      for(int i=0; i<nbinseta*3; i++){
        if(i == j){
          U_internal_to_physical_elements[j*nbinseta*3+i] = get<1>(getParameterNameAndScaling(i));  // M collumns non empty // TODO this is insufficient if you decorrelate M
        } else {
          U_internal_to_physical_elements[j*nbinseta*3+i] = 0.0; // the rest
        }
      }
    } else {std::cout << "\n" << "Miscounted U";}
  }

  U_internal_to_physical.SetMatrixArray(U_internal_to_physical_elements.GetArray());
  //U_internal_to_physical.Print();

  TMatrixD left_term(n_parameters, n_parameters), covariance_matrix_physical(n_parameters, n_parameters);
  
  left_term = TMatrixD(U_internal_to_physical,TMatrixD::kMult,Covariance_matrix);
  covariance_matrix_physical = TMatrixD(left_term, TMatrixD::kMultTranspose,U_internal_to_physical);
  
  //TMatrixD U_internal_to_physical_transpose(TMatrixD::kTransposed,U_internal_to_physical);
  //left_term = TMatrixD(U_internal_to_physical_transpose, TMatrixD::kMult, Covariance_matrix);
  //covariance_matrix_physical = TMatrixD(left_term, TMatrixD::kMult,U_internal_to_physical);
  
  for(int i=1; i<=n_parameters; i++){
    corr_physical_hist->GetXaxis()->SetBinLabel(i,get<2>(getParameterNameAndScaling(i-1)).c_str());
    corr_physical_hist->GetYaxis()->SetBinLabel(i,get<2>(getParameterNameAndScaling(n_parameters-i)).c_str());
    for(int j=1; j<=n_parameters; j++){
      bin = corr_physical_hist->GetBin(i, n_parameters-j+1);
      corr_physical_hist->SetBinContent(bin, covariance_matrix_physical(i-1,j-1) / A_e_M_errors[i-1] / A_e_M_errors[j-1]); 
    }
  }

  //cout << "Correlation: " << "\n";

  TMatrixTSym<double> Corr_matrix(n_parameters, corr_ar, "F"); //need option F to unroll collumn wise //TODO pass by ptr

  for(int i=0; i<n_parameters; i++){ 
    for(int j=0; j<n_parameters; j++){
      //cout << Corr_matrix(i,j) << " ";
    }
    //cout << "\n";
  }
  
  if (iter == total_iter-1){
    unique_ptr<TFile> f_a( TFile::Open("InOutputFiles/constants_fitted_correlation.root", "RECREATE") );
    
    TCanvas *c4 = new TCanvas("c4","c4",800,600);
    hessian_hist->SetStats(0);
    hessian_hist->Draw("COLZ");
    f_a->WriteObject(c4, "hessian_hist");

    TCanvas *c5 = new TCanvas("c5","c5",800,600);
    covariance_hist->SetStats(0);
    covariance_hist->Draw("COLZ");
    f_a->WriteObject(c5, "covariance_hist");

    corr_hist->SetStats(0);
    corr_hist->Draw("COLZ text");
    f_a->WriteObject(c5, "corr_hist");

    corr_physical_hist->SetStats(0);
    corr_physical_hist->Draw("COLZ text");
    f_a->WriteObject(c5, "corr_physical_hist");
  }
  //cout << "\n" << "Fitted A' [ ], e' [GeV], M [GeV^-1] parameters: " << "\n";
  //for (unsigned int i(0); i<upar.Params().size(); i++) {
  //cout <<"par[" << i << "]: " << get<0>(getParameterNameAndScaling(i)) << " fitted to: "<< min.UserState().Value(i) * get<1>(getParameterNameAndScaling(i)) << " +/- " << min.UserState().Error(i) * abs(get<1>(getParameterNameAndScaling(i))) << "\n";
  //}

  // Histogram with dummy parameters and fitted parameters

  TH1D *dummy_pars_A = new TH1D("dummy_pars_A", "A", nbinseta, -2.4, 2.4);
  TH1D *dummy_pars_e = new TH1D("dummy_pars_e", "e", nbinseta, -2.4, 2.4);
  TH1D *dummy_pars_M = new TH1D("dummy_pars_M", "M", nbinseta, -2.4, 2.4);
  
  TH1D *fitted_pars_A = new TH1D("fitted_pars_A", "A", nbinseta, -2.4, 2.4);
  TH1D *fitted_pars_e = new TH1D("fitted_pars_e", "e", nbinseta, -2.4, 2.4);
  TH1D *fitted_pars_M = new TH1D("fitted_pars_M", "M", nbinseta, -2.4, 2.4);

  TH1D *pull_A = new TH1D("pull_A", "A", nbinseta, -2.4, 2.4);
  TH1D *pull_e = new TH1D("pull_e", "e", nbinseta, -2.4, 2.4);
  TH1D *pull_M = new TH1D("pull_M", "M", nbinseta, -2.4, 2.4);

  for (int i=1; i<=nbinseta; i++){
    fitted_pars_A->SetBinContent(i, A_e_M_values[i-1]); 
    fitted_pars_A->SetBinError(i, A_e_M_errors[i-1]);
    A_from_iter[i-1] = A_e_M_values[i-1];

    fitted_pars_e->SetBinContent(i, A_e_M_values[i-1+nbinseta]);
    fitted_pars_e->SetBinError(i, A_e_M_errors[i-1+nbinseta]);
    e_from_iter[i-1] = A_e_M_values[i-1+nbinseta];

    fitted_pars_M->SetBinContent(i, A_e_M_values[i-1+2*nbinseta]);
    fitted_pars_M->SetBinError(i, A_e_M_errors[i-1+2*nbinseta]);
    M_from_iter[i-1] = A_e_M_values[i-1+2*nbinseta];

    dummy_pars_A->SetBinContent(i, A_values[i-1]);
    dummy_pars_A->SetBinError(i, 1e-10);
    dummy_pars_e->SetBinContent(i, e_values[i-1]);
    dummy_pars_e->SetBinError(i, 1e-10);
    dummy_pars_M->SetBinContent(i, M_values[i-1]);
    dummy_pars_M->SetBinError(i, 1e-10);

    pull_A->SetBinContent(i, (A_e_M_values[i-1] - A_values[i-1])/A_e_M_errors[i-1]);
    pull_A->SetBinError(i, 1e-10);
    pull_e->SetBinContent(i, (A_e_M_values[i-1+nbinseta] - e_values[i-1])/A_e_M_errors[i-1+nbinseta]);
    pull_e->SetBinError(i, 1e-10);
    pull_M->SetBinContent(i, (A_e_M_values[i-1+2*nbinseta] - M_values[i-1])/A_e_M_errors[i-1+2*nbinseta]);
    pull_M->SetBinError(i, 1e-10);
    
  }

  unique_ptr<TFile> f_for_plotting( TFile::Open("InOutputFiles/constants_fitted_for_plotting.root", "RECREATE") );
  f_for_plotting->WriteObject(fitted_pars_A, "fitted_pars_A");
  f_for_plotting->WriteObject(fitted_pars_e, "fitted_pars_e");
  f_for_plotting->WriteObject(fitted_pars_M, "fitted_pars_M");
  f_for_plotting->WriteObject(dummy_pars_A, "dummy_pars_A");
  f_for_plotting->WriteObject(dummy_pars_e, "dummy_pars_e");
  f_for_plotting->WriteObject(dummy_pars_M, "dummy_pars_M");
  f_for_plotting->WriteObject(pull_A, "pull_A");
  f_for_plotting->WriteObject(pull_e, "pull_e");
  f_for_plotting->WriteObject(pull_M, "pull_M");

  TCanvas *c1_c = new TCanvas("c1_c","c1_c",800,600);

  auto leg1_c = new TLegend(0.68, 0.78, 0.90, 0.90);
    
  leg1_c->SetFillStyle(0);
  leg1_c->SetBorderSize(0);
  leg1_c->SetTextSize(0.025);
  leg1_c->SetFillColor(10);
  leg1_c->SetNColumns(1);
  leg1_c->SetHeader("");
  
  //fitted_pars_A->SetMinimum(-0.00025);
  //fitted_pars_A->SetMaximum(0.0007);
  fitted_pars_A->GetYaxis()->SetNoExponent();
  fitted_pars_A->SetStats(0);
  fitted_pars_A->GetXaxis()->SetTitle("#eta");
  fitted_pars_A->GetYaxis()->SetTitle("Magnetic field correction");
  fitted_pars_A->Draw("");
  
  dummy_pars_A->SetLineColor(kRed);
  dummy_pars_A->SetMarkerStyle(kPlus);
  dummy_pars_A->SetMarkerColor(kRed);
  dummy_pars_A->Draw("SAME");
    
  leg1_c->AddEntry(dummy_pars_A, "input par", "l");
    
  leg1_c->AddEntry(fitted_pars_A, "fitted par", "l");
  leg1_c->Draw("SAME");

  f_control_const->WriteObject(c1_c, ("A_iter_" + to_string(iter)).c_str());
  leg1_c->Clear();
  c1_c->Clear();
    
  //TCanvas *c01_c = new TCanvas("c01_c","c01_c",800,600);

  pull_A->GetYaxis()->SetNoExponent();
  pull_A->SetStats(0);
  pull_A->GetXaxis()->SetTitle("#eta");
  pull_A->GetYaxis()->SetTitle("Pull");
  pull_A->Draw("");

  f_control_const->WriteObject(c1_c, ("pull_A_iter_" + to_string(iter)).c_str());
  c1_c->Clear();
  
  //TCanvas *c2_c = new TCanvas("c2_c","c2_c",800,600);

  //fitted_pars_e->SetMinimum(-0.001);
  //fitted_pars_e->SetMaximum(0.017);
  fitted_pars_e->GetYaxis()->SetNoExponent();
  fitted_pars_e->SetStats(0);
  fitted_pars_e->GetXaxis()->SetTitle("#eta");
  fitted_pars_e->GetYaxis()->SetTitle("Material correction (GeV)");
  fitted_pars_e->Draw("");


  dummy_pars_e->SetLineColor(kRed);
  dummy_pars_e->SetMarkerStyle(kPlus);
  dummy_pars_e->SetMarkerColor(kRed);
  dummy_pars_e->Draw("SAME");

  leg1_c->AddEntry(dummy_pars_e, "input par", "l");
  
  leg1_c->AddEntry(fitted_pars_e, "fitted par", "l");
  leg1_c->Draw("SAME");
  
  f_control_const->WriteObject(c1_c, ("e_iter_" + to_string(iter)).c_str());
  leg1_c->Clear();
  c1_c->Clear();
  
  //TCanvas *c02_c = new TCanvas("c02_c","c02_c",800,600);

  pull_e->GetYaxis()->SetNoExponent();
  pull_e->SetStats(0);
  pull_e->GetXaxis()->SetTitle("#eta");
  pull_e->GetYaxis()->SetTitle("Pull");
  pull_e->Draw("");
  
  f_control_const->WriteObject(c1_c, ("pull_e_iter_" + to_string(iter)).c_str());
  c1_c->Clear();

  
  //TCanvas *c3_c = new TCanvas("c3_c","c3_c",800,600);

  //fitted_pars_M->SetMinimum(-3e-5);
  fitted_pars_M->SetMaximum(4.5e-5);
  //fitted_pars_M->SetMaximum(fitted_pars_M->GetBinContent(fitted_pars_M->GetMaximumBin())*1.1);
  fitted_pars_M->SetStats(0);
  fitted_pars_M->GetXaxis()->SetTitle("#eta");
  fitted_pars_M->GetYaxis()->SetTitle("Misalignment correction (GeV^{-1})");
  fitted_pars_M->Draw("");

 
  dummy_pars_M->SetLineColor(kRed);
  dummy_pars_M->SetMarkerStyle(kPlus);
  dummy_pars_M->SetMarkerColor(kRed);
  dummy_pars_M->Draw("SAME");

  leg1_c->AddEntry(dummy_pars_M, "input par", "l");
 
  leg1_c->AddEntry(fitted_pars_M, "fitted par", "l");
  leg1_c->Draw("SAME");
  
  f_control_const->WriteObject(c1_c, ("M_iter_" + to_string(iter)).c_str());
  leg1_c->Clear();
  c1_c->Clear();
  
  //TCanvas *c03_c = new TCanvas("c03_c","c03_c",800,600);

  pull_M->GetYaxis()->SetNoExponent();
  pull_M->SetStats(0);
  pull_M->GetXaxis()->SetTitle("#eta");
  pull_M->GetYaxis()->SetTitle("Pull");
  pull_M->Draw("");
  f_control_const->WriteObject(c1_c, ("pull_M_iter_" + to_string(iter)).c_str());
  c1_c->Clear();
  
  cout<<"A_from_iter: "<<A_from_iter<<"\n";
  cout<<"e_from_iter: "<<e_from_iter<<"\n";
  cout<<"M_from_iter: "<<M_from_iter<<"\n";

  cout<<"Finished iter "<<iter<<"\n";
  }
  
  
  gettimeofday(&tv_stop, 0);
  t1 = (tv_stop.tv_sec - tv_start.tv_sec)*1000.0 + (tv_stop.tv_usec - tv_start.tv_usec)/1000.0;
  
  cout << "\n" << "# Calculation time: " << t1/60000.0 << " min" << "\n";  
  
  return 0;

}
  
