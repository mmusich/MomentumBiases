#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include <TMatrixD.h>
#include <TMatrixDSymfwd.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <boost/program_options.hpp>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

//#include <Eigen/Core>
//#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
using namespace ROOT::Minuit2;

typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

using namespace boost::program_options;

constexpr double MZ = 91.;
constexpr double GW = 2.5;


class TheoryFcn : public FCNGradientBase {
//class TheoryFcn : public FCNBase {

public:
  TheoryFcn(const int& debug, const int& seed, const int& bias, string fname)
    : errorDef_(1.0), debug_(debug), seed_(seed), bias_(bias)
  {

    ran_ = new TRandom3(seed);

    pt_edges_  = {25.0, 33.2011, 38.3067, 42.2411, 46.055, 55.0}; 
    eta_edges_ = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
    
    n_pt_bins_  = pt_edges_.size()-1;
    n_eta_bins_ = eta_edges_.size()-1;

    n_pars_ = 3*n_eta_bins_;
    
    for(unsigned int i = 0; i < pt_edges_.size(); i++){
      k_edges_.emplace_back(1./pt_edges_[i]);
    }
    for(unsigned int i = 0; i < n_pt_bins_; i++){
      kmean_vals_.emplace_back( 0.5*(k_edges_[i]+k_edges_[i+1]) );
    }
    kmean_val_ = 0.5*(kmean_vals_[n_pt_bins_-1] + kmean_vals_[0]);

    n_data_ = n_eta_bins_*n_eta_bins_*n_pt_bins_*n_pt_bins_;
    n_dof_ = 0;
    
    scales2_.reserve(n_data_);
    scales2Err_.reserve(n_data_);
    masks_.reserve(n_data_);
    for(unsigned int idata = 0; idata<n_data_; idata++){
      scales2_.push_back( 0.0 );
      scales2Err_.push_back( 0.0 );
      masks_.push_back(1);
    }

    x_vals_ = VectorXd(n_pars_);
    A_vals_ = VectorXd(n_eta_bins_);
    e_vals_ = VectorXd(n_eta_bins_);
    M_vals_ = VectorXd(n_eta_bins_);
    A_vals_prevfit_ = VectorXd(n_eta_bins_);
    e_vals_prevfit_ = VectorXd(n_eta_bins_);
    M_vals_prevfit_ = VectorXd(n_eta_bins_);

    int nbinseta = 24;
    
    for(int i=0; i<nbinseta; i++){ 
      A_vals_(i) = (-7.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 5.0)*0.0001;
      e_vals_(i) = (9.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 1.0)*0.001;
      M_vals_(i) = (6.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) - 2.0)*0.00001;
      A_vals_prevfit_(i) = (-7.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 5.0)*0.0001;
      e_vals_prevfit_(i) = (9.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) + 1.0)*0.001;
      M_vals_prevfit_(i) = (6.0*4.0/nbinseta/nbinseta*(i - nbinseta/2)*(i - nbinseta/2) - 2.0)*0.00001;
      x_vals_(i) = A_vals_(i);
      x_vals_(i+n_eta_bins_) = e_vals_(i);
      x_vals_(i+2*n_eta_bins_) = M_vals_(i);
      std::cout<< i <<": A: "<< x_vals_(i)<<", "<<"e: "<< x_vals_(i+n_eta_bins_)<<", "<<"M: "<< x_vals_(i+2*n_eta_bins_)<<"\n";
    }

    if(bias_>0){
      // bias for A out
      for(unsigned int i=0; i<n_eta_bins_; i++){
	double val = ran_->Uniform(-0.001, 0.001);
	if(bias_==2){
	  double mid_point = double(n_eta_bins_)*0.5;
	  val = (i-mid_point)*(i-mid_point)/mid_point/mid_point*0.001;
	}
	A_vals_(i) = val;
	x_vals_(i) = val;
      }
      // bias for e out
      for(unsigned int i=0; i<n_eta_bins_; i++){
	double val = ran_->Uniform(-0.0001/kmean_val_, 0.0001/kmean_val_);
	if(bias_==2){
	  double mid_point = double(n_eta_bins_)*0.5;
	  val = -(i-mid_point)*(i-mid_point)/mid_point/mid_point*0.0001;
	}
	e_vals_(i) = val;
	x_vals_(i+n_eta_bins_) = val;
      }
      // bias for M out
      for(unsigned int i=0; i<n_eta_bins_; i++){
	double val = ran_->Uniform(-0.001*kmean_val_, 0.001*kmean_val_);
	if(bias_==2){
	  double mid_point = double(n_eta_bins_)*0.5;
	  val = (i-mid_point)/mid_point*0.001;
	}
	M_vals_(i) = val;
	x_vals_(i+2*n_eta_bins_) = val;
      }            
    }

    else if(bias==-1){
      TFile* fin = TFile::Open(fname.c_str(), "READ");
      TH1D* h_scales = (TH1D*)fin->Get("h_scales");
      TH1D* h_masks = (TH1D*)fin->Get("h_masks");
      assert( h_scales->GetXaxis()->GetNbins() == n_data_);
      
      unsigned int n_unmasked_bins = 0;  
      for(unsigned int ibin=0;ibin<h_scales->GetXaxis()->GetNbins(); ibin++){
	scales2_[ibin]    = h_scales->GetBinContent(ibin+1)*h_scales->GetBinContent(ibin+1);
	scales2Err_[ibin] = 2*TMath::Abs(h_scales->GetBinContent(ibin+1))*h_scales->GetBinError(ibin+1);
	masks_[ibin] = h_masks->GetBinContent(ibin+1);
	if( masks_[ibin]>0.5 ) n_unmasked_bins++;
      }
      
      n_dof_ = n_unmasked_bins - n_pars_;
      n_data_ = n_unmasked_bins;
      
      fin->Close();
    }
    
    // generate initial set of data points;
    //generate_data();

    if(bias>=0)
      n_dof_ = n_data_ - n_pars_;

    U_ = MatrixXd(n_pars_,n_pars_);
    for(unsigned int i=0; i<n_pars_; i++){
      for(unsigned int j=0; j<n_pars_; j++){
	// block(A,A)
	if(i<n_eta_bins_ && j<n_eta_bins_)
	  U_(i,j) = i==j ? 1.0 : 0.0;
	// block(A,e)
	else if(i<n_eta_bins_ && (j>=n_eta_bins_ && j<2*n_eta_bins_) )
	  U_(i,j) = i==(j-n_eta_bins_) ? kmean_val_ : 0.0;
	// block(e,e)
	else if(i>=n_eta_bins_ && i<2*n_eta_bins_ && j>=n_eta_bins_ && j<2*n_eta_bins_)
	  U_(i,j) = i==j ? kmean_val_ : 0.0;
	// block(M,M)
	else if(i>=2*n_eta_bins_  && j>=2*n_eta_bins_)
	  U_(i,j) = i==j ? 1.0/kmean_val_ : 0.0;
	else U_(i,j) = 0.0;
      }
    }
    //cout << U_ << endl;
    
  }
  
  ~TheoryFcn() { delete ran_;}

  void generate_data();

  void set_seed(const int& seed){ ran_->SetSeed(seed);}

  double get_true_params(const unsigned int& i, const bool& external){
    if(external)
      return x_vals_(i);
    else
      return (U_*x_vals_)(i);
  }

  double get_A_prevfit(const unsigned int& i){
    return A_vals_prevfit_(i);
  }
  double get_e_prevfit(const unsigned int& i){
    return e_vals_prevfit_(i);
  }
  double get_M_prevfit(const unsigned int& i){
    return M_vals_prevfit_(i);
  }


  unsigned int get_n_params(){ return n_pars_;}
  unsigned int get_n_data(){ return n_data_;}
  unsigned int get_n_dof(){ return n_dof_;} 

  double get_U(const unsigned int& i, const unsigned int& j){
    return U_(i,j);
  }
  
  virtual double Up() const {return errorDef_;}
  virtual void SetErrorDef(double def) {errorDef_ = def;}

  virtual double operator()(const vector<double>&) const;
  virtual vector<double> Gradient(const vector<double>& ) const;
  virtual bool CheckGradient() const {return true;} 

private:

  vector<double> scales2_;
  vector<double> scales2Err_;
  vector<int> masks_;
  vector<float> pt_edges_;
  vector<double> k_edges_;
  vector<double> kmean_vals_;
  VectorXd A_vals_;
  VectorXd e_vals_;
  VectorXd M_vals_;
  VectorXd x_vals_;
  VectorXd A_vals_prevfit_;
  VectorXd e_vals_prevfit_;
  VectorXd M_vals_prevfit_;
  double kmean_val_;
  vector<float> eta_edges_;
  unsigned int n_pt_bins_;
  unsigned int n_eta_bins_;
  unsigned int n_data_;
  unsigned int n_pars_;
  unsigned int n_dof_;
  int debug_;
  int seed_;
  int bias_;
  double errorDef_;
  MatrixXd U_;
  TRandom3* ran_;
};

void TheoryFcn::generate_data(){
  //ran_->SetSeed(seed_);
  double chi2_start = 0.;
  unsigned int ibin = 0;
  for(unsigned int ieta_p = 0; ieta_p<n_eta_bins_; ieta_p++){
    for(unsigned int ipt_p = 0; ipt_p<n_pt_bins_; ipt_p++){
      double k_p = kmean_vals_[ipt_p];
      for(unsigned int ieta_m = 0; ieta_m<n_eta_bins_; ieta_m++){
	for(unsigned int ipt_m = 0; ipt_m<n_pt_bins_; ipt_m++){
	  double k_m = kmean_vals_[ipt_m];

	  // was 0.001
	  double ierr2_nom = 0.0001*(1+double(ieta_p)/n_eta_bins_)*(1+double(ieta_m)/n_eta_bins_);

	  //*(2-0.1*double(ipt_p)/n_pt_bins_)*(2-0.1*double(ipt_m)/n_pt_bins_);
	  double ierr2 = ran_->Gaus(ierr2_nom,  ierr2_nom*0.1);
	  while(ierr2<=0.){
	    ierr2 = ran_->Gaus(ierr2_nom,  ierr2_nom*0.1);
	  }
	  
	  double iscale2_bias =
	    (1.0 + A_vals_(ieta_p) + e_vals_(ieta_p)*k_p - M_vals_(ieta_p)/k_p)*
	    (1.0 + A_vals_(ieta_m) + e_vals_(ieta_m)*k_m + M_vals_(ieta_m)/k_m);
	  double iscale2 = ran_->Gaus(iscale2_bias, ierr2);

	  //if(ibin<3) cout << iscale2 << endl;
	  scales2_[ibin]    = iscale2 ;
	  scales2Err_[ibin] =  ierr2 ;
	  double dchi2 = (scales2_[ibin]-1.0)/scales2Err_[ibin];
	  //cout << dchi2*dchi2 << endl;
	  chi2_start += dchi2*dchi2 ;
	  ibin++;
	}
      }
    }
  }
  cout << "Inistial chi2 = " << chi2_start << " / " << n_data_ << " ndof has prob " << TMath::Prob(chi2_start, n_data_ ) <<  endl;
  return;
}

double TheoryFcn::operator()(const vector<double>& par) const {

  double val = 0.0;
  const unsigned int npars = par.size();

  unsigned int ibin = 0;
  for(unsigned int ieta_p = 0; ieta_p < n_eta_bins_; ieta_p++){
    double A_p = par[ieta_p];
    double e_p = par[ieta_p+n_eta_bins_];
    double M_p = par[ieta_p+2*n_eta_bins_]; 
    for(unsigned int ipt_p = 0; ipt_p < n_pt_bins_; ipt_p++){   
      double k_p = kmean_vals_[ipt_p];
      double p_term = (1.0 + A_p + e_p*(k_p-kmean_val_)/kmean_val_ - M_p/k_p*kmean_val_ );
      for(unsigned int ieta_m = 0; ieta_m < n_eta_bins_; ieta_m++){
	double A_m = par[ieta_m];
	double e_m = par[ieta_m+n_eta_bins_];
	double M_m = par[ieta_m+2*n_eta_bins_];
	for(unsigned int ipt_m = 0; ipt_m < n_pt_bins_; ipt_m++){	  
	  double k_m = kmean_vals_[ipt_m];
	  double m_term = (1.0 + A_m + e_m*(k_m-kmean_val_)/kmean_val_ + M_m/k_m*kmean_val_);
	  double ival = (scales2_[ibin] - p_term*m_term)/scales2Err_[ibin];
	  double ival2 = ival*ival;
	  if(masks_[ibin])
	    val += ival2;
	  ibin++;
	}
      }
    }
  }
  val /= n_dof_;

  val -= 1.0;
  
  //val = 0.0;
  //for(unsigned int ipar=0; ipar<par.size(); ipar++) val += (par[ipar]-0.001)*(par[ipar]-0.001);
  //cout << val << endl;
  
  return val;
}

vector<double> TheoryFcn::Gradient(const vector<double> &par ) const {

  //cout << "Using gradient" << endl; 
  vector<double> grad(par.size(), 0.0);

  for(unsigned int ipar = 0; ipar < par.size(); ipar++){
    unsigned int ieta     = ipar % n_eta_bins_;
    unsigned int par_type = ipar / n_eta_bins_;    
    //cout << "ipar " << ipar << ": " << ieta << ", " << par_type << endl;
    double grad_i = 0.0;    
    unsigned int ibin = 0;
    for(unsigned int ieta_p = 0; ieta_p < n_eta_bins_; ieta_p++){
      double A_p = par[ieta_p];
      double e_p = par[ieta_p+n_eta_bins_];
      double M_p = par[ieta_p+2*n_eta_bins_]; 
      for(unsigned int ipt_p = 0; ipt_p < n_pt_bins_; ipt_p++){   
	double k_p = kmean_vals_[ipt_p];
	double p_term = 0.;
	if(ieta_p != ieta) p_term = (1.0 + A_p + e_p*(k_p-kmean_val_)/kmean_val_ - M_p/k_p*kmean_val_);
	else{
	  if(par_type==0)       p_term = 1.0;
	  else if( par_type==1) p_term = (k_p-kmean_val_)/kmean_val_;
	  else                  p_term = -1./k_p*kmean_val_;
	}
	for(unsigned int ieta_m = 0; ieta_m < n_eta_bins_; ieta_m++){
	  double A_m = par[ieta_m];
	  double e_m = par[ieta_m+n_eta_bins_];
	  double M_m = par[ieta_m+2*n_eta_bins_];
	  for(unsigned int ipt_m = 0; ipt_m < n_pt_bins_; ipt_m++){	  
	    double k_m = kmean_vals_[ipt_m];
	    double m_term = 0.;
	    if(ieta_m != ieta) m_term = (1.0 + A_m + e_m*(k_m-kmean_val_)/kmean_val_ + M_m/k_m*kmean_val_);
	    else{
	      if(par_type==0)       m_term = 1.0;
	      else if( par_type==1) m_term = (k_m-kmean_val_)/kmean_val_;
	      else                  m_term = +1./k_m*kmean_val_;
	    }
	    double ival = -2*(scales2_[ibin] - (1.0 + A_p + e_p*(k_p-kmean_val_)/kmean_val_ - M_p/k_p*kmean_val_)*
			      (1.0 + A_m + e_m*(k_m-kmean_val_)/kmean_val_ + M_m/k_m*kmean_val_))
	      /scales2Err_[ibin]/scales2Err_[ibin];

	    double term = 0.0;
	    if(ieta_p==ieta || ieta_m==ieta){
	      if(ieta_p!=ieta_m) term = p_term*m_term;
	      else{
		if(par_type==0)
		  term = 1.0*(1.0 + A_m + e_m*(k_m-kmean_val_)/kmean_val_ + M_m/k_m*kmean_val_) +
		    (1.0 + A_p + e_p*(k_p-kmean_val_)/kmean_val_ - M_p/k_p*kmean_val_)*1.0;
		else if(par_type==1)
		  term = (k_p-kmean_val_)/kmean_val_ * (1.0 + A_m + e_m*(k_m-kmean_val_)/kmean_val_ + M_m/k_m*kmean_val_) +
		    (1.0 + A_p + e_p*(k_p-kmean_val_)/kmean_val_ - M_p/k_p*kmean_val_) * (k_m-kmean_val_)/kmean_val_;
		else
		  term = -1.0/k_p*kmean_val_ * (1.0 + A_m + e_m*(k_m-kmean_val_)/kmean_val_ + M_m/k_m*kmean_val_) +
		    (1.0 + A_p + e_p*(k_p-kmean_val_)/kmean_val_ - M_p/k_p*kmean_val_) * 1.0/k_m*kmean_val_;
	      }
	    }
	    //cout << "ival " << ival << "," << term << endl; 
	    double ig = ival*term;
	    ig /= n_dof_;
	    //cout << "ibin " << ibin << " += " << ig << endl;
	    if(masks_[ibin])
	      grad_i += ig;
	    ibin++;
	  }
	}
      }
    }
    //cout << "\t" << ipar << ": " << grad_i << endl;
    //grad_i = 2*(par[ipar]-0.001);
    grad[ipar] = grad_i;
  }

  return grad; 
}


int massfit()
//int main(int argc, char* argv[])
{

  TStopwatch sw;
  sw.Start();

  //ROOT::EnableImplicitMT();
  /*
  variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("nevents",     value<long>()->default_value(1000), "number of events")
	("lumi",        value<long>()->default_value(1000), "number of events")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("bias",        value<int>()->default_value(0), "bias")
	("infile",      value<std::string>()->default_value("massscales"), "run type")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  long lumi       = vm["lumi"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string infile = vm["infile"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int bias        = vm["bias"].as<int>();
  int seed        = vm["seed"].as<int>();
  */

  long nevents = 1, lumi = 1000;
  int bias = -1, seed = 4357;
  std::string tag("closure"), run("closure"), infile("massscales");
  
  TFile* fout = TFile::Open(("/home/users/alexe/workingarea/MomentumBiases/OtherFitter/massfit_"+tag+"_"+run+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "tree");
  
  double edm, fmin, prob;
  int isvalid, hasAccurateCovar, hasPosDefCovar;
  tree->Branch("edm", &edm, "edm/D");
  tree->Branch("fmin", &fmin, "fmin/D");
  tree->Branch("prob", &prob, "prob/D");
  tree->Branch("isvalid", &isvalid, "isvalid/I");
  tree->Branch("hasAccurateCovar", &hasAccurateCovar, "hasAccurateCovar/I");
  tree->Branch("hasPosDefCovar", &hasPosDefCovar, "hasPosDefCovar/I");
  //fFCN->set_seed(seed);

  int debug = 0;
  string infname = "InOutputFiles/mass_fits_control_histos_smear_beta_val.root";
  TheoryFcn* fFCN = new TheoryFcn(debug, seed, bias, infname);  
  fFCN->SetErrorDef(1.0 / fFCN->get_n_dof());
  unsigned int n_parameters = fFCN->get_n_params();
  MatrixXd U(n_parameters,n_parameters);
  for (int i=0; i<n_parameters; i++){
    for (int j=0; j<n_parameters; j++){
      U(i,j) = fFCN->get_U(i,j);
    }
  }
  MatrixXd Uinv = U.inverse();
  
  vector<double> tparIn0(n_parameters);
  vector<double> tparIn(n_parameters);
  vector<double> tparInErr(n_parameters);
  vector<double> tparOut0(n_parameters);
  vector<double> tparOut(n_parameters);
  vector<double> tparOutErr(n_parameters);

  for (int i=0; i<n_parameters/3; i++){
    tree->Branch(Form("A%d",i),       &tparOut[i],    Form("A%d/D",i));
    tree->Branch(Form("A%d_true",i),  &tparOut0[i],   Form("A%d_true/D",i));
    tree->Branch(Form("A%d_err",i),   &tparOutErr[i], Form("A%d_err/D",i));
    tree->Branch(Form("A%d_in",i),    &tparIn[i],    Form("A%d_in/D",i));
    tree->Branch(Form("A%d_intrue",i),&tparIn0[i],    Form("A%d_intrue/D",i));
    tree->Branch(Form("A%d_inerr",i), &tparInErr[i], Form("A%d_inerr/D",i));
  }
  for (int i=0; i<n_parameters/3; i++){
    tree->Branch(Form("e%d",i),        &tparOut[i+n_parameters/3],    Form("e%d/D",i));
    tree->Branch(Form("e%d_true",i),   &tparOut0[i+n_parameters/3],   Form("e%d_true/D",i));
    tree->Branch(Form("e%d_err",i),    &tparOutErr[i+n_parameters/3], Form("e%d_err/D",i));
    tree->Branch(Form("e%d_in",i),     &tparIn[i+n_parameters/3],     Form("e%d_in/D",i));
    tree->Branch(Form("e%d_intrue",i), &tparIn0[i+n_parameters/3],    Form("e%d_intrue/D",i));
    tree->Branch(Form("e%d_inerr",i),  &tparInErr[i+n_parameters/3],  Form("e%d_inerr/D",i));
  }
  for (int i=0; i<n_parameters/3; i++){
    tree->Branch(Form("M%d",i),        &tparOut[i+2*n_parameters/3],    Form("M%d/D",i));
    tree->Branch(Form("M%d_true",i),   &tparOut0[i+2*n_parameters/3],   Form("M%d_true/D",i));
    tree->Branch(Form("M%d_err",i),    &tparOutErr[i+2*n_parameters/3], Form("M%d_err/D",i));
    tree->Branch(Form("M%d_in",i),     &tparIn[i+2*n_parameters/3],     Form("M%d_in/D",i));
    tree->Branch(Form("M%d_intrue",i), &tparIn0[i+2*n_parameters/3],    Form("M%d_in/D",i));
    tree->Branch(Form("M%d_inerr",i),  &tparInErr[i+2*n_parameters/3],  Form("M%d_inerr/D",i));
  }

  TH1D* h_A_vals_nom  = new TH1D("h_A_vals_nom", "A nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_e_vals_nom  = new TH1D("h_e_vals_nom", "e nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_M_vals_nom  = new TH1D("h_M_vals_nom", "M nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Ain_vals_nom  = new TH1D("h_Ain_vals_nom", "(A+e#bar{k})", n_parameters/3, 0, n_parameters/3);
  TH1D* h_ein_vals_nom  = new TH1D("h_ein_vals_nom", "e/#bar{k} nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Min_vals_nom  = new TH1D("h_Min_vals_nom", "M#bar{k} nominal", n_parameters/3, 0, n_parameters/3);

  TH1D* h_A_vals_fit  = new TH1D("h_A_vals_fit", "#hat{A}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_e_vals_fit  = new TH1D("h_e_vals_fit", "#hat{e}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_M_vals_fit  = new TH1D("h_M_vals_fit", "#hat{M}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Ain_vals_fit  = new TH1D("h_Ain_vals_fit", "(#hat{A}+#hat{e}#bar{k})", n_parameters/3, 0, n_parameters/3);
  TH1D* h_ein_vals_fit  = new TH1D("h_ein_vals_fit", "#hat{e}/#bar{k}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Min_vals_fit  = new TH1D("h_Min_vals_fit", "#hat{M}#bar{k}", n_parameters/3, 0, n_parameters/3);
  
  TH1D* h_A_vals_prevfit  = new TH1D("h_A_vals_prevfit", "#hat{A}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_e_vals_prevfit  = new TH1D("h_e_vals_prevfit", "#hat{e}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_M_vals_prevfit  = new TH1D("h_M_vals_prevfit", "#hat{M}", n_parameters/3, 0, n_parameters/3);

  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  double tolerance(0.001);
  int verbosity = int(nevents<2); 
  ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);
  
  for(unsigned int itoy=0; itoy<nevents; itoy++){

    if(itoy%10==0) cout << "Toy " << itoy << " / " << nevents << endl;

    if(bias>=0)
      fFCN->generate_data();
    
    MnUserParameters upar;
    double start=0.0, par_error=0.01;
    for (int i=0; i<n_parameters/3; i++){
      upar.Add(Form("A%d",i), start, par_error);
    }
    for (int i=0; i<n_parameters/3; i++){
      upar.Add(Form("e%d",i), start, par_error);
    }
    for (int i=0; i<n_parameters/3; i++){
      upar.Add(Form("M%d",i), start, par_error);      
    }

    MnMigrad migrad(*fFCN, upar, 1);    

    //fFCN->set_seed(seed);

    cout << "\tMigrad..." << endl;
    FunctionMinimum min = migrad(maxfcn, tolerance);

    edm = double(min.Edm());
    fmin = double(min.Fval());
    prob = TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() );
    isvalid = int(min.IsValid());
    hasAccurateCovar = int(min.HasAccurateCovar());
    hasPosDefCovar = int(min.HasPosDefCovar());

    cout << "\tHesse..." << endl;
    MnHesse hesse(1);
    hesse(*fFCN, min);

    cout << "\t => final chi2/ndf: " << min.Fval()+1 << " (prob: " << TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() ) << ")" <<  endl;;
    
    MatrixXd Vin(n_parameters,n_parameters);    
    for(unsigned int i = 0 ; i<n_parameters; i++){    
      for(unsigned int j = 0 ; j<n_parameters; j++){
	Vin(i,j) = i>j ?
	  min.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    MatrixXd Vout = Uinv*Vin*Uinv.transpose();

    VectorXd xin(n_parameters);
    VectorXd xinErr(n_parameters);     
    for(unsigned int i = 0 ; i<n_parameters; i++){
      xin(i)    = min.UserState().Value(i) ;
      xinErr(i) = min.UserState().Error(i) ;
    }
    VectorXd x = Uinv*xin;
    VectorXd xErr(n_parameters);
    for(unsigned int i = 0 ; i<n_parameters; i++){
      xErr(i) = TMath::Sqrt(Vout(i,i));
    }    
    for(unsigned int i = 0 ; i<n_parameters; i++){
      tparIn[i]     = xin(i);
      tparIn0[i]    = fFCN->get_true_params(i, false) ;
      tparInErr[i]  = xinErr(i);
      tparOut[i]    = x(i);
      tparOut0[i]   = fFCN->get_true_params(i, true) ;
      tparOutErr[i] = xErr(i);
      int ip = i%(n_parameters/3);
      if(i<n_parameters/3){
	h_A_vals_fit->SetBinContent(ip+1, x(i));
	h_A_vals_fit->SetBinError(ip+1, xErr(i));
	h_A_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, true));
	h_Ain_vals_fit->SetBinContent(ip+1, xin(i));
	h_Ain_vals_fit->SetBinError(ip+1, xinErr(i));
	h_Ain_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, false));
	h_A_vals_prevfit->SetBinContent(ip+1, fFCN->get_A_prevfit(ip) + x(i));
      }
      else if(i>=n_parameters/3 && i<2*n_parameters/3){
	h_e_vals_fit->SetBinContent(ip+1, x(i));
	h_e_vals_fit->SetBinError(ip+1, xErr(i));
	h_e_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, true));
	h_ein_vals_fit->SetBinContent(ip+1, xin(i));
	h_ein_vals_fit->SetBinError(ip+1, xinErr(i));
	h_ein_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, false));
	h_e_vals_prevfit->SetBinContent(ip+1, fFCN->get_e_prevfit(ip) + x(i));
      }
      else{
	h_M_vals_fit->SetBinContent(ip+1, x(i));
	h_M_vals_fit->SetBinError(ip+1, xErr(i));
	h_M_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, true));
	h_Min_vals_fit->SetBinContent(ip+1, xin(i));
	h_Min_vals_fit->SetBinError(ip+1, xinErr(i));
	h_Min_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, false));
	h_M_vals_prevfit->SetBinContent(ip+1, fFCN->get_M_prevfit(ip) + x(i));
      }
      //cout << "Param " << i << ": " << x(i) << " +/- " << xErr(i) << ". True value is " << fFCN->get_true_params(i, true) << endl;
    }

    TH2D* hcov = new TH2D(Form("hcov_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);
    TH2D* hcor = new TH2D(Form("hcor_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);  
    TH2D* hcovin = new TH2D(Form("hcovin_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);
    TH2D* hcorin = new TH2D(Form("hcorin_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);  

    for(unsigned int i = 0 ; i<n_parameters; i++){    
      hcov->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
      hcor->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
      hcovin->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
      hcorin->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
      for(unsigned int j = 0 ; j<n_parameters; j++){
	double covin_ij = Vin(i,j);
	double corin_ij = Vin(i,j)/TMath::Sqrt(Vin(i,i)*Vin(j,j)); 
	double cov_ij = Vout(i,j);
	double cor_ij = Vout(i,j)/TMath::Sqrt(Vout(i,i)*Vout(j,j)); 
	hcov->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	hcor->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	hcovin->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	hcorin->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	hcovin->SetBinContent(i+1, j+1, covin_ij);
	hcorin->SetBinContent(i+1, j+1, corin_ij);
	hcov->SetBinContent(i+1, j+1, cov_ij);
	hcor->SetBinContent(i+1, j+1, cor_ij);
      }
    }

    tree->Fill();

    hcor->SetMinimum(-1.0);
    hcor->SetMaximum(+1.0);
    hcorin->SetMinimum(-1.0);
    hcorin->SetMaximum(+1.0);
    fout->cd();
    if(itoy<1){
      hcor->Write();
      hcov->Write();
      hcorin->Write();
      hcovin->Write();
    }

    if(verbosity){
      cout << "Data points: " << fFCN->get_n_data() << endl;
      cout << "Number of parameters: " << fFCN->get_n_params() << endl;
      cout << "chi2/ndf: " << min.Fval()+1 << " (prob: " << TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() ) << ")" <<  endl;;
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
    }
  }

  fout->cd();
  tree->Write();

  TH1D* hpulls = new TH1D("hpulls", "", n_parameters, 0, n_parameters);
  TH1D* hsigma = new TH1D("hsigma", "", n_parameters, 0, n_parameters);
  for (int i=0; i<n_parameters; i++){
    TH1D* h = new TH1D(Form("h%d", i), "", 100,-3,3);
    int ip = i%(n_parameters/3);
    if(i<n_parameters/3){
      tree->Draw(Form("(A%d - A%d_true)/A%d_err>>h%d", ip, ip, ip, i), "", "");
      hpulls->GetXaxis()->SetBinLabel(i+1, Form("A%d", ip));
      h_A_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("A%d", ip));
      h_Ain_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("Ain%d", ip));
      h_A_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("A%d", ip));
      h_Ain_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("Ain%d", ip));
    }
    else if(i>=n_parameters/3 && i<2*n_parameters/3){
      tree->Draw(Form("(e%d - e%d_true)/e%d_err>>h%d", ip, ip, ip, i), "", "");
      hpulls->GetXaxis()->SetBinLabel(i+1, Form("e%d", ip));
      h_e_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("e%d", ip));
      h_ein_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("ein%d", ip));
      h_e_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("e%d", ip));
      h_ein_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("ein%d", ip));
    }
    else{
      tree->Draw(Form("(M%d - M%d_true)/M%d_err>>h%d", ip, ip, ip, i), "", "");
      hpulls->GetXaxis()->SetBinLabel(i+1, Form("M%d", ip));
      h_M_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("M%d", ip));
      h_Min_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("Min%d", ip));
      h_M_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("M%d", ip));
      h_Min_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("Min%d", ip));
    }
    //cout << i << "-->" << h->GetMean() << endl;
    float pull_i = h->GetMean();
    float pull_i_err = 0.;
    float sigma_i = 0.;
    float sigma_i_err = 0.;
    if(h->GetEntries()>10){
      h->Fit("gaus", "Q");
      TF1* gaus = (TF1*)h->GetFunction("gaus");
      if(gaus==0){
	cout << "no func" << endl;
	continue;
      }
      pull_i = gaus->GetParameter(1);
      pull_i_err = gaus->GetParError(1);
      sigma_i = gaus->GetParameter(2);
      sigma_i_err = gaus->GetParError(2);
    }
    hpulls->SetBinContent(i+1, pull_i);
    hpulls->SetBinError(i+1, pull_i_err);
    hsigma->SetBinContent(i+1, sigma_i);
    hsigma->SetBinError(i+1, sigma_i_err);
    delete h;
  }
  hpulls->Write();
  hsigma->Write();

  h_A_vals_fit->Write();
  h_e_vals_fit->Write();
  h_M_vals_fit->Write();
  h_A_vals_prevfit->Write();
  h_e_vals_prevfit->Write();
  h_M_vals_prevfit->Write();
  h_Ain_vals_fit->Write();
  h_ein_vals_fit->Write();
  h_Min_vals_fit->Write();
  h_A_vals_nom->Write();
  h_e_vals_nom->Write();
  h_M_vals_nom->Write();
  h_Ain_vals_nom->Write();
  h_ein_vals_nom->Write();
  h_Min_vals_nom->Write();
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 
  //for(auto r : rans) delete r;

  return 1;
}
