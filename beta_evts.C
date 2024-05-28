using namespace std;

int beta_evts(){
  
  std::unique_ptr<TFile> inputFile( TFile::Open("mass_fits_control_histos_smear_beta_val.root") );
  std::unique_ptr<TH1D> beta(inputFile->Get<TH1D>("beta"));
  std::unique_ptr<TH1D> evts(inputFile->Get<TH1D>("bin_occupancy")); 

  int n_data_points, n_data_points_evts, no_evts;
  double beta_sig;
  string label;

  n_data_points = beta->GetEntries();
  n_data_points_evts = evts->GetEntries();
  
  cout << "beta "<<n_data_points<<"evts "<<n_data_points_evts;
  if(n_data_points == n_data_points_evts){
  
    TH2D *beta_evts = new TH2D("beta_evts", "#beta/error vs events in k bin", 30, 0.0, 10.0, 50, evts->GetBinContent(evts->GetMinimumBin()), evts->GetBinContent(evts->GetMaximumBin()));
    beta_evts->GetXaxis()->SetTitle("#beta/error");
    beta_evts->GetYaxis()->SetTitle("evts");
    
    for(int i=0; i<n_data_points; i++){
      beta_sig = abs(beta->GetBinContent(i+1)) / beta->GetBinError(i+1);
      no_evts = evts->GetBinContent(i+1);
      beta_evts->Fill(beta_sig, no_evts);
    }

    TFile f("beta_evts.root","recreate");

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    beta_evts->SetStats(0);
    beta_evts->Draw("COLZ");

    f.WriteObject(c1, "beta_evts");
    f.Close();
  }
  return 0; 
  
}
