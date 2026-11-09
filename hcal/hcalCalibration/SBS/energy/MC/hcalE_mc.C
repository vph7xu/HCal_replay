//SSeeds 3.2.23 GMn: Script to evaluate HCal sampling fraction and energy spectra by kinematic from MC data with primary cluster constraint.
// Update 5.29.23: added class structure functionality. mc_data_path will need updated for other users.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TF1.h"
#include "../../include/sbs.h"

//#include <bits/stdc++.h>

//Main.
void hcalE_mc( const char *experiment = "gen", Int_t config=3 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  std::string date = util::getDate();

  TChain *C = new TChain("T");

  //set up paths and variables, first set from simulations 032123, second set from 091523
  //std::string outdir_path = gSystem->Getenv("OUT_DIR");
  //std::string mc_data_path = Form("/lustre19/expphy/volatile/halla/sbs/seeds/sim091523/%s_sbs%d*",experiment,config);
  std::string mc_data_path = Form("/lustre24/expphy/volatile/halla/sbs/vimukthi/MC/GEN%d_H2/GEN%d_H2_elastic*",config,config);
  //std::string mc_fout_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_mc_class_%s_conf%d.root",experiment,config);
  std::string mc_fout_path = Form("./hcalE_mc_class_%s_conf%d.root",experiment,config);

  // Declare outfile
  TFile *fout = new TFile( mc_fout_path.c_str(), "RECREATE" );

  // Declare chain for potentially many root files
  C->Add( mc_data_path.c_str() );

  // Get beam energy by configuration
  SBSconfig config_parameters(experiment,config);    

  Double_t beamE = config_parameters.GetEbeam();

  std::cout << endl << "Setup parameters loaded and chains linked." << std::endl;

  // Set up MC structure
  MC *G = new MC(C);

  // Get total number of events
  Int_t Nevents = C->GetEntries();

  //Declare diagnostic histograms (as sparse as possible)
  TH1D *hSampFrac_det = new TH1D( "hSampFrac_det","HCal Det Sum E Dep / (Beam E - P Track E)", 400, 0., 1. );
  TH1D *hSampFrac_clus = new TH1D( "hSampFrac_clus","HCal Cluster E Dep / (Beam E - P Track E)", 400, 0., 1. );
  TH1D *hep = new TH1D( "hep", "Event Electron Momentum", 500, 0., 5. );
  TH2D *hcEvID = new TH2D( "hcEvID", "HCal Cluster Element E vs Cell", 288, 0, 288, 400, 0., 1.);
  TH1D *csize = new TH1D( "csize", "P Cluster size, hit thresh 0.01 GeV", 30, 0, 30 );
  TH1D *bhits = new TH1D( "bhits", "Number of blocks over 1MeV thresh", 30, 0, 30 );
  TH1D *hHCALe = new TH1D( "hHCALe","HCal Detector E", 400, 0., 1. );
  TH1D *hHCALe_clus = new TH1D( "hHCALe_clus","HCal Cluster E", 400, 0., 1. );
  TH2D *hcEvID_1b = new TH2D( "hcEvID_1b", "HCal One Block Cluster E vs Cell", 288, 0, 288, 400, 0., 10.);
  TH1D *hNKE = new TH1D( "hNKE","Nucleon KE", 300, 0., 15. );

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;

  for( Long64_t nevent = 1; nevent <Nevents; nevent++){
      
    if ( nevent%1000==0 ) cout << "MC configuration " << config << ", entry: " << nevent << "/" << Nevents << " \r";
      cout.flush();

      G->GetEntry( nevent ); 
      ///////////////////////////////////////////////////
      //General cuts
      Double_t xptrack = (*(G->PTrack_momx))[0]; //primary track momentum
      if( xptrack<0 ) continue;
      //cout << xptrack << endl;

      Int_t nhits = G->Harm_HCalScint_hit_nhits;
      if( nhits<1 ) continue; //Cut on hits in hcal
      //cout << nhits << endl;

      Double_t PSE = G->Earm_BBPSTF1_det_esum;
      if( PSE<0.2 ) continue; //Cut on BBCal preshower energy
      //cout << PSE << endl;

      ///////////////////////////////////////////////////

      Double_t hcalE = G->Harm_HCalScint_det_esum;
      hHCALe->Fill(hcalE);
      Double_t ep = G->ev_ep;

      Double_t detSF = hcalE/(beamE-ep);
      hSampFrac_det->Fill(detSF);

      Int_t blockhits = 0;
      Double_t clustermax = 0.;
      Int_t highblock = -1;
      Int_t highblockR = -1;
      Int_t highblockC = -1;
      for( Int_t h=0; h<nhits; h++ ){
	Double_t blockE = (*(G->Harm_HCalScint_hit_sumedep))[h];
	if( blockE < 0.01 ) continue; //Min E to be considered for a cluster in SBS-offline
	Int_t cell = (*(G->Harm_HCalScint_hit_cell))[h];
	Int_t row = (*(G->Harm_HCalScint_hit_row))[h];
	Int_t col = (*(G->Harm_HCalScint_hit_col))[h];

	blockhits++;
	if( blockE>clustermax ){
	  clustermax=blockE;
	  highblock=cell;
	  highblockR=row;
	  highblockC=col;
	}
      }
      if( highblock==-1 ) continue; //Cut event if no cluster seed of sufficient E is found

      if( highblockR==0 || 
	  highblockR==23 || 
	  highblockC==0 || 
	  highblockC==11 ) continue; //All events with primary cluster element on edge blocks cut
      
      //Make check on block position to build primary cluster
      Double_t clusterE = 0.;
      Int_t clustersize = 0;
      for( Int_t h=0; h<nhits; h++ ){
	Double_t blockE = (*(G->Harm_HCalScint_hit_sumedep))[h];
	if( blockE < 0.01 ) continue; //Min E to be considered for a cluster in SBS-offline
	Int_t cell = (*(G->Harm_HCalScint_hit_cell))[h];
	Int_t row = (*(G->Harm_HCalScint_hit_row))[h];
	Int_t col = (*(G->Harm_HCalScint_hit_col))[h];

	if( abs(highblockR-row)>3 ) continue; //Cluster row cut
	if( abs(highblockC-col)>3 ) continue; //Cluster col cut
	clustersize++;
	clusterE+=blockE;
	hcEvID->Fill(cell,blockE);
      }
      if( nhits==1 ) hcEvID_1b->Fill( (*(G->Harm_HCalScint_hit_cell))[0], clusterE );

      if( clusterE<0.01 ) continue; //Mirror min cluster energy required to be considered in calibration

      hHCALe_clus->Fill(clusterE);
      hSampFrac_clus->Fill(clusterE/(beamE-ep));
      hep->Fill(ep);
      csize->Fill(clustersize);
      bhits->Fill(blockhits);

  }

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}

//int MID = (*(G->SDTrack_MID))[(*(G->Harm_HCalScint_hit_sdtridx))[ihit]];
