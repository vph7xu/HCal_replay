//SSeeds 5.15.23 - Script adapted from ecal.C at https://github.com/sebastianseeds/HCal_replay to calibrate energy with elastic events from all available data. Assumes tdc/adct alignment scripts have been written and new alignments are available. Check sbs.h for class details. Updated to include rootfile directory and hcal vertical offset in runlist csv
//NOTE: This script is configured to calibrate by configuration/kinematic. One could, in principle, calibrate the energy over all configurations. This method is chosen to be sensitive to expected changes in hcal sampling fraction and E/sigma ratios prediced by MC. Timestamps are dynamically handled to address mismatches. Any timestamps occurring within a configuration are assumed to be hardware/HV related and will result in multiple calibration sets for that configuration.
//ADDITIONAL NOTE: requires $DB_DIR and $OUT_DIR paths set correctly in current environment and DB directory structure to match the github model from SBS-replay.

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

const Double_t minEventPerCell = 10; //minimum elastic events detected in cell to be considered for calibration
const Double_t highDelta = 0.01; //max allowed discrepancy factor between energy deposited in block and expected
const Int_t hcal_first_channel = 0; //first channel in hcal index
const Int_t bins_magnet_percent = 21; //SBS field at 2100A, so divide by this factor to get percent
const Double_t start_magnet_percent = 0.;
const Double_t end_magnet_percent = 105.;
const Int_t bins_dxdy = 250;
const Double_t start_dx = -4.;
const Double_t end_dx = 3.;
const Double_t start_dy = -1.25;
const Double_t end_dy = 1.25;
const Double_t bins_W2 = 1000;
const Double_t start_W2 = -2.0;
const Double_t end_W2 = 3.0;
const Double_t bins_SFE = 400;
const Double_t start_SFE = 0.;
const Double_t end_SFE = 1.4;
const Double_t bins_coin = 200;
const Double_t start_coin = -50.;
const Double_t end_coin = 50.;
const Int_t linecount = 25;
const Int_t atimeNsig = 2;

bool altpath = false; //Added for manual paths

//Main <experiment> <configuration> <quasi-replay-option> <replay-pass> <target-option>; qreplay should only be performed after new offsets obtained
void ecal( const char *experiment = "gen", Int_t config = 3, bool qreplay = true, Int_t pass = 2, bool h2only = false, bool new_atime = false, bool verbose = true ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();

  // outfile path
  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";
  std::string outdir_path = "/work/halla/sbs/vimukthi/HCal_replay/hcal";//gSystem->Getenv("OUT_DIR");
  std::string ecal_path = outdir_path + Form("/hcalCalibration/SBS/energy/ecal%s_%s_conf%d_qr%d_pass%d.root",h2opt.c_str(),experiment,config,(Int_t)qreplay,pass);

  if(altpath)
    ecal_path = outdir_path + Form("/hcalCalibration/SBS/energy/ecal_newclus_%s_%s_conf%d_qr%d_pass%d.root",h2opt.c_str(),experiment,config,(Int_t)qreplay,pass);

  //Set up path variables and output files
  string db_path = gSystem->Getenv("DB_DIR");
  TFile *fout = new TFile( ecal_path.c_str(), "RECREATE" );
  std::string new_adcgain_path = Form("parameters/adcgaincoeff_%s%s_conf%d_pass%d.txt",experiment,h2opt.c_str(),config,pass);

  if(altpath)
    new_adcgain_path = Form("parameters/adcgaincoeff_newclus%s%s_conf%d_pass%d.txt",experiment,h2opt.c_str(),config,pass);

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s_p%d/",experiment,pass); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  cout<<"debug : 00000"<<endl;
  
  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
 
  cout<<"debug : 000"<<endl;

  //define structures for holding calibration data and indices
  calset gain_coeff[hcal::gNstamp];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;

  //define structures for holding cut report data and indices	    
  reportset report_set[hcal::gNmag];
  Int_t report_set_size = 0;
  Int_t report_set_idx = 0;
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree");   
  
  // Timing
  Double_t pblkid_out; //hcal cluster, primary block id
  Double_t tdc_out; //hcal cluster tdc time, tree
  Double_t atime_out; //hcal cluster adc time, tree
  Double_t cointime_out;
  Double_t hodotmean_out; //hodoscope cluster mean tdc time
  Double_t pblkcol_out;
  Double_t pblkrow_out;

  // Physics
  Double_t dx_out; //hcal cluster dx
  Double_t dy_out; //hcal cluster dy
  Double_t HCALx_out;
  Double_t HCALy_out;
  Double_t W2_out; //Invariant mass squared W2
  Double_t Q2_out; //Inverse momentum squared Q2
  Double_t hcale_out; //hcal cluster energy
  Double_t sfrac_out; //hcal cluster energy
  Double_t thetapq_pout; //proton thetapq
  Double_t thetapq_nout; //neutron thetapq
  Double_t nu_out; //hadron KE
  Double_t ep_out; //track reconstructed e' momentum
  Int_t nblk_out; //total blocks in primary cluster
  Int_t pid_out; //particle id (-1:neither,0:ambiguous,1:proton,2:neutron)
  Int_t mag_out; //sbs magnetic field strength (percent)
  Int_t run_out; //run number
  Int_t tar_out; //target, LH2 or LD2
  Int_t calset_out; //calibration set in chronological order
  Int_t badclus_out; //0: all blocks pass coin, 1: at least one block fails coin
  Int_t failedglobal_out;

  //cluster block variables
  Double_t cblkid_out[hcal::maxHCalBlk];
  Double_t cblkatime_out[hcal::maxHCalBlk];
  Double_t cblktime_out[hcal::maxHCalBlk];
  Double_t cblke_out[hcal::maxHCalBlk];

  // Set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "tdc", &tdc_out, "tdc/D" );
  P->Branch( "atime", &atime_out, "atime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "cointime", &cointime_out, "cointime/D");
  P->Branch( "pblkcol", &pblkcol_out, "pblkcol/D");
  P->Branch( "pblkrow", &pblkrow_out, "pblkrow/D");
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "sfrac", &sfrac_out, "sfrac/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_p/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_n/D" );
  P->Branch( "nblk", &nblk_out, "nblk_out/I" );
  P->Branch( "pid", &pid_out, "pid_out/I" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );
  P->Branch( "calset", &calset_out, "calset_out/I" );
  P->Branch( "badclus", &badclus_out, "badclus_out/I" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal_out/I" );

  P->Branch( "cblkid", &cblkid_out, Form("cblkid[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblkatime", &cblkatime_out, Form("cblkatime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblktime", &cblktime_out, Form("cblktime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblke", &cblke_out, Form("cblke[%d]/D",hcal::maxHCalBlk) );

  P->Branch("HCALx", &HCALx_out,"HCALx/D");
  P->Branch("HCALy", &HCALy_out,"HCALy/D");


  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();
  // Double_t hcal_v_offset;
  // if( pass>1 )
  //   hcal_v_offset = hcal::HCalvoff_p2;
  // else
  //   hcal_v_offset = hcal::HCalvoff;

  //Add quality plots
  //Declare histograms for modification later where additional calibration sets are necessary
  //troubleshooting vimukthi
  TH1D *h_cointime = new TH1D("h_cointime","cointime",2500,-50,200);


  TH1D *hcoin[hcal::gNstamp];
  TH1D *hW2[hcal::gNstamp];

  TH1D *hE[hcal::gNstamp];
  TH1D *hEblk[hcal::gNstamp];
  TH1D *hSF[hcal::gNstamp];

  TH2D *hdxdy[hcal::gNstamp];
  TH2D *hXY[hcal::gNstamp];
  TH2D *hdxvmag_h[hcal::gNstamp];
  TH2D *hdyvmag_h[hcal::gNstamp];
  TH2D *hdxvmag_d[hcal::gNstamp];
  TH2D *hdyvmag_d[hcal::gNstamp];
  TH2D *hdxvmag_h3[hcal::gNstamp];
  TH2D *hdyvmag_h3[hcal::gNstamp];
  TH2D *hEvX[hcal::gNstamp];
  TH2D *hEvY[hcal::gNstamp];
  TH2D *hSFvID[hcal::gNstamp];
  TH2D *hSFvrow[hcal::gNstamp];
  TH2D *hSFvcol[hcal::gNstamp];

  for( Int_t set=0; set<hcal::gNstamp; set++ ){
    
    hE[set] = new TH1D(Form("hE_set%d",set),
		       "null",
		       bins_SFE,
		       start_SFE,
		       end_SFE);
    
    hcoin[set] = new TH1D(Form("hcoin_set%d",set),
		       "null",
		       bins_coin,
		       start_coin,
		       end_coin);

    hW2[set] = new TH1D(Form("hW2_set%d",set),
		    "null",
		    bins_W2,
		    start_W2,
		    end_W2);

    hEblk[set] = new TH1D(Form("hEblk_set%d",set),
			  "null",
			  bins_SFE,
			  start_SFE,
			  end_SFE);
    
    hSF[set] = new TH1D(Form("hSF_set%d",set),
			"null",
			bins_SFE,
			start_SFE,
			end_SFE);

    hdxdy[set] = new TH2D(Form("hdxdy%d",set),
		    "null",
		    200,
		    -4,
		    4,
		    200,
		    -4,
		    4);

    hXY[set] = new TH2D(Form("hXY_set%d",set),
			      "null",
			      600,
			      -3,
			      3,
			      600,
			      -3,
			      3);

    hdxvmag_h[set] = new TH2D(Form("hdxvmag_h_set%d",set),
			      "null",
			      bins_magnet_percent,
			      start_magnet_percent,
			      end_magnet_percent,
			      bins_dxdy,
			      start_dx,
			      end_dx);
    
    hdyvmag_h[set] = new TH2D(Form("hdyvmag_h_set%d",set),
			      "null",
			      bins_magnet_percent,
			      start_magnet_percent,
			      end_magnet_percent,
			      bins_dxdy,
			      start_dy,
			      end_dy);
    
    hdxvmag_d[set] = new TH2D(Form("hdxvmag_d_set%d",set),
			      "null",
			      bins_magnet_percent,
			      start_magnet_percent,
			      end_magnet_percent,
			      bins_dxdy,
			      start_dx,
			      end_dx);
    
    hdyvmag_d[set] = new TH2D(Form("hdyvmag_d_set%d",set),
			      "null",
			      bins_magnet_percent,
			      start_magnet_percent,
			      end_magnet_percent,
			      bins_dxdy,
			      start_dy,
			      end_dy);
    
    hdxvmag_h3[set] = new TH2D(Form("hdxvmag_h3_set%d",set),
			       "null",
			       bins_magnet_percent,
			       start_magnet_percent,
			       end_magnet_percent,
			       bins_dxdy,
			       start_dx,
			       end_dx);
    
    hdyvmag_h3[set] = new TH2D(Form("hdyvmag_h3_set%d",set),
			       "null",
			       bins_magnet_percent,
			       start_magnet_percent,
			       end_magnet_percent,
			       bins_dxdy,
			       start_dy,
			       end_dy);
   
    hEvX[set] = new TH2D(Form("hEvX_set%d",set),
			 "null",
			 hcal::maxHCalRows,
			 hcal::posHCalXi,
			 hcal::posHCalXf,
			 bins_SFE,
			 start_SFE,
			 end_SFE);
   
    hEvY[set] = new TH2D(Form("hEvY_set%d",set),
			 "null",
			 hcal::maxHCalRows,
			 hcal::posHCalYi,
			 hcal::posHCalYf,
			 bins_SFE,
			 start_SFE,
			 end_SFE);

    hSFvID[set] = new TH2D(Form("hSFvID_set%d",set),
			   "null",
			   hcal::maxHCalChan,
			   hcal_first_channel,
			   hcal::maxHCalChan,
			   bins_SFE,
			   start_SFE,
			   end_SFE);

    hSFvrow[set] = new TH2D(Form("hSFvrow_set%d",set),
			   "null",
			   hcal::maxHCalRows,
			   hcal_first_channel,
			   hcal::maxHCalRows,
			   bins_SFE,
			   start_SFE,
			   end_SFE);

    hSFvcol[set] = new TH2D(Form("hSFvcol_set%d",set),
			   "null",
			   hcal::maxHCalCols,
			   hcal_first_channel,
			   hcal::maxHCalCols,
			   bins_SFE,
			   start_SFE,
			   end_SFE);

  }

  //reporting indices
  Int_t target_change_index;
  Double_t config_sampling_fraction;
  Double_t config_e_sigma_ratio;
  bool first = true;
  Double_t test_adcg_coeff[hcal::maxHCalChan] = {1.};
      
  long postcut1_nevent[hcal::gNstamp] = {0}, postcut2_nevent[hcal::gNstamp] = {0}, postcut3_nevent[hcal::gNstamp] = {0}, postcut4_nevent[hcal::gNstamp] = {0}, postcut5_nevent[hcal::gNstamp] = {0}, postcut6_nevent[hcal::gNstamp] = {0};

  //TEST
  vector<int> elastics_per_run;
  vector<int> elastic_runs;
  int total_elastics_allruns = 0;

  //Main loop over runs
  for (Int_t r=0; r<nruns; r++) {

    //Get run experimental parameters
    std::string current_timestamp = runs[r].adcg_ts;

    cout<<"#### current timestamp : "<<current_timestamp<<" ####"<<endl;

    Int_t current_runnumber = runs[r].runnum;

    std::string current_target = runs[r].target;
    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    Int_t mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

    //Get run paths
    //std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment,pass,config,targ_uppercase.c_str());

    //add new path to pass 2 QA data
    std::string exper = experiment;
    // if( exper.compare("gmn")==0 && pass==2 )
    //   rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2/QA_FINAL/SBS%d/LH2/rootfiles/",config);

    //Read in rootfile directory and get path
    std::string rootfile_dir = runs[r].rootfile_dir;
    std::string rootfile_path = rootfile_dir + Form("/*%d*",current_runnumber);

    cout<<"debug : 0000"<<endl;

    //Read in hcal vertical offset per run
    Double_t hcal_v_offset = runs[r].hcal_v_offset;

    //Get target configuration
    SBStarget target_parameters(current_target);
    Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
        
    //check target and continue on deuterium of only calibrating with hydrogen
    if( h2only && target_index!=1 )
      continue;

    //Record old gain and offset parameters by tstamp from database (assumes one file for all timestamps)
    std::string old_db_path = db_path + "/db_sbs.hcal.dat";
    std::string new_adct_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
    if( exper.compare("gmn")==0 && pass==2 ){
      new_adct_path = Form("../timing/parameters/adctsetoffsets_newclus_%s_conf%d_pass%d_0_to_0_exclude_0_to_0.txt",experiment,config,pass);  //Use if timing is not subdivided
      //rootfile_path = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2/QA_FINAL/SBS%d/LH2/rootfiles/*%d*", config, current_runnumber); //only good for GMN
    }

    std::string db_gain_variable = "sbs.hcal.adc.gain";
    std::string db_adct_variable_old = "sbs.hcal.adc.timeoffset";
    std::string db_adct_variable_new = "(inverted) sbs.hcal.adc.timeoffset";
    Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t old_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t new_adcg_coeff[hcal::maxHCalChan] = {1.};
    //Load new and old adctime parameters only if new alignments are available and necessary
    if(new_atime){
      util::readDB( old_db_path, runs[r].adct_ts, db_adct_variable_old, old_adct_offsets );
      //get new adct offsets
      std::string adct_active_timestamp;
      util::tsCompare( config_parameters.GetSBSTimestamp(), runs[r].adct_ts, adct_active_timestamp );
      std::cout << "Loading parameters from timestamp: " << adct_active_timestamp << std::endl;
      util::readDB( new_adct_path, adct_active_timestamp, db_adct_variable_new, new_adct_offsets );
    }
    if( qreplay ){ //Check which timestamp is more current, allowing for hardware changes within a config
      std::string active_timestamp;
      util::tsCompare( current_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
      util::readDB( new_adcgain_path, active_timestamp, db_gain_variable, new_adcg_coeff );
    }

    //Look through available calibration sets by gain coeff timestamp for duplicate entries. Correct index where necessary.
    bool found_ts = false;
    for( Int_t cs=0; cs<=Ncal_set_size; cs++ ){
      if( gain_coeff[cs].timestamp == current_timestamp ){
	found_ts = true;
	Ncal_set_idx = cs;
      }
    }
    
    // If no, add new calibration set and correct the index.
    if( !found_ts ){
      Ncal_set_idx = Ncal_set_size;

      if( Ncal_set_size == hcal::gNstamp-1 ){
	cout << "Error: Allowable calibration sets exceed maximum of " << hcal::gNstamp << ". Check timestamps and adjust limits as necessary." << endl;
	return;
      }

      gain_coeff[Ncal_set_size].timestamp = current_timestamp.c_str();
      gain_coeff[Ncal_set_size].good_events = 0;

      //resize the matrices and vectors
      gain_coeff[Ncal_set_size].Ma.ResizeTo(hcal::maxHCalChan,hcal::maxHCalChan);
      gain_coeff[Ncal_set_size].Ma_oneblock.ResizeTo(hcal::maxHCalChan,hcal::maxHCalChan);
      gain_coeff[Ncal_set_size].ba.ResizeTo(hcal::maxHCalChan);
      gain_coeff[Ncal_set_size].ba_oneblock.ResizeTo(hcal::maxHCalChan);

      cout<<"##### current timestamp : "<<current_timestamp<<" ####"<<endl;

      util::readDB( old_db_path, current_timestamp, db_gain_variable, gain_coeff[Ncal_set_size].old_param );      

      //printout old param
      cout<<"####printing old gain coeff####"<<endl;
      for (int i = 0 ; i<288; i++){
      	cout<<gain_coeff[0].old_param[i]<<" ";
	if ((i+1)%12==0){cout<<endl;}
      }


      //Build quality histograms as necessary from array
      hE[Ncal_set_size]->SetTitle(Form("HCal Cluster E, timestamp: %s; GeV",current_timestamp.c_str()));
      hcoin[Ncal_set_size]->SetTitle(Form("BBCal Shower Atime - HCal Atime, timestamp: %s; ns",current_timestamp.c_str()));
      hW2[Ncal_set_size]->SetTitle(Form("W2 , timestamp: %s; GeV2)",current_timestamp.c_str()));
      hEblk[Ncal_set_size]->SetTitle(Form("HCal Cluster Block E, timestamp: %s; GeV",current_timestamp.c_str()));
      hSF[Ncal_set_size]->SetTitle(Form("HCal Cluster Sampling Fraction, timestamp: %s; %%",current_timestamp.c_str()));
      hXY[Ncal_set_size]->SetTitle(Form("HCal X vs Y (LH2), timestamp: %s; y_{HCAL}; x_{HCAL} (m)",current_timestamp.c_str()));
      hdxdy[Ncal_set_size]->SetTitle(Form("HCal dx vs dy, timestamp: %s; dy_{HCAL}; dx_{HCAL} (m)",current_timestamp.c_str()));
      hdxvmag_h[Ncal_set_size]->SetTitle(Form("Delta X vs Field Setting (LH2), timestamp: %s; field (percent); x_{HCAL} - x_{exp} (m)",current_timestamp.c_str()));
      hdyvmag_h[Ncal_set_size]->SetTitle(Form("Delta Y vs Field Setting (LH2), timestamp: %s; field (percent); y_{HCAL} - y_{exp} (m)",current_timestamp.c_str()));
      hdxvmag_d[Ncal_set_size]->SetTitle(Form("Delta X vs Field Setting (LD2), timestamp: %s; field (percent); x_{HCAL} - x_{exp} (m)",current_timestamp.c_str()));
      hdyvmag_d[Ncal_set_size]->SetTitle(Form("Delta Y vs Field Setting (LD2), timestamp: %s; field (percent); y_{HCAL} - y_{exp} (m)",current_timestamp.c_str()));
      hdxvmag_h3[Ncal_set_size]->SetTitle(Form("Delta X vs Field Setting (He3), timestamp: %s; field (percent); x_{HCAL} - x_{exp} (m)",current_timestamp.c_str()));
      hdyvmag_h3[Ncal_set_size]->SetTitle(Form("Delta Y vs Field Setting (He3), timestamp: %s; field (percent); y_{HCAL} - y_{exp} (m)",current_timestamp.c_str()));
      hEvX[Ncal_set_size]->SetTitle(Form("HCal Cluster E vs X, timestamp: %s; m; GeV",current_timestamp.c_str()));
      hEvY[Ncal_set_size]->SetTitle(Form("HCal Cluster E vs Y, timestamp: %s; m; GeV",current_timestamp.c_str()));
      hSFvID[Ncal_set_size]->SetTitle(Form("HCal Cluster Sampling Fraction vs channel, timestamp: %s; channel; %%",current_timestamp.c_str()));
      hSFvrow[Ncal_set_size]->SetTitle(Form("HCal Cluster Sampling Fraction vs row, timestamp: %s; row; %%",current_timestamp.c_str()));
      hSFvcol[Ncal_set_size]->SetTitle(Form("HCal Cluster Sampling Fraction vs col, timestamp: %s; col; %%",current_timestamp.c_str()));

      Ncal_set_size++;
      
    }

    //quick check on gain coefficients
    //Double_t test_adcg_coeff[hcal::maxHCalChan] = {1.};
    bool different = false;
    if( first )
      std::cout << "First calibration set, Ncal_set_idx: " << Ncal_set_idx << std::endl << std::endl;  
    for( int r=0; r<24; r++ ){
      for ( int c=0; c<12; c++){
	int i = r*12+c;
	if( first ){
	  test_adcg_coeff[i] = new_adcg_coeff[i];
	}else{
	  if( test_adcg_coeff[i] != new_adcg_coeff[i] )
	    different = true;
	}
      }
    }

    if( different ){
      std::cout << "ADC gain coeffients have changed on calibration set " << Ncal_set_idx << std::endl << std::endl;
      for( int r=0; r<24; r++ ){
	for ( int c=0; c<12; c++){
	  int i = r*12+c;
	  test_adcg_coeff[i] = new_adcg_coeff[i];
	}
      }
    }
    std::cout<<"debug 0 "<<endl;
    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    std::vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,1/*verb*/,cut);
    std::cout<<"debug 1 "<<endl;
    if( different || first ){
      std::cout << cut[0];
    }

    std::cout<<"debug : "<<endl; 
    //Remove first calibration set trigger for future processing
    if( first )
      first = false;

    //TEST
    elastic_runs.push_back(current_runnumber);

    // Setting up chain and branch addresses
    C = new TChain("T");
    if(verbose)
      cout << "Adding file to chain: " << rootfile_path << endl;
    C->Add(rootfile_path.c_str());

    C->SetBranchStatus("*",0);
    Double_t SBStr_n[hcal::maxTracks],SBStr_p[hcal::maxTracks];    
    Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
    Double_t BBtr_vz[hcal::maxTracks];
    Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, BBsh_atime;	
    Double_t HCALx, HCALy, HCALe;
    Double_t pblkrow, pblkcol, nblk, nclus;
    Int_t Ncid;
    Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
    Double_t cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
    Double_t HODOtmean;

    // Speed up processing by switching on only useful branches
    C->SetBranchStatus( "*", 0 );
    C->SetBranchStatus( "sbs.hcal.x", 1 );
    C->SetBranchStatus( "sbs.hcal.y", 1 );
    C->SetBranchStatus( "sbs.hcal.e", 1 );
    C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    C->SetBranchStatus( "sbs.hcal.colblk", 1 );
    C->SetBranchStatus( "sbs.hcal.nblk", 1 );
    C->SetBranchStatus( "sbs.hcal.nclus", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    //C->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
    //added by vph
    C->SetBranchStatus("sbs.tr.n",1);
    C->SetBranchStatus("sbs.tr.p",1);
    C->SetBranchStatus("sbs.tr.vz",1);

    C->SetBranchStatus( "bb.tr.n", 1 );
    C->SetBranchStatus( "bb.tr.px", 1 );
    C->SetBranchStatus( "bb.tr.py", 1 );
    C->SetBranchStatus( "bb.tr.pz", 1 );    
    C->SetBranchStatus( "bb.tr.vz", 1 );
    C->SetBranchStatus( "bb.tr.p", 1 );
    C->SetBranchStatus( "bb.ps.e", 1 );
    C->SetBranchStatus( "bb.ps.x", 1 );
    C->SetBranchStatus( "bb.ps.y", 1 );
    C->SetBranchStatus( "bb.sh.e", 1 );
    C->SetBranchStatus( "bb.sh.x", 1 );
    C->SetBranchStatus( "bb.sh.y", 1 );
    C->SetBranchStatus( "bb.sh.atimeblk", 1 );
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "bb.gem.track.nhits", 1 );
    C->SetBranchStatus( "bb.etot_over_p", 1 );
    C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

    // Linking memory
    C->SetBranchAddress( "sbs.hcal.x", &HCALx );
    C->SetBranchAddress( "sbs.hcal.y", &HCALy );
    C->SetBranchAddress( "sbs.hcal.e", &HCALe );
    C->SetBranchAddress( "sbs.hcal.rowblk", &pblkrow );
    C->SetBranchAddress( "sbs.hcal.colblk", &pblkcol );
    C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
    C->SetBranchAddress( "sbs.hcal.nclus", &nclus ); // Total number of clusters
    C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
    C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
    C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
    C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
    //C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // Array of block ADC gain coeff
    C->SetBranchAddress( "sbs.hcal.clus.id", cid );
    C->SetBranchAddress( "sbs.hcal.clus.row", crow );
    C->SetBranchAddress( "sbs.hcal.clus.col", ccol );
    C->SetBranchAddress( "sbs.hcal.clus.e", ce );
    C->SetBranchAddress( "sbs.hcal.clus.x", cx );
    C->SetBranchAddress( "sbs.hcal.clus.y", cy );
    C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
    //added by vph
    C->SetBranchAddress("sbs.tr.n", &SBStr_n);
    C->SetBranchAddress("sbs.tr.p", &SBStr_p);

    C->SetBranchAddress( "bb.tr.n", &BBtr_n );
    C->SetBranchAddress( "bb.tr.px", BBtr_px );
    C->SetBranchAddress( "bb.tr.py", BBtr_py );
    C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
    C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
    C->SetBranchAddress( "bb.tr.p", BBtr_p );
    C->SetBranchAddress( "bb.ps.e", &BBps_e );
    C->SetBranchAddress( "bb.ps.x", &BBps_x );
    C->SetBranchAddress( "bb.ps.y", &BBps_y );
    C->SetBranchAddress( "bb.sh.e", &BBsh_e );
    C->SetBranchAddress( "bb.sh.x", &BBsh_x );
    C->SetBranchAddress( "bb.sh.y", &BBsh_y ); 
    C->SetBranchAddress( "bb.sh.atimeblk", &BBsh_atime ); 
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

    //globalcut enables
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );

    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();

    //Add globalcut and elastic cuts for reporting
    if( target_change_index != target_index ){
      std::cout << target_parameters;
      target_change_index = target_index;
    }

    //Look through available report sets by magnetic field setting and target for duplicate entries. Correct index where necessary.
    bool found_ri = false;
    for( Int_t rs=0; rs<=report_set_size; rs++ ){
      if( report_set[rs].mag == mag && report_set[rs].targetidx == target_index ){
	found_ri = true;
	report_set_idx = rs;
      }
    }

    if( !found_ri ){
      report_set_idx = report_set_size;

      report_set[report_set_idx].gcut = cut[0].gcut;
      report_set[report_set_idx].target = current_target;
      report_set[report_set_idx].mag = mag;
      report_set[report_set_idx].targetidx = target_index;
      report_set[report_set_idx].W2mean = cut[0].W2_mean;
      report_set[report_set_idx].W2sigma = cut[0].W2_sig;
      report_set[report_set_idx].dxmean_n = cut[0].dx0_n;
      report_set[report_set_idx].dxmean_p = cut[0].dx0_p;
      report_set[report_set_idx].dxsigma_n = cut[0].dx_sig_n;
      report_set[report_set_idx].dxsigma_p = cut[0].dx_sig_p;
      report_set[report_set_idx].dymean = cut[0].dy0;
      report_set[report_set_idx].dysigma = cut[0].dy_sig;
      report_set[report_set_idx].atimemean = cut[0].atime0;
      report_set[report_set_idx].atimesigma = cut[0].atime_sig;
      report_set[report_set_idx].minEv = minEventPerCell;
      report_set[report_set_idx].highdelta = highDelta;
    
      report_set_size++;

    }

    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );
    
    // Set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
    Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

    long nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    //keep track of elastics overall
    int total_elastics = 0;

    //Main loop over events in run
    //
    //

    cout<<"####printing old gain coeff####"<<endl;
    for (int j = 0 ; j<2; j++){
      cout<<"Nset : "<<j<<endl;
      for (int i = 0 ; i<288; i++){
        cout<<gain_coeff[j].old_param[i]<<" ";
        if ((i+1)%12==0){cout<<endl;}
      }
    }


    while (C->GetEntry(nevent++)) {

      cout << "Analyzing run " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
      cout.flush();

      ///////
      //Single-loop elastic globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;

      //cout<<"here here"<<endl;

      if( failedglobal ) 
      	continue;

      postcut1_nevent[Ncal_set_idx]++;

      ///////
      //HCal Active Area Cut
      bool failedactivearea = 
	pblkrow==0 || 
	pblkrow==23 || 
	pblkcol==0 || 
	pblkcol==11;

      if( failedactivearea ) 
      	continue; //All events with primary cluster element on edge blocks cut

      postcut2_nevent[Ncal_set_idx]++;

      ///////
      //HCal primary cluster coincidence time cut (using adctime while hcal tdc suspect, new offsets)
      Int_t pblkid = cblkid[0]-1; //define primary block, primary cluster ID

      Double_t natime = cblkatime[0]+old_adct_offsets[pblkid]-new_adct_offsets[pblkid]; //new hcal atime
      Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
      Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time
      Double_t cointime = -BBsh_atime+natime;

      //troubleshooting vimukthi
      h_cointime->Fill(cointime);


      //std::cout<<"cointime : "<<cointime<<endl;
      
      bool failedcoin = abs(cointime-atime0)>atimeNsig*atimesig;
      //use internal time of hcal prior to bbcal timing alignments, pass0/1
      if( (exper.compare("gmn")==0 /*|| exper.compare("gen")==0*/ )&& pass<2 ) //adding gen to this as well
	failedcoin = abs(natime-atime0)>atimeNsig*atimesig;

      hcoin[Ncal_set_idx]->Fill(cointime);

      if( failedcoin ) 
      	continue; //All events where adctime outside of reasonable window cut
	
      postcut3_nevent[Ncal_set_idx]++;

      //ADC arrays reset per event. 
      Double_t A[hcal::maxHCalChan] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
      Double_t A_oneblock[hcal::maxHCalChan] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev

      ///////
      //Physics calculations
      //correct beam energy with vertex information, primary track
      Double_t ebeam_c = ebeam - ( (BBtr_vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );

      TVector3 vertex( 0., 0., BBtr_vz[0] );

      //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t precon = BBtr_p[0] + Eloss_outgoing;

      //set up four-momenta with some empty for various calculation methods
      TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
      TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
      TLorentzVector ptarg; //target momentum
      ptarg.SetPxPyPzE( 0., 0., 0., M_avg );
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      TLorentzVector pN; //N' momentum
      
      //simple calculations for e' and N'
      Double_t etheta = acos( pe.Pz() / pe.E() );
      Double_t ephi = atan2( pe.Py(), pe.Px() );
      Double_t pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
      Double_t phNexp = ephi + hcal::PI;
      Double_t Q2, W2;

      //e' p reconstruction with track angles (not momentum)
      Double_t nu = pbeam.E() - pcent;
      Double_t pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
      Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
      TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
      pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
      Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
      W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;

      //Calculate h-arm quantities
      vector<Double_t> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
      //Double_t KE_p = nu; //For elastics total predicted by earm
      Double_t KE_p = sqrt(SBStr_p[0]*SBStr_p[0]+(0.938*0.938))-0.938; //added by vph: trying something
      Double_t SFrac = HCALe/KE_p; //Measured
      Double_t dx = HCALx - xyhcalexp[0];
      Double_t dy = HCALy - xyhcalexp[1];
      TVector3 neutdir = (hcalpos - vertex).Unit();
      Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + hcal::sbsdipolegap/2.0) );
      TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
      Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
      Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );
	  
      hXY[Ncal_set_idx]->Fill(HCALy,HCALx);

      if( current_target.compare("lh2")==0 ){
	hdxvmag_h[Ncal_set_idx]->Fill( mag, dx );
	hdyvmag_h[Ncal_set_idx]->Fill( mag, dy );
      }
      if( current_target.compare("ld2")==0 ){
	hdxvmag_d[Ncal_set_idx]->Fill( mag, dx );
	hdyvmag_d[Ncal_set_idx]->Fill( mag, dy );
      }
      if( current_target.compare("he3")==0 ){
	hdxvmag_h3[Ncal_set_idx]->Fill( mag, dx );
	hdyvmag_h3[Ncal_set_idx]->Fill( mag, dy );
      }

      //Target Energy in HCal for calibrations
      Double_t hcal_samp_frac = cut[0].hcal_sf; config_sampling_fraction = hcal_samp_frac;
      Double_t hcal_esratio = cut[0].hcal_es; config_e_sigma_ratio = hcal_esratio;
      Double_t KE_exp = KE_p*hcal_samp_frac/hcal_esratio; //Expected E in HCal, equivalent to KE*modified_samp_frac

      ///////
      //BigBite/SBS Acceptance matching. Redundant here with dxdy cuts. Removed.
      bool failedaccmatch = 
	xyhcalexp[1] > hcal::posHCalYf ||
	xyhcalexp[1] < hcal::posHCalYi ||
	xyhcalexp[0] > hcal::posHCalXf ||
	xyhcalexp[0] < hcal::posHCalXi;

      // if( failedaccmatch ) 
      // 	continue;

      ///////
      //dy elastic cut
      Double_t dy0 = cut[0].dy0;
      Double_t dysig = cut[0].dy_sig;
      bool faileddy = abs(dy-dy0)>3*dysig;

      hdxdy[Ncal_set_idx]->Fill(dy,dx);

      if( faileddy ) 
      	continue;

      postcut4_nevent[Ncal_set_idx]++;

      ///////
      //PID
      Int_t pid = -1;
      util::checkPID(current_target,cut[0].dx0_p,cut[0].dx0_n,cut[0].dx_sig_p,cut[0].dx_sig_n,dx,pid);

      ///////
      //dx elastic cut
      if( pid==-1 ) 
      	continue;

      postcut5_nevent[Ncal_set_idx]++;

      ////////////////////////////
      //Primary W2 cut on elastics
      Double_t W2_mean = cut[0].W2_mean;
      Double_t W2_sig = cut[0].W2_sig;
      bool failedW2 = fabs(W2-W2_mean)>W2_sig;

      hW2[Ncal_set_idx]->Fill(W2);

      if( failedW2 ) 
      	continue; //Observed mean W2 cut on elastic peak
      
      postcut6_nevent[Ncal_set_idx]++;
      
      ///////
      //Get corrected primary cluster energy
      Double_t clusE = 0.0;
      Double_t clusblkE = 0.0;
      Int_t badclus = 0;
      //for( Int_t blk = 0; blk<Ncid; blk++ ){
      for( Int_t blk = 0; blk<nblk; blk++ ){
	Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	Double_t blke = cblke[blk];
	
	Double_t blknatime = cblkatime[blk]+old_adct_offsets[blkid]-new_adct_offsets[blkid]; //new atime

	////////
	//Cluster block coincidence adc time check, 5sig estimate
	if( abs(blknatime-atime0)>atimeNsig*atimesig ){
	  badclus = 1;
	  //continue; //May wish to remove this cut. Should investigate expected spread in time in MC
	}
	
	//Correct the block energy with new gain coeff on quasi-replay
	if( qreplay )
	  blke = cblke[blk]/gain_coeff[Ncal_set_idx].old_param[blkid]*new_adcg_coeff[blkid];

	clusE += blke;
	if( blk==0 ) clusblkE += blke;
	
	hEblk[Ncal_set_idx]->Fill( blke );

	////////
	//Fill output tree HCal cluster block variables
	cblkid_out[blk] = blkid;
	cblkatime_out[blk] = cblkatime[blk]+old_adct_offsets[blkid]-new_adct_offsets[blkid];
	cblktime_out[blk] = cblktime[blk];
	cblke_out[blk] = blke;
      }

      //TEST
      total_elastics++;
      total_elastics_allruns++;

      Double_t sampling_fraction = clusE/KE_p; //This must be KE_p NOT KE_exp

      hE[Ncal_set_idx]->Fill( clusE );
      hSF[Ncal_set_idx]->Fill( sampling_fraction );
      hEvX[Ncal_set_idx]->Fill( HCALx, clusE );
      hEvY[Ncal_set_idx]->Fill( HCALy, clusE );
      hSFvID[Ncal_set_idx]->Fill( cblkid[0], sampling_fraction );
      hSFvrow[Ncal_set_idx]->Fill( crow[0], sampling_fraction );
      hSFvcol[Ncal_set_idx]->Fill( ccol[0], sampling_fraction );

      ///////
      //Fill output tree general variables
      pblkid_out = cblkid[0];
      tdc_out = cblktime[0];
      atime_out = natime;
      cointime_out = cointime;
      hcale_out = HCALe;
      pblkcol_out = ccol[0];
      pblkrow_out = crow[0];
      dx_out = dx;
      dy_out = dy;
      HCALx_out = HCALx;
      HCALy_out = HCALy;
      W2_out = W2;
      Q2_out = Q2;
      nu_out = nu;
      hcale_out = clusE;
      sfrac_out = sampling_fraction;
      thetapq_pout = thetapq_p;
      thetapq_nout = thetapq_n;
      nblk_out = nblk;
      pid_out = pid;
      mag_out = mag;
      run_out = current_runnumber;
      tar_out = target_index;
      calset_out = Ncal_set_idx;
      badclus_out = badclus;
      hodotmean_out = HODOtmean;
      failedglobal_out = failedglobal;

      P->Fill();

      ////////
      //Proceed only if calibrating
      if( qreplay ) 
	continue;

      gain_coeff[Ncal_set_idx].good_events++;

      cout<<"good events : "<<gain_coeff[Ncal_set_idx].good_events;
      clusE = 0.0; //reset cluster energy
      Double_t cluspC = 0.0;
      for( Int_t blk = 0; blk<nblk; blk++ ){
      
	Int_t blkid = int(cblkid[blk])-1;
	Double_t blke = cblke[blk];
	//cout<<gain_coeff[Ncal_set_idx].old_param[blkid]<<endl;
	Double_t blkpC = cblke[blk]/gain_coeff[Ncal_set_idx].old_param[blkid]; //divide off the old coefficient
	Double_t blknatime = cblkatime[blk]+old_adct_offsets[blkid]-new_adct_offsets[blkid]; //new atime

	////////
	//Cluster block coincidence time cut (using adct for now)
	bool badblock = false;
	if( abs(blknatime-atime0)>atimeNsig*atimesig )
	  badblock=true;
	    
	//if( badblock ) //May wish to remove this cut until adc time is better calibrated
	   //continue;

	clusE += blke;
	cluspC += blkpC;
	A[blkid] += blkpC;
	if( nblk==1 ){
	  A_oneblock[blkid] += blkpC;
	  gain_coeff[Ncal_set_idx].err_ev_oneblock[blkid] += pow( ((1.+0.5*blkpC)/blkpC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  gain_coeff[Ncal_set_idx].NEV_oneblock[blkid]++;
	}
	gain_coeff[Ncal_set_idx].err_ev[blkid] += pow( ((1.+blkpC)/blkpC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean
	gain_coeff[Ncal_set_idx].NEV[blkid]++;
	    
      }

      //Build the 288x288 matrix as clearly as possible
      for(Int_t icol = 0; icol<hcal::maxHCalChan; icol++){ //matrix column here
	gain_coeff[Ncal_set_idx].ba(icol) += A[icol];

	if( Ncid==1 ) 
	  gain_coeff[Ncal_set_idx].ba_oneblock(icol) += A_oneblock[icol];

	for(Int_t irow = 0; irow<hcal::maxHCalChan; irow++){ //matrix row here
	  gain_coeff[Ncal_set_idx].Ma(icol,irow) += A[icol]*A[irow]/KE_exp;
	  if( Ncid==1 ){	    
	    gain_coeff[Ncal_set_idx].Ma_oneblock(icol,irow) += A_oneblock[icol]*A_oneblock[irow]/KE_exp;
	  } //endloop oneblock
	} //endloop over matrix element rows
      } //endloop over matrix element cols

    }//endloop over events

    //TEST
    elastics_per_run.push_back(total_elastics);

    // getting ready for the next run
    C->Reset();

  }//endloop over runs

  cout << "Ended loop over runs. " << Ncal_set_size << " calibration sets for " << experiment << " configuration " << config << ". Proceeding to final analysis." << endl; 

  //Declare coeff outfile
  ofstream GainCoeff;
  if( qreplay )
    new_adcgain_path = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay
    
  //Loop over all independent data sets
  for( Int_t s=0; s<Ncal_set_size; s++ ){
    
    //check if timestamp for calibration set is newer than config timestamp. Use the newer of the two.
    std::string active_timestamp;
    util::tsCompare(sbs_timestamp,gain_coeff[s].timestamp,active_timestamp);

    //Reject the bad cells and normalize the oneblock check
    Int_t badcell[hcal::maxHCalChan];
    Int_t badcell_oneblock[hcal::maxHCalChan];
    Int_t cellBad = 0;
    Double_t y[hcal::maxHCalChan] = {0.0}; // For easy TGraphErrors build
  
    for(Int_t i=0; i<hcal::maxHCalChan; i++){
      if( qreplay ){
	break;
      }

      badcell[i] = 0;
      y[i] = i;
    
      //Do not change ADC gain coeff if insufficient events or energy dep in cell
      if( gain_coeff[s].NEV[i] < minEventPerCell || gain_coeff[s].Ma(i,i) < highDelta*gain_coeff[s].ba(i) ){ 

	cellBad = 1;

	Double_t elemRatio = gain_coeff[s].Ma(i,i)/gain_coeff[s].ba(i);

	gain_coeff[s].ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
	gain_coeff[s].Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
	for(Int_t j=0; j<hcal::maxHCalChan; j++){
	  if( j != i ){
	    gain_coeff[s].Ma(i,j) = 0.0;
	    gain_coeff[s].Ma(j,i) = 0.0;
	  }
	}

	badcell[i] = 1;
      } 

      badcell_oneblock[i] = 0;

      //Do not change ADC gain coeff if insufficient events or energy dep in cell, oneblock
      if( gain_coeff[s].NEV_oneblock[i] < minEventPerCell || gain_coeff[s].Ma_oneblock(i,i) < highDelta*gain_coeff[s].ba_oneblock(i) ){ 

	cellBad = 1;

	Double_t elemRatio = gain_coeff[s].Ma_oneblock(i,i)/gain_coeff[s].ba_oneblock(i);

	gain_coeff[s].ba_oneblock(i) = 1.0;  // Set RHS vector for cell i to 1.0 
	gain_coeff[s].Ma_oneblock(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
	for(Int_t j=0; j<hcal::maxHCalChan; j++){
	  if( j != i ){
	    gain_coeff[s].Ma_oneblock(i,j) = 0.0; 
	    gain_coeff[s].Ma_oneblock(j,i) = 0.0;
	  }
	}
	badcell_oneblock[i] = 1;
      }   

      //Calculate error per block (element of A)
      gain_coeff[s].err[i] = sqrt(gain_coeff[s].err_ev[i]/gain_coeff[s].NEV[i]);
      gain_coeff[s].err_oneblock[i] = sqrt(gain_coeff[s].err_ev_oneblock[i]/gain_coeff[s].NEV_oneblock[i]);
    } //endloop over channels
  
    if( cellBad!=0 && !qreplay ) cout << "Bad cells detected on set " << gain_coeff[s].timestamp << ". See report for details." << endl << endl;

    //If iteration == 0, invert the matrix, solve for ratios
    if( !qreplay ){
      TMatrixD M_inv = gain_coeff[s].Ma.Invert();
      TMatrixD M_inv_oneblock = gain_coeff[s].Ma_oneblock.Invert();
      TVectorD Coeff = M_inv*gain_coeff[s].ba; // Stays unmodified for reference
      TVectorD Coeff_oneblock = M_inv_oneblock*gain_coeff[s].ba_oneblock; // Stays unmodified for reference

      for(Int_t i=0; i<hcal::maxHCalChan; i++){

	if(badcell_oneblock[i]==0){
	  gain_coeff[s].new_param_oneblock[i]=Coeff_oneblock[i];
	}else{
	  gain_coeff[s].new_param_oneblock[i]=gain_coeff[s].old_param[i];
	}

	if(badcell[i]==0&&Coeff[i]>0){ //Only update the coefficent if coefficient passes quality checks
	  gain_coeff[s].new_param[i]=Coeff[i];
	  gain_coeff[s].new_param_divide[i]=gain_coeff[s].new_param[i]/gain_coeff[s].new_param_oneblock[i];
	}else{
	  gain_coeff[s].new_param[i]=gain_coeff[s].old_param[i]; // If the cell calibration is bad, use the old coefficient
	  gain_coeff[s].new_param_divide[i]=-1.0;
	}
      }

      Double_t yErr[hcal::maxHCalChan] = {0.0};

      TGraphErrors *ccgraph = new TGraphErrors( hcal::maxHCalChan, y, gain_coeff[s].new_param, yErr, gain_coeff[s].err ); 
      ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
      ccgraph->GetYaxis()->SetLimits(0.0,0.25);
      ccgraph->SetTitle(Form("Calibration Coefficients, %s",active_timestamp.c_str()));
      ccgraph->GetXaxis()->SetTitle("Channel");
      ccgraph->GetYaxis()->SetTitle("Unitless");
      ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
      ccgraph->Write(Form("constants_%d",s));

      TGraphErrors *ccgraph_oneBlock = new TGraphErrors( hcal::maxHCalChan, y, gain_coeff[s].new_param_oneblock, yErr, gain_coeff[s].err_oneblock ); 
      ccgraph_oneBlock->GetXaxis()->SetLimits(0.0,288.0);  
      ccgraph_oneBlock->GetYaxis()->SetLimits(0.0,0.25);
      ccgraph_oneBlock->SetTitle(Form("Calibration Coefficients One Block, %s",active_timestamp.c_str()));
      ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
      ccgraph_oneBlock->GetYaxis()->SetTitle("Unitless");
      ccgraph_oneBlock->SetMarkerStyle(21);
      ccgraph_oneBlock->Write(Form("constants_oneblock_%d",s));

      TGraphErrors *ccgraph_divide = new TGraphErrors( hcal::maxHCalChan, y, gain_coeff[s].new_param_divide, yErr, yErr ); 
      ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
      ccgraph_divide->SetTitle(Form("Calibration Coefficients / OneBlock Coeff, %s",active_timestamp.c_str()));
      ccgraph_divide->GetXaxis()->SetTitle("Channel");
      ccgraph_divide->GetYaxis()->SetTitle("Unitless");
      ccgraph_divide->SetMarkerStyle(21);
      ccgraph_divide->Write(Form("constants_divide_%d",s));

      if( s==0 )
	GainCoeff.open( new_adcgain_path );
      else{
	GainCoeff.open( new_adcgain_path, fstream::app ); //append to the end of the file after first cal set
	GainCoeff << endl << endl;
      }

      //Console/txt outs
      Int_t cell = 0;

      GainCoeff << "#HCal gain coefficients from SBS-" << config << " obtained " << date.c_str() << " timestamp set " << active_timestamp << endl;
      GainCoeff << "Total number of events available to calibrate for this set: " << gain_coeff[s].good_events << endl;
      if( active_timestamp.compare("none")!=0 )
	GainCoeff << active_timestamp << endl;
      GainCoeff << "sbs.hcal.adc.gain =" << endl;

      cout << "Gain Coefficients" << endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for( Int_t c=0; c<hcal::maxHCalCols; c++ ){
	  GainCoeff << gain_coeff[s].new_param[cell] << "  ";
	  cout << gain_coeff[s].new_param[cell] << "  ";
	  cell++;
	}
	GainCoeff << endl;
	cout << endl;
      }

      cell = 0;

      cout << endl << "One Block Std:" << endl;
      GainCoeff << endl << "#One Block Std = " << endl;

      for( Int_t r = 0; r<hcal::maxHCalRows; r++){
	for( Int_t c = 0; c<hcal::maxHCalCols; c++){
	  GainCoeff << gain_coeff[s].new_param_oneblock[cell] << "  ";
	  cout << gain_coeff[s].new_param_oneblock[cell] << "  ";
	  cell++;
	}
	GainCoeff << endl;
	cout << endl;
      }

      cell = 0;

      cout << endl << "Total Number of events available for calibration: " << endl << endl;
      GainCoeff << endl << "#Number of events available for calibration = " << endl;

      for( Int_t r = 0; r<hcal::maxHCalRows; r++){
	for( Int_t c = 0; c<hcal::maxHCalCols; c++){
	  GainCoeff << gain_coeff[s].NEV[cell] << "  ";
	  cout << gain_coeff[s].NEV[cell] << "  ";
	  cell++;
	}
	GainCoeff << endl;
	cout << endl;
      }
      GainCoeff << endl << endl;
      GainCoeff.close();
    } //end if iter==0

  }//endloop over calibration sets

  //Add output report canvas
  TCanvas *c1[report_set_size];
  Int_t report_height = 1800;

  for( Int_t s=0; s<report_set_size; s++ ){

    c1[s] = new TCanvas(Form("ecalreport_report%d",s), Form("Configuration/Cut Information, Report %d",s), report_height, 1000);
  
    // Set margin.
    c1[s]->SetLeftMargin(0.01);

    // Create a TText object.
    TText *t = new TText();

    // Set text align to the left (horizontal alignment = 1).
    t->SetTextAlign(11);
    t->SetTextSize(0.02);

    //make an array of strings
    std::string report[linecount] = {
      "General HCal Energy Calibration Info",
      Form("Experiment: %s, Configuration: %d, Pass: %d", experiment, config, pass),
      Form("Creation Date: %s", date.c_str() ),
      Form("Target: %s", report_set[s].target.c_str() ),
      Form("SBS Field: %d%%", report_set[s].mag ),
      "",
      "Elastic Cuts",
      Form("Global Elastic Cuts: %s", report_set[s].gcut.c_str() ),
      Form("W2 mean (GeV): %f", report_set[s].W2mean),
      Form("W2 sigma (GeV): %f", report_set[s].W2sigma),
      Form("dx mean, neutron (m): %f", report_set[s].dxmean_n),
      Form("dx mean, proton (m): %f", report_set[s].dxmean_p),
      Form("dx sigma, neutron (m): %f", report_set[s].dxsigma_p),
      Form("dx sigma, proton (m): %f", report_set[s].dxsigma_p),
      Form("dy mean (m): %f", report_set[s].dymean),
      Form("dy sigma (m): %f", report_set[s].dysigma),
      Form("adc time mean (ns): %f", report_set[s].atimemean),
      Form("adc time sigma (ns): %f", report_set[s].atimesigma),
      "",
      "Other Cuts/Information",
      Form("Minimum Ev per Cell : %d", report_set[s].minEv),
      Form("Minimum Energy Deposited in Cell (factor, vs expectation) : %0.2f", report_set[s].highdelta),
      Form("Sampling Fraction Target from Monte Carlo: %f", config_sampling_fraction),
      Form("Observed Energy to Energy Sigma Ratio: %f", config_e_sigma_ratio),
      "HCal Active Area (Projected Nucleon 1 row/col Within HCal Acceptance)"
    };
    // Loop to write the lines to the canvas.
    for( Int_t i = 0; i<linecount; i++ ) {
      // Vertical position adjusted according to line number.
      Double_t verticalPosition = 0.95 - i * 0.03;
      t->DrawTextNDC(0.1, verticalPosition, report[i].c_str());
    }

    c1[s]->Write();    
  }

  //clean up
  for( Int_t unused_set_index=Ncal_set_size; unused_set_index<hcal::gNstamp; unused_set_index++ ){
    delete hE[unused_set_index];
    delete hEblk[unused_set_index];
    delete hSF[unused_set_index];
    delete hdxvmag_h[unused_set_index];
    delete hdyvmag_h[unused_set_index];
    delete hdxvmag_d[unused_set_index];
    delete hdyvmag_d[unused_set_index];
    delete hdxvmag_h3[unused_set_index];
    delete hdyvmag_h3[unused_set_index];
    delete hEvX[unused_set_index];
    delete hEvY[unused_set_index];
    delete hSFvID[unused_set_index];
    delete hSFvrow[unused_set_index];
    delete hSFvcol[unused_set_index];
  }

  fout->Write();

  st->Stop();

  //TEST
  if( elastics_per_run.size() != elastic_runs.size() ){
    cout << "TEST ERROR: vector size mismatch" << endl;
    return;
  }

  //TEST
  for( int i=0; i<elastic_runs.size(); ++i ){
    cout << "run:" << elastic_runs[i] << "  events: " << elastics_per_run[i] << endl;

  }
  cout << "Total elastics: " << total_elastics_allruns << endl;

  for( Int_t s=0; s<report_set_size; ++s ){
    cout << "Cut detail set " << s << ", post cut 1:2:3:4:5:6, " << postcut1_nevent[s] << ":" << postcut2_nevent[s] << ":" << postcut3_nevent[s] << ":" << postcut4_nevent[s] << ":" << postcut5_nevent[s] << ":" << postcut6_nevent[s] << endl;
  }

  cout << endl << endl << "Analysis complete. Output file located " << ecal_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    
  
  //h_cointime->Draw();

}
