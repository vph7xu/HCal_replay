//sseeds 11.23.23 - Updated parsing script to apply wide and certain globalcuts and to parse back branches to inform only calibration cuts for hcal energy.

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

const Int_t maxClus = 35; //larger than max per tree
const Int_t maxBlk = 25;
const Int_t maxTrack = 1000;

//MAIN
void hydrogen_cut_diagnostic( const char *experiment = "gen", Int_t kine=3, Int_t pass=1, bool verb=false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //One set of data files in json config shared between pass0/1 per kinematic
  if( pass==0 )
    pass=1;

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s_p%d/",experiment,pass); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,kine,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available

  // outfile path
  std::string out_dir = "/work/halla/sbs/vimukthi/HCal_replay/hcal";//gSystem->Getenv("OUT_DIR");
  std::string out_path = out_dir + Form("/hcal_calibrations/SBS/energy/diagnostics/caldiag_%s_sbs%d_pass%d.root",experiment,kine,pass);

  //set up output files
  TFile *fout = new TFile( out_path.c_str(), "RECREATE" );

  //set up diagnostic histograms
  TH2D *hW2mag = new TH2D( "hW2mag", "W^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 2 );
  TH2D *hQ2mag = new TH2D( "hQ2mag", "Q^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 4 ); //can generalize
  TH2D *hdxmag = new TH2D( "hdxmag","dx vs sbsmag; \%; x_{HCAL}-x_{expect} (m)", 20, 0, 100, 800, -4, 4 );
  TH2D *hdymag = new TH2D( "hdymag","dy vs sbsmag; \%; y_{HCAL}-y_{expect} (m)", 20, 0, 100, 800, -4, 4 );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  TChain *TS = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // new output tree vars
  Double_t dx_out;
  Double_t dy_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t nu_out;
  Double_t NT_out;
  Double_t NT_corr_out;
  Double_t eoverp_out;
  Double_t thetapq_out;
  Double_t charge_ts_out;
  Double_t epcent_out;
  Double_t eprecon_out;
  Double_t adc_coin_out; //BBCal/HCal ADC coincidence time
  Int_t run_out;
  Int_t tar_out; //0:LH2, 1:LD2
  Int_t mag_out;
  Int_t event_out;
  Int_t trig_out;
  Double_t charge_avg_out;
  Double_t comp_event_frac_out;
  Double_t event_frac_out;
  Double_t acc_charge_out;

  // relevant old output tree vars
  Double_t bb_tr_vx_out;
  Double_t bb_tr_vy_out;
  Double_t bb_tr_vz_out;
  Double_t bb_tr_p_out;
  Double_t bb_tr_px_out;
  Double_t bb_tr_py_out;
  Double_t bb_tr_pz_out;
  Double_t bb_tr_n_out;
  Double_t e_kine_Q2_out;
  Double_t e_kine_W2_out;
  Double_t e_kine_nu_out;
  Double_t bb_ps_e_out;
  Double_t bb_ps_rowblk_out;
  Double_t bb_ps_colblk_out;
  Double_t bb_sh_e_out;
  Double_t bb_sh_rowblk_out;
  Double_t bb_sh_colblk_out;
  Double_t bb_sh_atimeblk_out;
  Double_t bb_sh_nclus_out;
  Double_t bb_hodotdc_clus_tmean_out;
  Double_t bb_gem_track_nhits_out;
  Double_t bb_etot_over_p_out;
  Double_t sbs_hcal_e_out;
  Double_t sbs_hcal_x_out;
  Double_t sbs_hcal_y_out;
  Double_t sbs_hcal_rowblk_out;
  Double_t sbs_hcal_colblk_out;
  Double_t sbs_hcal_atimeblk_out;
  Double_t sbs_hcal_tdctimeblk_out;
  Double_t sbs_hcal_clus_id_out[maxClus];
  Double_t sbs_hcal_clus_e_out[maxClus];
  Double_t sbs_hcal_clus_x_out[maxClus];
  Double_t sbs_hcal_clus_y_out[maxClus];
  Double_t sbs_hcal_clus_tdctime_out[maxClus];
  Double_t sbs_hcal_clus_atime_out[maxClus];
  Double_t sbs_hcal_nclus_out;
  Double_t sbs_hcal_nblk_out;
  Double_t sbs_hcal_clus_blk_id_out[maxBlk];
  Double_t sbs_hcal_clus_blk_e_out[maxBlk];
  Double_t sbs_hcal_clus_blk_x_out[maxBlk];
  Double_t sbs_hcal_clus_blk_y_out[maxBlk];
  Double_t sbs_hcal_clus_blk_atime_out[maxBlk];
  Double_t sbs_hcal_clus_blk_tdctime_out[maxBlk];
  Int_t Ndata_sbs_hcal_clus_blk_id_out;
  Int_t Ndata_sbs_hcal_clus_id_out;

  // set new output tree branches
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "NT", &NT_out, "NT/D" );
  P->Branch( "NT_corr", &NT_corr_out, "NT_corr/D" );
  P->Branch( "eoverp", &eoverp_out, "eoverp/D" );
  P->Branch( "thetapq", &thetapq_out, "thetapq/D" );
  P->Branch( "charge_ts", &charge_ts_out, "charge_ts/D" );
  P->Branch( "epcent", &epcent_out, "epcent/D" );
  P->Branch( "eprecon", &eprecon_out, "eprecon/D" );
  P->Branch( "adc_coin", &adc_coin_out, "adc_coin/D" );
  P->Branch( "mag", &mag_out, "mag/I" );
  P->Branch( "run", &run_out, "run/I" );
  P->Branch( "tar", &tar_out, "tar/I" );
  P->Branch( "event", &event_out, "event/I" );
  P->Branch( "trig", &trig_out, "trig/I" );
  P->Branch( "charge_avg", &charge_avg_out, "charge_avg/D" );
  P->Branch( "comp_event_frac", &comp_event_frac_out, "comp_event_frac/D" );
  P->Branch( "event_frac", &event_frac_out, "event_frac/D" );

  // set relevant old output tree branches
  P->Branch( "bb.tr.vx", &bb_tr_vx_out, "bb.tr.vx/D" );
  P->Branch( "bb.tr.vy", &bb_tr_vy_out, "bb.tr.vy/D" );
  P->Branch( "bb.tr.vz", &bb_tr_vz_out, "bb.tr.vz/D" );
  P->Branch( "bb.tr.p", &bb_tr_p_out, "bb.tr.p/D" );
  P->Branch( "bb.tr.px", &bb_tr_px_out, "bb.tr.px/D" );
  P->Branch( "bb.tr.py", &bb_tr_py_out, "bb.tr.py/D" );
  P->Branch( "bb.tr.pz", &bb_tr_pz_out, "bb.tr.pz/D" );
  P->Branch( "bb.tr.n", &bb_tr_n_out, "bb.tr.n/D" );
  P->Branch( "e.kine.Q2", &e_kine_Q2_out, "e.kine.Q2/D" );
  P->Branch( "e.kine.W2", &e_kine_W2_out, "e.kine.W2/D" );
  P->Branch( "e.kine.nu", &e_kine_nu_out, "e.kine.nu/D" );
  P->Branch( "bb.ps.e", &bb_ps_e_out, "bb.ps.e/D" );
  P->Branch( "bb.ps.rowblk", &bb_ps_rowblk_out, "bb.ps.rowblk/D" );
  P->Branch( "bb.ps.colblk", &bb_ps_colblk_out, "bb.ps.colblk/D" );
  P->Branch( "bb.sh.e", &bb_sh_e_out, "bb.sh.e/D" );
  P->Branch( "bb.sh.rowblk", &bb_sh_rowblk_out, "bb.sh.rowblk/D" );
  P->Branch( "bb.sh.colblk", &bb_sh_colblk_out, "bb.sh.colblk/D" );
  P->Branch( "bb.sh.atimeblk", &bb_sh_atimeblk_out, "bb.sh.atimeblk/D" );
  P->Branch( "bb.sh.nclus", &bb_sh_nclus_out, "bb.sh.nclus/D" );
  P->Branch( "bb.hodotdc.clus.tmean", &bb_hodotdc_clus_tmean_out, "bb.hodotdc.clus.tmean/D" );
  P->Branch( "bb.gem.track.nhits", &bb_gem_track_nhits_out, "bb.gem.track.nhits/D" );
  P->Branch( "bb.etot_over_p", &bb_etot_over_p_out, "bb.etot_over_p/D" );
  P->Branch( "sbs.hcal.e", &sbs_hcal_e_out, "sbs.hcal.e/D" );
  P->Branch( "sbs.hcal.x", &sbs_hcal_x_out, "sbs.hcal.x/D" );
  P->Branch( "sbs.hcal.y", &sbs_hcal_y_out, "sbs.hcal.y/D" );
  P->Branch( "sbs.hcal.rowblk", &sbs_hcal_rowblk_out, "sbs.hcal.rowblk/D" );
  P->Branch( "sbs.hcal.colblk", &sbs_hcal_colblk_out, "sbs.hcal.colblk/D" );
  P->Branch( "sbs.hcal.atimeblk", &sbs_hcal_atimeblk_out, "sbs.hcal.atimeblk/D" );
  P->Branch( "sbs.hcal.clus.id", &sbs_hcal_clus_id_out, "sbs.hcal.clus.id/D" );
  P->Branch( "sbs.hcal.clus.e", &sbs_hcal_clus_e_out, "sbs.hcal.clus.e/D" );
  P->Branch( "sbs.hcal.clus.x", &sbs_hcal_clus_x_out, "sbs.hcal.clus.x/D" );
  P->Branch( "sbs.hcal.clus.y", &sbs_hcal_clus_y_out, "sbs.hcal.clus.y/D" );
  P->Branch( "sbs.hcal.clus.tdctime", &sbs_hcal_clus_tdctime_out, "sbs.hcal.clus.tdctime/D" );
  P->Branch( "sbs.hcal.clus.atime", &sbs_hcal_clus_atime_out, "sbs.hcal.clus.atime/D" );
  P->Branch( "sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id_out, "sbs.hcal.clus_blk.id/D" );
  P->Branch( "sbs.hcal.clus_blk.e", &sbs_hcal_clus_blk_e_out, "sbs.hcal.clus_blk.e/D" );
  P->Branch( "sbs.hcal.clus_blk.x", &sbs_hcal_clus_blk_x_out, "sbs.hcal.clus_blk.x/D" );
  P->Branch( "sbs.hcal.clus_blk.y", &sbs_hcal_clus_blk_y_out, "sbs.hcal.clus_blk.y/D" );
  P->Branch( "sbs.hcal.clus_blk.tdctime", &sbs_hcal_clus_blk_tdctime_out, "sbs.hcal.clus_blk.tdctime/D" );
  P->Branch( "sbs.hcal.clus_blk.atime", &sbs_hcal_clus_blk_atime_out, "sbs.hcal.clus_blk.atime/D" );
  P->Branch( "sbs.hcal.nclus", &sbs_hcal_nclus_out, "sbs.hcal.nclus/D" );
  P->Branch( "sbs.hcal.nblk", &sbs_hcal_nblk_out, "sbs.hcal.nblk/D" );
  P->Branch( "Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id_out, "Ndata.sbs.hcal.clus.id/I" );
  P->Branch( "Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id_out, "Ndata.sbs.hcal.clus_blk.id/I" );

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,kine);    
  cout << config_parameters;
  Double_t hcaltheta = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbthr = config_parameters.GetBBtheta_rad(); //in radians

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";

  for (Int_t irun=0; irun<nruns; irun++) {
      
    // only analyze hydrogen data to minimize better isolate elastic selection cuts
    std::string targ = runs[irun].target;
    if( targ.compare("lh2")!=0 )
      continue;
    
    Int_t t = 0; //only lh2

    std::string targ_uppercase = targ; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    
    // get additional information from run object
    Int_t runnum = runs[irun].runnum;
    Int_t mag = runs[irun].sbsmag / 21; //convert to percent where max field is at 2100A
    Double_t ebeam = runs[irun].ebeam; //get beam energy per run
    Double_t charge = runs[irun].charge; //get total charge per run
    
    //Read in rootfile directory and get path
    std::string rootfile_dir = runs[irun].rootfile_dir;
    std::string rootfile_path = rootfile_dir + Form("/*%d*",runnum);

    //Get target parameters
    SBStarget target_parameters(targ);
    Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
    
    //Read in hcal vertical offset per run
    Double_t hcal_v_offset = runs[irun].hcal_v_offset;
    
    //analyze bcm from TSsbs for charge accumulated per recorded event number
    typedef struct{
      Double_t enumber;
      Double_t clock_diff;
      Double_t avg_current;
      Double_t acc_charge;
    } SCALARAGG;
    
    vector<SCALARAGG> Sagg;
    
    //Create chain for epics/scalar tree
    TS = new TChain("TSsbs");
    TS->Add( rootfile_path.c_str() );
    
    //Declare dummies and counters for needed vars
    Double_t Sevnum = 0., Sdnew = 0., Sclk = 0., clkLast = 0., curLast = 0., curAvg = 0., clkDiff = 0.;
    Int_t DumTN = 0, DumCurTN = 0, DumEvLast = 0, DumEvOff = 0, DumRun = 0, evDiff = 0, evTCoff = 0;
    Long64_t CurEv = 0;
    Long64_t NumEv = TS->GetEntries();
    bool treechange = false;
    bool begin = true;
    Double_t chargeLast = 0;
    Double_t chargeAvgLast = 0;
    //Switch on and set addresses for branches
    TS->SetBranchStatus( "*", 0 );
    TS->SetBranchStatus( "evNumber", 1 );
    TS->SetBranchStatus( "sbs.bcm.dnew.cnt", 1 );
    TS->SetBranchStatus( "sbs.104kHz_CLK.cnt", 1 );
    TS->SetMakeClass(1);
    TS->SetBranchAddress( "evNumber", &Sevnum );
    TS->SetBranchAddress( "sbs.bcm.dnew.cnt", &Sdnew );
    TS->SetBranchAddress( "sbs.104kHz_CLK.cnt", &Sclk );
    //Loop over events (dt=2s for epics tree)
    //Get rate (Sdnew/3318/clkLast) then multiply by the time (clkLast), without the offset, 3318 is cnts/C
    while( TS->GetEntry( CurEv++ ) ){
      
      Double_t clk = Sclk/103700;
      if(begin)
	clkLast=clk;
      clkDiff = clk-clkLast;
      clkLast = clk;
      
      Double_t cur = Sdnew/3318;
      if(begin){
	curLast=cur;
	begin=false;
      }
      curAvg = (cur+curLast)/2;
      curLast = cur;
      
      chargeLast = curAvg*clkDiff*10E-9; //charge = current*time, convert from ns->s
      
      //Fill the scalar data structure with charge and clock corrections applied
      SCALARAGG thisSCALAR = { Sevnum, clkDiff, curAvg, chargeLast };
      Sagg.push_back( thisSCALAR );
      
    }

    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    vector<caldiag> cut;
    util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ,mag,verb,cut);  
    
    std::string gcut = cut[0].gcut;
    Double_t W2min = cut[0].W2_min;
    Double_t W2max = cut[0].W2_max;
    Double_t dy0 = cut[0].dy0;
    Double_t dysig = cut[0].dy_sig;
    Double_t sffac = cut[0].sf_corr;
    
    C = new TChain("T");
    C->Add(rootfile_path.c_str());
    
    //cout << endl << rootfile_path << endl << endl;

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    
    C->SetMakeClass(1); //For event header branches

    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, nclus, nblk;

    C->SetBranchStatus( "sbs.hcal.idblk", 1 );
    C->SetBranchStatus( "sbs.hcal.e", 1 );
    C->SetBranchStatus( "sbs.hcal.x", 1 );
    C->SetBranchStatus( "sbs.hcal.y", 1 );
    C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    C->SetBranchStatus( "sbs.hcal.colblk", 1 );
    C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );
    C->SetBranchStatus( "sbs.hcal.atimeblk", 1 );
    C->SetBranchStatus( "sbs.hcal.index", 1 );
    C->SetBranchStatus( "sbs.hcal.nclus", 1 );
    C->SetBranchStatus( "sbs.hcal.nblk", 1 );

    C->SetBranchAddress( "sbs.hcal.idblk", &hcalid );
    C->SetBranchAddress( "sbs.hcal.e", &hcale );
    C->SetBranchAddress( "sbs.hcal.x", &hcalx );
    C->SetBranchAddress( "sbs.hcal.y", &hcaly );
    C->SetBranchAddress( "sbs.hcal.rowblk", &hcalr );
    C->SetBranchAddress( "sbs.hcal.colblk", &hcalc );
    C->SetBranchAddress( "sbs.hcal.tdctimeblk", &hcaltdc );
    C->SetBranchAddress( "sbs.hcal.atimeblk", &hcalatime );
    C->SetBranchAddress( "sbs.hcal.index", &hcalidx );
    C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
    C->SetBranchAddress( "sbs.hcal.nblk", &nblk );

    // HCal cluster branches
    Double_t hcalcid[maxClus], hcalce[maxClus], hcalcx[maxClus], hcalcy[maxClus], hcalctdctime[maxClus], hcalcatime[maxClus], hcalcrow[maxClus], hcalccol[maxClus];
    Int_t Nhcalcid;

    C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 );

    C->SetBranchAddress( "sbs.hcal.clus.id", hcalcid );
    C->SetBranchAddress( "sbs.hcal.clus.e", hcalce );
    C->SetBranchAddress( "sbs.hcal.clus.x", hcalcx );
    C->SetBranchAddress( "sbs.hcal.clus.y", hcalcy );
    C->SetBranchAddress( "sbs.hcal.clus.tdctime", hcalctdctime );
    C->SetBranchAddress( "sbs.hcal.clus.atime", hcalcatime );
    C->SetBranchAddress( "sbs.hcal.clus.row", hcalcrow );
    C->SetBranchAddress( "sbs.hcal.clus.col", hcalccol );
    C->SetBranchAddress( "sbs.hcal.clus.id", hcalcid );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Nhcalcid );

    // HCal cluster blk branches
    Double_t hcalcbid[maxClus], hcalcbe[maxClus], hcalcbx[maxClus], hcalcby[maxClus], hcalcbtdctime[maxClus], hcalcbatime[maxClus];
    Int_t Nhcalcbid;
    
    C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.x", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.y", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    C->SetBranchStatus( "Ndata.sbs.hcal.clus_blk.id", 1 );
   
    C->SetBranchAddress( "sbs.hcal.clus_blk.id", hcalcbid );
    C->SetBranchAddress( "sbs.hcal.clus_blk.e", hcalcbe );
    C->SetBranchAddress( "sbs.hcal.clus_blk.x", hcalcbx );
    C->SetBranchAddress( "sbs.hcal.clus_blk.y", hcalcby );
    C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", hcalcbtdctime );
    C->SetBranchAddress( "sbs.hcal.clus_blk.atime", hcalcbatime );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus_blk.id", &Nhcalcbid );

    // bbcal clus var
    Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, nclusSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;

    C->SetBranchStatus( "bb.sh.e", 1 );
    C->SetBranchStatus( "bb.sh.x", 1 );
    C->SetBranchStatus( "bb.sh.y", 1 );
    C->SetBranchStatus( "bb.sh.rowblk", 1 );
    C->SetBranchStatus( "bb.sh.colblk", 1 );
    C->SetBranchStatus( "bb.sh.idblk", 1 );
    C->SetBranchStatus( "bb.sh.atimeblk", 1 );
    C->SetBranchStatus( "bb.sh.nclus", 1 );
    C->SetBranchStatus( "bb.ps.e", 1 );
    C->SetBranchStatus( "bb.ps.rowblk", 1 );
    C->SetBranchStatus( "bb.ps.colblk", 1 );
    C->SetBranchStatus( "bb.ps.idblk", 1 );
    C->SetBranchStatus( "bb.ps.atimeblk", 1 );

    C->SetBranchAddress( "bb.sh.e", &eSH );
    C->SetBranchAddress( "bb.sh.x", &xSH );
    C->SetBranchAddress( "bb.sh.y", &ySH );
    C->SetBranchAddress( "bb.sh.rowblk", &rblkSH );
    C->SetBranchAddress( "bb.sh.colblk", &cblkSH );
    C->SetBranchAddress( "bb.sh.idblk", &idblkSH );
    C->SetBranchAddress( "bb.sh.atimeblk", &atimeSH );
    C->SetBranchAddress( "bb.sh.nclus", &nclusSH );
    C->SetBranchAddress( "bb.ps.e", &ePS );
    C->SetBranchAddress( "bb.ps.rowblk", &rblkPS );
    C->SetBranchAddress( "bb.ps.colblk", &cblkPS );
    C->SetBranchAddress( "bb.ps.idblk", &idblkPS );
    C->SetBranchAddress( "bb.ps.atimeblk", &atimePS );

    // hodoscope cluster mean time
    Int_t Nhodotmean; 
    Double_t hodotmean[maxClus];
      
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "Ndata.bb.hodotdc.clus.tmean", 1 );
      
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", hodotmean );
    C->SetBranchAddress( "Ndata.bb.hodotdc.clus.tmean", &Nhodotmean );

    // track branches
    Double_t ntrack, p[maxTrack],px[maxTrack],py[maxTrack],pz[maxTrack],xtr[maxTrack],ytr[maxTrack],thtr[maxTrack],phtr[maxTrack];
    Double_t vx[maxTrack],vy[maxTrack],vz[maxTrack];
    Double_t xtgt[maxTrack],ytgt[maxTrack],thtgt[maxTrack],phtgt[maxTrack];

    C->SetBranchStatus( "bb.tr.n", 1 );
    C->SetBranchStatus( "bb.tr.p", 1 );
    C->SetBranchStatus( "bb.tr.px", 1 );
    C->SetBranchStatus( "bb.tr.py", 1 );
    C->SetBranchStatus( "bb.tr.pz", 1 );    
    C->SetBranchStatus( "bb.tr.x", 1 );    
    C->SetBranchStatus( "bb.tr.y", 1 );    
    C->SetBranchStatus( "bb.tr.th", 1 );    
    C->SetBranchStatus( "bb.tr.ph", 1 );    
    C->SetBranchStatus( "bb.tr.vx", 1 );
    C->SetBranchStatus( "bb.tr.vy", 1 );
    C->SetBranchStatus( "bb.tr.vz", 1 );
    C->SetBranchStatus( "bb.tr.tg_x", 1 );
    C->SetBranchStatus( "bb.tr.tg_y", 1 );
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );
    
    C->SetBranchAddress( "bb.tr.n", &ntrack );
    C->SetBranchAddress( "bb.tr.p", p );
    C->SetBranchAddress( "bb.tr.px", px );
    C->SetBranchAddress( "bb.tr.py", py );
    C->SetBranchAddress( "bb.tr.pz", pz );    
    C->SetBranchAddress( "bb.tr.x", xtr );    
    C->SetBranchAddress( "bb.tr.y", ytr );    
    C->SetBranchAddress( "bb.tr.th", thtr );    
    C->SetBranchAddress( "bb.tr.ph", phtr );    
    C->SetBranchAddress( "bb.tr.vx", vx );
    C->SetBranchAddress( "bb.tr.vy", vy );
    C->SetBranchAddress( "bb.tr.vz", vz );
    C->SetBranchAddress( "bb.tr.tg_x", xtgt );
    C->SetBranchAddress( "bb.tr.tg_y", ytgt );
    C->SetBranchAddress( "bb.tr.tg_th", thtgt );
    C->SetBranchAddress( "bb.tr.tg_ph", phtgt );

    // ekine branches
    Double_t ekineQ2, ekineW2, ekineeps, ekinenu, ekineqx, ekineqy, ekineqz;

    C->SetBranchStatus( "e.kine.Q2", 1 );
    C->SetBranchStatus( "e.kine.W2", 1 );
    C->SetBranchStatus( "e.kine.epsilon", 1 );
    C->SetBranchStatus( "e.kine.nu", 1 );
    C->SetBranchStatus( "e.kine.q_x", 1 );
    C->SetBranchStatus( "e.kine.q_y", 1 );
    C->SetBranchStatus( "e.kine.q_z", 1 );

    C->SetBranchAddress( "e.kine.Q2", &ekineQ2 );
    C->SetBranchAddress( "e.kine.W2", &ekineW2 );
    C->SetBranchAddress( "e.kine.epsilon", &ekineeps );
    C->SetBranchAddress( "e.kine.nu", &ekinenu );
    C->SetBranchAddress( "e.kine.q_x", &ekineqx );
    C->SetBranchAddress( "e.kine.q_y", &ekineqy );
    C->SetBranchAddress( "e.kine.q_z", &ekineqz );

    // fEvtHdr branches
    UInt_t rnum, gevnum, trigbits;

    C->SetBranchStatus( "fEvtHdr.fRun", 1 );
    C->SetBranchStatus( "fEvtHdr.fEvtNum", 1 );
    C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );

    C->SetBranchAddress( "fEvtHdr.fRun", &rnum );
    C->SetBranchAddress( "fEvtHdr.fEvtNum", &gevnum );
    C->SetBranchAddress( "fEvtHdr.fTrigBits", &trigbits );

    // other bb branches
    Double_t gemNhits, eop;

    C->SetBranchStatus( "bb.gem.track.nhits", 1 );
    C->SetBranchStatus( "bb.etot_over_p", 1 );

    C->SetBranchAddress( "bb.gem.track.nhits", &gemNhits );
    C->SetBranchAddress( "bb.etot_over_p", &eop );

    // define very wide globalcut
    TCut GCut = gcut.c_str();

    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    // get experimental quantities by run
    // set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
    Double_t BdL = hcal::sbsmaxfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbthr) * target_rho * target_dEdx;

    // set nucleon default for LH2
    std::string nucleon = "p"; 

    // event indices
    long nevent = 0, npassed = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    Double_t TSlastEv = -1.;
    Double_t TSnextEv = Sagg[0].enumber;
    Double_t TSlastt = Sagg[0].clock_diff;
    Double_t TSlastCur = Sagg[0].avg_current;
    Double_t TSlastC = chargeAvgLast; //Just use last average for first few events wherever necessary.
    Double_t TSstructidx=-1;

    while (C->GetEntry(nevent++)) {
	
      std::cout << "Processing run " << runnum << "(" << irun << " / " << nruns << "), event " << nevent << " / " << nevents << "\r";
      std::cout.flush();
	
      //Access TS charge struct members
      if( gevnum>=TSnextEv ){
	TSstructidx++;
	TSlastEv = Sagg[TSstructidx].enumber;
	TSnextEv = Sagg[TSstructidx+1].enumber;
	TSlastt = Sagg[TSstructidx].clock_diff;
	TSlastCur = Sagg[TSstructidx].avg_current;
	TSlastC = Sagg[TSstructidx].acc_charge;
      }

      ///////
      //Single-loop globalcut method. Save pass/fail for output tree.
      bool failedglobal = false;

      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
	//cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
      }
      failedglobal = GlobalCut->EvalInstance(0) == 0;
      if( failedglobal )
	continue;

      ///////
      //Physics calculations
      //correct beam energy with vertex information
      Double_t ebeam_c = ebeam - ( (vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );

      TVector3 vertex( 0., 0., vz[0] );

      //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t precon = p[0] + Eloss_outgoing;

      //set up four-momenta with some empty for various calculation methods
      TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
      TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect
      TLorentzVector ptarg; ptarg.SetPxPyPzE( 0., 0., 0., M_avg ); //target momentum
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      TLorentzVector pN; //N' momentum
      
      //simple calculations for e' and N'
      Double_t etheta = acos( pe.Pz() / pe.E() );
      Double_t ephi = atan2( pe.Py(), pe.Px() );
      Double_t pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
      Double_t phNexp = ephi + hcal::PI;

      //e' p reconstruction with track angles (not momentum)
      Double_t nu = pbeam.E() - pcent;
      Double_t pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
      Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
      TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
      pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
      Double_t Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
      Double_t W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;

      /////////////////////
      //W2 elastic cut
      bool failedW2 = W2>W2max || W2<W2min;
      if(failedW2)
	continue;

      npassed++;

      //Calculate h-arm quantities
      vector<Double_t> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1]; //from primary blk
      Double_t KE_p = nu; //For elastics total predicted by earm
      Double_t SFrac = hcale/KE_p; //Measured
      Double_t dx = hcalx - xyhcalexp[0];
      Double_t dy = hcaly - xyhcalexp[1];
      TVector3 neutdir = (hcalpos - vertex).Unit();

      // Proton deflection not yet reliable - perhaps never will be given non-uniformity of field
      // Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + hcal::sbsdipolegap/2.0) );
      // TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
      // Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );

      Double_t thetapq = acos( neutdir.Dot( qv.Unit() ) );
      Double_t eoverp = (ePS + eSH) / p[0];

      //Target Energy in HCal for calibrations
      Double_t KE_exp = KE_p;
      Double_t KE_exp_corr = KE_p*sffac;

      //charge 
      Double_t comp_ev_fraction = (Double_t)npassed/(Double_t)nevent;
      Double_t ev_fraction = (Double_t)npassed/(Double_t)nevents;
      Double_t accumulated_charge = charge*ev_fraction;
      
      Double_t cointime = hcalatime - atimeSH;

      //Fill diagnostic histos
      hQ2mag->Fill( mag, Q2 );
      hW2mag->Fill( mag, W2 );
      hdxmag->Fill( mag, dx );
      hdymag->Fill( mag, dy );

      //Fill new output tree     
      dx_out = dx;
      dy_out = dy;
      W2_out = W2;
      Q2_out = Q2;
      nu_out = nu;
      NT_out = KE_exp;
      NT_corr_out = KE_exp_corr;
      eoverp_out = eoverp;
      thetapq_out = thetapq;
      charge_ts_out = TSlastC;
      epcent_out = pcent;
      eprecon_out = precon;
      adc_coin_out = cointime;
      mag_out = mag;
      run_out = runnum;
      tar_out = t;
      event_out = (Int_t)gevnum;
      trig_out = (Int_t)trigbits;
      charge_avg_out = charge;
      comp_event_frac_out = comp_ev_fraction;
      event_frac_out = ev_fraction;
	
      //Fill old output tree
      bb_tr_vx_out = vx[0];
      bb_tr_vy_out = vy[0];
      bb_tr_vz_out = vz[0];
      bb_tr_p_out = p[0];
      bb_tr_px_out = px[0];
      bb_tr_py_out = py[0];
      bb_tr_pz_out = pz[0];
      bb_tr_n_out = ntrack;
      e_kine_Q2_out = ekineQ2;
      e_kine_W2_out = ekineW2;
      e_kine_nu_out = ekinenu;
      bb_ps_e_out = ePS;
      bb_ps_rowblk_out = rblkPS;
      bb_ps_colblk_out = cblkPS;
      bb_sh_e_out = eSH;
      bb_sh_rowblk_out = rblkSH;
      bb_sh_colblk_out = cblkSH;
      bb_sh_atimeblk_out = atimeSH;
      bb_sh_nclus_out = nclusSH;
      bb_hodotdc_clus_tmean_out = hodotmean[0];
      bb_gem_track_nhits_out = gemNhits;
      bb_etot_over_p_out = eop;
      sbs_hcal_e_out = hcale;
      sbs_hcal_x_out = hcalx;
      sbs_hcal_y_out = hcaly;
      sbs_hcal_rowblk_out = hcalr;
      sbs_hcal_colblk_out = hcalc;
      sbs_hcal_atimeblk_out = hcalatime;
      sbs_hcal_tdctimeblk_out = hcaltdc;
      sbs_hcal_nclus_out = nclus;
      Ndata_sbs_hcal_clus_id_out = Nhcalcid;
      for( Int_t c=0; c<Nhcalcid; ++c ){
	Double_t cid = hcalcid[c];
	Double_t ce = hcalce[c];
	Double_t cx = hcalcx[c];
	Double_t cy = hcalcy[c];
	Double_t ctdc = hcalctdctime[c];
	Double_t catime = hcalcatime[c];

	sbs_hcal_clus_id_out[c] = cid;
	sbs_hcal_clus_e_out[c] = ce;
	sbs_hcal_clus_x_out[c] = cx;	  
	sbs_hcal_clus_y_out[c] = cy;
	sbs_hcal_clus_tdctime_out[c] = ctdc;
	sbs_hcal_clus_atime_out[c] = catime;
      }
      sbs_hcal_nblk_out = nblk;
      Ndata_sbs_hcal_clus_blk_id_out = Nhcalcbid;
      for( Int_t b=0; b<Nhcalcbid; ++b ){
	Double_t cbid = hcalcbid[b];
	Double_t cbe = hcalcbe[b];
	Double_t cbx = hcalcbx[b];
	Double_t cby = hcalcby[b];
	Double_t cbtdc = hcalcbtdctime[b];
	Double_t cbatime = hcalcbatime[b];

	sbs_hcal_clus_blk_id_out[b] = cbid;
	sbs_hcal_clus_blk_e_out[b] = cbe;
	sbs_hcal_clus_blk_x_out[b] = cbx;	  
	sbs_hcal_clus_blk_y_out[b] = cby;
	sbs_hcal_clus_blk_tdctime_out[b] = cbtdc;
	sbs_hcal_clus_blk_atime_out[b] = cbatime;

      }

      P->Fill();

    }//end event loop

    // getting ready for the next run
    C->Reset();

  }//end run loop


  fout->Write();

  std::cout << std::endl << "Diagnostic complete. Output written to " << out_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}
