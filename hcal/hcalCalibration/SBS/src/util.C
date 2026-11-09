#include "../include/util.h"
#include "../include/hcal.h"

namespace util {

  ////Console
  // returns string with current month day and year (date)
  string getDate(){
    time_t now = time(0);
    tm ltm = *localtime(&now);
  
    string yyyy = to_string(1900 + ltm.tm_year);
    string mm = to_string(1 + ltm.tm_mon);
    string dd = to_string(ltm.tm_mday);
    string date = mm + '_' + dd + '_' + yyyy;
  
    return date;
  }

  // reads a set of parameters into a vector from a txt file at <const_path>. Assumes endl delimiter
  void readParam( std::string const_path, vector<Double_t> &param ) {

    ifstream const_file; const_file.open( const_path );

    //Get all params
    if(const_file.is_open()){
  
      Int_t n1=0;
      Double_t d1;
      string readline;
    
      while( getline( const_file, readline ) ){
	if( readline.at(0) == '#' ) {
	  continue;
	}
	istringstream iss( readline );
	iss >> d1;     
	param.push_back(d1);
	n1++;
      }
      
    }else{
      cerr << "Error on [util::readParam] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
    const_file.close();
  }

  // Counts the number of calibration sets for a given calibration parameter
  Int_t countSets( std::string const_path,
		   std::string type ){
    
    Int_t Nsets = 0;

    ifstream const_file; const_file.open( const_path.c_str() );

    //Get all params at datatype
    if(const_file.is_open()){

      string readline;

      while( getline( const_file, readline ) ){
	
	TString Tline = (TString)readline;
	
	if( Tline.BeginsWith(type) )
	  Nsets++;

      }

      return Nsets;

    }else{
      cerr << "Error on [util::countSets] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  } 
  

  // Reads database-style file for 288 (hcal::maxHCalChan) elements after a given moniker (type)
  // A return character must exist after the list of parameters or error will be thrown
  void readDB( std::string const_path, 
	       std::string timestamp, 
	       std::string type,
	       Double_t param[] ) {

    ifstream const_file; const_file.open( const_path.c_str() );

    //Get all params at datatype
    if(const_file.is_open()){

      Int_t n0=0;
      Double_t d0;
      std::string readline;
      bool read_param = false;
      bool found_tstamp = false;

      if( timestamp.compare("none")==0 ){
	cout<<"here none == 0"<<endl;
	found_tstamp = true;
      }

      while( getline( const_file, readline ) ){
	//cout<<"got the line, in the while loop"<<endl;
	if( n0==hcal::maxHCalChan ){
	  cout<<"limit reached"<<endl;	  
	  const_file.close();
	  return;
	}

	TString Tline = (TString)readline;
  
	if( Tline.BeginsWith(timestamp) && !found_tstamp ){
	  cout<<"time stamp found"<<endl;
	  found_tstamp = true;
	  continue;
	}

	if( Tline.BeginsWith(type) && !read_param && found_tstamp ){
	  //cout<<"here"<<endl;
	  read_param = true;
	  continue;
	}
      
	if( found_tstamp && read_param ){

	  istringstream iss( readline );
	  while( iss >> d0 ){
	    
	    param[n0] = d0;

	    n0++;
	    //std::cout << "N values: " << n0 << "/" << hcal::maxHCalChan << ". Value: " << d0 << std::endl;
	  }
	}
      }
      
      cerr << "Error on [util::readDB] Did not find required number of elements in param[] at " << const_path << endl;
      throw;

    }else{
      cerr << "Error on [util::readDB] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  }

  // Overload for tdc calibration constant (1 element, same line)
  void readDB( std::string const_path, 
	       std::string timestamp, 
	       std::string type,
	       Double_t &param ) {

    ifstream const_file; const_file.open( const_path.c_str() );
    
    //Get all params at datatype
    if(const_file.is_open()){

      Double_t d0;
      string readline;
      bool found_tstamp = false;

      if( timestamp.compare("none")==0 )
	found_tstamp = true;

      while( getline( const_file, readline ) ){

	TString Tline = (TString)readline;
	
	if( Tline.BeginsWith(timestamp) && !found_tstamp ){	  
	  found_tstamp = true;
	  continue;
	}

	if( found_tstamp && Tline.BeginsWith(type) ){
	  string temp;
	  bool val = false;
	  for( Int_t l=0; l<readline.length(); l++){
	    if( val==true )
	      temp+=readline[l];
	    if( readline[l] == '=' ) 
	      val=true;
	  }
	  param = stof( temp );
	  const_file.close();
	  return;
	}
      }
      cerr << "Error on [util::readDB] timestamp/type does not exist in file at " << const_path << endl;
      throw;

    }else{
      cerr << "Error on [util::readDB] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  }


  // Overload for other constants (N elements, same line)
  void readDB( std::string const_path, 
	       std::string timestamp, 
	       std::string type,
	       vector<Double_t> &param ) {

    ifstream const_file; const_file.open( const_path.c_str() );
    
    //Get all params at datatype
    if(const_file.is_open()){

      Double_t d0;
      string readline;
      bool found_tstamp = false;

      if( timestamp.compare("none")==0 )
	found_tstamp = true;

      while( getline( const_file, readline ) ){

	TString Tline = (TString)readline;
	
	if( Tline.BeginsWith(timestamp) && !found_tstamp ){	  
	  found_tstamp = true;
	  continue;
	}

	if( found_tstamp && Tline.BeginsWith(type) ){
	  string temp;
	  bool val = false;
	  for( Int_t l=0; l<readline.length(); l++){
	    if( val==true )
	      temp+=readline[l];
	    if( readline[l] == '=' ) 
	      val=true;
	  }
	  istringstream iss( temp );
	  while( iss >> d0 )
	    param.push_back(d0);

	  const_file.close();
	  return;
	}
      }
      cerr << "Error on [util::readDB] timestamp/type does not exist in file at " << const_path << endl;
      throw;

    }else{
      cerr << "Error on [util::readDB] Parameter file doesn't exist at " << const_path << endl;
      throw;
    }
  }

  // compares two timestamps and returns the one that is later in time. Format must be similar
  void tsCompare( std::string tsA, 
  		  std::string tsB, 
  		  std::string &tsLate ){
    
    std::string tempA;
    std::string tempB;

    if( tsA.compare(tsB)==0 ){ //if they are the same, return tsA arb.
      tsLate=tsA;
      return;
    }else if(tsA.compare("none")==0){ //if they are not, check if A is none, then return the other
      tsLate=tsB;
      return;
    }else if(tsB.compare("none")==0){ //vice versa
      tsLate=tsA;
      return;
    }else{ //else, get new strings deleting "-------[ " preamble as it can be of different length
      bool Atoi = false;
      for( Int_t l=0; l<tsA.size(); l++ ){
	char a = tsA[l];
	if( Atoi )
	  tempA+=a;
	if( a==' ' )
	  Atoi = true;
      }
      bool Btoi = false;
      for( Int_t l=0; l<tsB.size(); l++ ){
	char b = tsB[l];
	if( Btoi )
	  tempB+=b;
	if( b==' ' )
	  Btoi = true;
      }
    
      Int_t l=0;
      while( tempA[l++] ){ //compare the two parsed strings and return newest
	
    	char a = tempA[l];
    	char b = tempB[l]; 

	if( (Int_t)a>(Int_t)b ){
	  tsLate = tsA;
	  return;
	}else if( (Int_t)a<(Int_t)b ){
	  tsLate = tsB;
	  return;
	}else
	  continue;
      }

    }

    cerr << "Erroron [util::tsCompare] String comparison failed" << endl;
    throw; //comparison should be exhaustive

  }
  
  //Reads run list .csv file - may wish to build in selections and exclusions in an overload
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   std::string experiment,     // experiment {gmn,gen,genrp,gep,etc.}
		   Int_t &nruns,               // No. of runs to analyze
		   Int_t sbsconf,              // SBS configuration
		   Int_t replay_pass,          // replay pass
		   Int_t verbose,              // verbosity
		   vector<calrun> &crun)       // Output: Vector of crun structs
  {
    // Define the name of the relevant run spreadsheet
    if ( experiment.compare("gmn")==0 && replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = Form("%s%sruns_pass%d.csv",runsheet_dir.c_str(),experiment.c_str(),replay_pass);

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;     
      string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;

	  //cout << temptoken << endl;

	  temp.push_back(temptoken);
	}
	// add relevant info to calrun objects
	if (stoi(temp[0]) == sbsconf ) {
	  if (crun.size() >= nruns) break;
	  calrun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  crun.push_back(temp_cr);
	}

	temp.clear();
      }
      //Update nruns with total no. of runs to analyze
      nruns = crun.size();
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << crun[0];
	std::cout << "Last run info:" << std::endl << crun[nruns-1];
      }
    }else{
      cerr << "Error on [util::ReadRunList], run spreadsheet doesn't exist" << endl;
      throw;
    }
    run_data.close();
  }

  //Read cut list .csv file. sbsconf + target + field MUST uniquely specify cuts loaded here
  void ReadCutList(std::string cutsheet_dir,  // Dir. path containing CSV files with cut info
		   std::string experiment,     // experiment {gmn,gen,genrp,gep,etc.}
		   Int_t sbsconfig,            // SBS configuration
		   Int_t calib_set,            // Calibration set
		   Int_t replay_pass,          // replay pass
		   std::string target,         // target
		   Int_t field,                // sbs field in percent
		   Int_t verbose,              // verbosity
		   vector<calcut> &ccut)       // Output: Vector of ccut structs
  {
    // Define the name of the relevant cut spreadsheet
    if ( experiment.compare("gmn")==0 && replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string cut_spreadsheet = Form("%s%scuts_pass%d.csv",cutsheet_dir.c_str(),experiment.c_str(),replay_pass);

    // Reading the spreadsheet
    ifstream cut_data; cut_data.open(cut_spreadsheet);
    string readline;
    bool max_pushback = false;
    if(cut_data.is_open()){  
      std::cout << "Reading cut info from: "<< cut_spreadsheet 
		<< std::endl << std::endl;   
      string skip_header; getline(cut_data, skip_header); // skipping column header
      while(getline(cut_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  
	  //cout << temptoken << endl;

	  temp.push_back(temptoken);
	}
	// select cut based on config, target, and field
	if (stoi(temp[0]) == sbsconfig &&
	    stoi(temp[1]) == calib_set &&
	    temp[2].compare(target)==0 &&
	    stoi(temp[3]) == field ) {
	  if( max_pushback ) //throw error if more than one element added to cut array
	    throw "Error on [util::ReadCutList], config/target/field not uniquely specified in cut spreadsheet";
	  calcut temp_cc;
	  temp_cc.SetDataCutSheet(temp);
	  ccut.push_back(temp_cc);
	  max_pushback = true;
	  cout<<"got to the loop"<<endl;
	}

	temp.clear();
      }
      std::cout<<"DEBUG"<<endl;
      //Update nruns with total no. of runs to analyze
      if (verbose == 1) {
	std::cout << "First cut info:" << std::endl << ccut[0];
	std::cout << "Last cut info:" << std::endl << ccut.back();
      }
    }else{
      cerr << "Error on [util::ReadCutList], cut spreadsheet doesn't exist" << endl;
      throw;
    }
    cut_data.close();
    std::cout<<"DEBUG 1"<<endl;
  }


  //Read diagnostic cut list .csv file. sbsconf + target + field MUST uniquely specify cuts loaded here
  void ReadDiagnosticCutList(std::string cutsheet_dir,  // Dir. path containing CSV files with cut info
			     std::string experiment,     // experiment {gmn,gen,genrp,gep,etc.}
			     Int_t sbsconfig,            // SBS configuration
			     std::string target,         // target
			     Int_t field,                // sbs field in percent
			     Int_t verbose,              // verbosity
			     vector<caldiag> &cdcut)       // Output: Vector of cdcut structs
  {
    // Define the name of the relevant cut spreadsheet
    std::string cut_spreadsheet = Form("%s%scuts_diagnostic.csv",cutsheet_dir.c_str(),experiment.c_str());

    // Reading the spreadsheet
    ifstream cut_data; cut_data.open(cut_spreadsheet);
    string readline;
    bool max_pushback = false;
    if(cut_data.is_open()){   
      std::cout << "Reading diagnostic cut info from: "<< cut_spreadsheet 
		<< std::endl << std::endl;
      string skip_header; getline(cut_data, skip_header); // skipping column header
      while(getline(cut_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;

	  //cout << temptoken << endl;

	  temp.push_back(temptoken);
	}
	// select cut based on config, target, and field
	if (stoi(temp[0]) == sbsconfig &&
	    temp[1].compare(target)==0 &&
	    stoi(temp[2]) == field ) {
	  if( max_pushback ) //throw error if more than one element added to cut array
	    throw "Error on [util::ReadCutList], config/target/field not uniquely specified in cut spreadsheet";
	  caldiag temp_cc;
	  temp_cc.SetDataDiagnosticCutSheet(temp);
	  cdcut.push_back(temp_cc);
	  max_pushback = true;
	}

	temp.clear();
      }
      //Update nruns with total no. of runs to analyze
      if (verbose == 1) {
	std::cout << "First cut info:" << std::endl << cdcut[0];
	std::cout << "Last cut info:" << std::endl << cdcut.back();
      }
    }else{
      cerr << "Error on [util::ReadCutList], cut spreadsheet doesn't exist" << endl;
      throw;
    }
    cut_data.close();
  }

  //takes diagnostic and supplemental structs and creates many cut and parameter reports by targ/field
  void diagnosticReport(const std::vector<caldiag>& data, const suppset& supplement,std::string report_path = "report_out.root") {
    // Create a new ROOT file to save the canvases
    TFile *file = new TFile(report_path.c_str(), "RECREATE");

    for (size_t i = 0; i < data.size(); ++i) {

      //Assuming standard sizes, these numbers produce an output with reasonable margins
      TCanvas *canvas = new TCanvas(Form("c%zu", i+1), Form("Canvas %zu", i+1), 1000, 900);
      canvas->cd();
      double yPos = 0.9; // Starting Y position for text, top of the canvas
      double yStep = 0.025; // Step size for each line of text
      double textSize = 0.02; //size of text
      double textSize_h = 0.025; //size of header text

      std::stringstream ss;
      ss << "General Set ADCt Alignment Info";
      TLatex *genLatex_l1 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l1->SetTextColor(kRed-5);
      genLatex_l1->SetTextSize(textSize_h);
      genLatex_l1->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Experiment: " << supplement.exper 
	 << ", Configuration: " << supplement.kine 
	 << ", Pass: " << supplement.pass;
      TLatex *headerLatex_l1 = new TLatex(0.1, yPos, ss.str().c_str());
      headerLatex_l1->SetTextColor(kBlack);
      headerLatex_l1->SetTextSize(textSize);
      headerLatex_l1->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Creation Date: " << supplement.date;
      TLatex *headerLatex_l2 = new TLatex(0.1, yPos, ss.str().c_str());
      headerLatex_l2->SetTextColor(kBlack);
      headerLatex_l2->SetTextSize(textSize);
      headerLatex_l2->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Run Range: " << supplement.runb << " - " << supplement.rune;
      TLatex *headerLatex_l3 = new TLatex(0.1, yPos, ss.str().c_str());
      headerLatex_l3->SetTextColor(kBlack);
      headerLatex_l3->SetTextSize(textSize);
      headerLatex_l3->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Exclusion Range: " << supplement.runexb << " - " << supplement.runexe;
      TLatex *headerLatex_l4 = new TLatex(0.1, yPos, ss.str().c_str());
      headerLatex_l4->SetTextColor(kBlack);
      headerLatex_l4->SetTextSize(textSize);
      headerLatex_l4->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Target(s) Used: " << supplement.targ;
      TLatex *headerLatex_l5 = new TLatex(0.1, yPos, ss.str().c_str());
      headerLatex_l5->SetTextColor(kBlack);
      headerLatex_l5->SetTextSize(textSize);
      headerLatex_l5->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "";
      TLatex *genLatex_l2 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l2->SetTextColor(kRed-5);
      genLatex_l2->SetTextSize(textSize);
      genLatex_l2->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Electron Arm Elastic Cuts";
      TLatex *genLatex_l3 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l3->SetTextColor(kRed-5);
      genLatex_l3->SetTextSize(textSize_h);
      genLatex_l3->Draw();
      yPos -= yStep;

      // Now display the data from the caldiag struct
      ss.str("");
      ss << "Target: " << data[i].target;
      TLatex *caldiagLatex_la = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_la->SetTextColor(kBlue);
      caldiagLatex_la->SetTextSize(textSize);
      caldiagLatex_la->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "SBS Field: " << data[i].field;
      TLatex *caldiagLatex_lb = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_lb->SetTextColor(kBlue);
      caldiagLatex_lb->SetTextSize(textSize);
      caldiagLatex_lb->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Global Elastic Cuts: " << data[i].gcut;
      TLatex *caldiagLatex_l1 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l1->SetTextColor(kBlue);
      caldiagLatex_l1->SetTextSize(textSize);
      caldiagLatex_l1->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "W2 Min: " << data[i].W2_min;
      TLatex *caldiagLatex_l2 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l2->SetTextColor(kBlue);
      caldiagLatex_l2->SetTextSize(textSize);
      caldiagLatex_l2->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "W2 Max: " << data[i].W2_max;
      TLatex *caldiagLatex_l3 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l3->SetTextColor(kBlue);
      caldiagLatex_l3->SetTextSize(textSize);
      caldiagLatex_l3->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "";
      TLatex *genLatex_l4 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l4->SetTextColor(kRed-5);
      genLatex_l4->SetTextSize(textSize);
      genLatex_l4->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Hadron Arm Elastic Cuts (" << supplement.spotsig << " sigma)";
      TLatex *genLatex_l5 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l5->SetTextColor(kRed-5);
      genLatex_l5->SetTextSize(textSize_h);
      genLatex_l5->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "dy mean: " << data[i].dy0;
      TLatex *caldiagLatex_l4 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l4->SetTextColor(kBlue);
      caldiagLatex_l4->SetTextSize(textSize);
      caldiagLatex_l4->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "dy sigma: " << data[i].dy_sig;
      TLatex *caldiagLatex_l5 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l5->SetTextColor(kBlue);
      caldiagLatex_l5->SetTextSize(textSize);
      caldiagLatex_l5->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "dx proton mean: " << data[i].dx0_p;
      TLatex *caldiagLatex_l6 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l6->SetTextColor(kBlue);
      caldiagLatex_l6->SetTextSize(textSize);
      caldiagLatex_l6->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "dx proton sigma: " << data[i].dx_sig_p;
      TLatex *caldiagLatex_l7 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l7->SetTextColor(kBlue);
      caldiagLatex_l7->SetTextSize(textSize);
      caldiagLatex_l7->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "dx neutron mean: " << data[i].dx0_n;
      TLatex *caldiagLatex_l8 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l8->SetTextColor(kBlue);
      caldiagLatex_l8->SetTextSize(textSize);
      caldiagLatex_l8->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "dx neutron sigma: " << data[i].dx_sig_n;
      TLatex *caldiagLatex_l9 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l9->SetTextColor(kBlue);
      caldiagLatex_l9->SetTextSize(textSize);
      caldiagLatex_l9->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "";
      TLatex *genLatex_l6 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l6->SetTextColor(kRed-5);
      genLatex_l6->SetTextSize(textSize);
      genLatex_l6->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Other Cuts";
      TLatex *genLatex_l7 = new TLatex(0.1, yPos, ss.str().c_str());
      genLatex_l7->SetTextColor(kRed-5);
      genLatex_l7->SetTextSize(textSize);
      genLatex_l7->Draw();
      yPos -= yStep;

      ss.str("");
      ss << "Minimum Events Per Cell: " << data[i].min_ev;
      TLatex *caldiagLatex_l10 = new TLatex(0.1, yPos, ss.str().c_str());
      caldiagLatex_l10->SetTextColor(kBlue);
      caldiagLatex_l10->SetTextSize(textSize);
      caldiagLatex_l10->Draw();
      yPos -= yStep;

      // Write the canvas to the file
      canvas->Write();
      canvas->Update();
    }

    // Close the file
    file->Close();
  }

  //Geometry/physics
  //sets hcal axes by kinematic (from hall coordinates)
  void sethcalaxes( Double_t sbstheta_rad,                           // SBS angle in radian 
		    vector<TVector3> &hcal_axes ) {
    TVector3 hcal_zaxis( sin(-sbstheta_rad), 0, cos(-sbstheta_rad) ); // Clock-wise rotation about Y axis
    TVector3 hcal_xaxis( 0, -1, 0 );                                  // -Y axis of Hall CoS = X axis of hcal CoS
    TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
    hcal_axes.push_back(hcal_xaxis);
    hcal_axes.push_back(hcal_yaxis);
    hcal_axes.push_back(hcal_zaxis);
  }

  //gets expected location of scattered nucleon assuming straight line projections from BB track
  void getxyhcalexpect( TVector3 vertex, TVector3 pNhat, TVector3 hcal_origin, 
			vector<TVector3> hcal_axes, vector<Double_t> &xyhcalexpect ) {
    // intersection of a ray with a plane in hcal coordinates
    Double_t sintersect = ( hcal_origin - vertex).Dot(hcal_axes[2] ) / ( pNhat.Dot(hcal_axes[2]) );
    // ray from Hall origin onto the face of hcal where the nucleon hit
    TVector3 hcal_intersect = vertex + sintersect*pNhat; 

    Double_t xexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[0] );
    Double_t yexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[1] );

    xyhcalexpect.push_back( xexpect_hcal );
    xyhcalexpect.push_back( yexpect_hcal );
  }
  
  //checks particle id with information from observed dx distributions.
  void checkPID( std::string target, Double_t dx0p, Double_t dx0n, Double_t dxsigp, 
		 Double_t dxsign, Double_t dx, Int_t &pid ){
    bool isproton = abs(dx-dx0p)<3*dxsigp;
    bool isneutron = abs(dx-dx0n)<3*dxsign;
    bool isambiguous = isproton && isneutron;
    
    if( isproton ) 
      pid = 1;
    else if( isneutron )
      pid = 2;
    else if( isambiguous ){
      if( target.compare("lh2")==0 || target.compare("LH2") ) //probably not useful since isproton comes first
	pid = 1;
      else
	pid = 0;
    }
    else
      pid = -1; //neither
  }

  // checks if a point is within proton or neutron spot. rotation angle in rad, if ever applicable
  bool Nspotcheck(Double_t dy, Double_t dx, Double_t dy_mean, Double_t dx_mean, Double_t dy_sigma, Double_t dx_sigma, Double_t rotationAngle=0) {

    // Caculate semimajor and semiminor axes
    Double_t dyAxis = dy_sigma;
    Double_t dxAxis = dx_sigma;

    // Apply rotation angle if applicable
    Double_t cosAngle = std::cos(rotationAngle);
    Double_t sinAngle = std::sin(rotationAngle);
    Double_t dyRot = (dy - dy_mean) * cosAngle + (dx - dx_mean) * sinAngle;
    Double_t dxRot = (dx - dx_mean) * cosAngle - (dy - dy_mean) * sinAngle;

    // Check if point is within the ellipse equation
    Double_t result = ((dyRot * dyRot) / (dyAxis * dyAxis)) + ((dxRot * dxRot) / (dxAxis * dxAxis));
    return result <= 1.0;
  }

  //assign score to potential cluster based on probability density of gaussian fit to atime and maximum energy, exclude position info. Here, gfit is a vector with coin atime fit parameters, val_2 is the cluster coin atime, val_1 is the cluster energy, and max1 is the max cluster energy for the event
  Double_t assignScore( double val_1, double val_2, double max1, const std::vector<double> &gfit ) {

    if( gfit.size()!=3 ){
      cout << "ERROR: size of dxfit vector not equal to expected number of gaussian parameters (3)." << endl;
      return 0.;
    }

    // Build a TF1 with fit provided fit parameters
    TF1 *gauss = new TF1("gauss", "gaus", 0, 100);
    gauss->SetParameter(0,gfit[0]);
    gauss->SetParameter(1,gfit[1]);
    gauss->SetParameter(2,gfit[2]);
    
    Double_t x_value = val_2;
    Double_t density = gauss->Eval(x_value);

    // Compute the score based on val 2
    Double_t score = density / gauss->GetParameter(0);  // Normalize by the peak of the Gaussian

    // Include a product component based on val 1 (linear weight)
    score *= val_1 / max1;

    if(score==0){
      //cout << "WARNING: score is zero. val_1=" << val_1 << ", val_2= " << val_2 << " density= " << density << " 1 comp= " << val_1 / max1 << endl;
      score=1e-38; //write very small number to score to exclude it without breaking the sorting later
    }

    delete gauss; //prevent memory leak

    return score;
  }
  
  //Histogram utility functions
  //Gets the maximum bin on a user passed range
  Int_t get_max_bin_on_range( TH1D *h, Double_t xlow, Double_t xhigh ){

    Int_t bin1 = h->FindBin(xlow);  // convert xlow to bin number
    Int_t bin2 = h->FindBin(xhigh); // convert xhigh to bin number

    // Initial maximum and its position
    Double_t maxVal = h->GetBinContent(bin1);
    Int_t maxBin = bin1;

    // Loop over bins in the range to find the maximum
    for (Int_t i = bin1; i <= bin2; i++) {
      Double_t binContent = h->GetBinContent(i);
      if (binContent > maxVal) {
	maxVal = binContent;
	maxBin = i;
      }
    }

    return maxBin;
  }

  //Just gets max bin ignoring overflow and underflow
  double get_max_bin_content_ignore_overflow(TH1D* hist) {
    double maxContent = 0;
    // Start from bin 1 to ignore underflow, end at GetNbinsX() to ignore overflow
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
      double binContent = hist->GetBinContent(i);
      if (binContent > maxContent) {
	maxContent = binContent;
      }
    }
    return maxContent;
  }

  //Fits a gaussian to a distribution twice, first course, then fine
  std::vector<Double_t> fitGaussianAndGetFineParams(TH1D* hist, Double_t sig, Double_t low = -1e38, Double_t high = 1e38) {
    
    std::vector<Double_t> params(3); // Vector to store amplitude, mean, and sigma

    if (!hist) {
      std::cerr << "Histogram is null!" << std::endl;
      return params;
    }

    // Find the bin numbers corresponding to the specified range. If default values are passed, real edge bins returned by FindBin().
    Int_t binLow = std::max(hist->FindBin(low),1);
    Int_t binHigh = std::min(hist->FindBin(high),hist->GetNbinsX());

    // Initial values for maximum content and bin
    Double_t maxContent = 0;
    Int_t maxBin = -1;

    // Iterate over the bins in the range
    for (Int_t i = binLow; i <= binHigh; ++i) {
      Double_t content = hist->GetBinContent(i);
      if (content > maxContent) {
	maxContent = content;
	maxBin = i;
      }
    }

    Double_t xMax = hist->GetXaxis()->GetBinCenter(maxBin);

    // Define the fit range
    Double_t fitMin = xMax - sig;
    Double_t fitMax = xMax + sig;

    // Fit the histogram within the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
    hist->Fit(gausFit, "RQ"); // "R" for fit range, "Q" for quiet mode (no print)

    // Get the mean from the fit
    Double_t amplitude = gausFit->GetParameter(0); // Parameter 0 is the amplitude of the Gaussian 
    Double_t mean = gausFit->GetParameter(1); // Parameter 1 is the mean of the Gaussian
    Double_t sigma = gausFit->GetParameter(2); // Parameter 2 is the std dev of the Gaussian

    // Clean up
    delete gausFit;

    // Fit the histogram within the fine range
    Double_t fit2Min = mean - 1*sigma;
    Double_t fit2Max = mean + 1*sigma;
    TF1 *gausFit_fine = new TF1("gausFit_fine", "gaus", fit2Min, fit2Max);
    gausFit_fine->SetParameters(amplitude,mean,sigma);
    hist->Fit(gausFit_fine, "RQ"); // "R" for fit range, "Q" for quiet mode (no print)

    // Store the parameters in the vector
    params[0] = gausFit_fine->GetParameter(0); // Amplitude
    params[1] = gausFit_fine->GetParameter(1); // Mean
    params[2] = gausFit_fine->GetParameter(2); // Sigma

    return params;
  }

  //function to get gaussian fit mean and std dev per x-bin for later tgraph plotting
  void sliceHisto( TH2D *h2, Int_t Nslices, Double_t fwhm, Int_t min_ev, vector<Double_t> &cell, vector<Double_t> &mean, vector<Double_t> &err, Double_t meanXmin = -1e38, Double_t meanXmax = 1e38, Double_t noise_floor = 1e38 ){
  
    TH1D *cellslice[Nslices];

    for( Int_t i=0; i<Nslices; i++ ){
    
      Double_t cellval = (Double_t)i+0.5;
      Double_t meanval = 0.;
      Double_t errval = 0.;

      cellslice[i] = h2->ProjectionY(Form("cellslice_%d",i+1),i+1,i+1);

      // Calculate the integral (total number of entries) in the range
      double sliceN = cellslice[i]->Integral(1, cellslice[i]->GetNbinsX());
      if( sliceN<min_ev ){ //continue if too sparse to get distribution
	//cout << "WARNING: Insufficient statistics on hist " << h2->GetName() << " cell " << cellval << " with " << sliceN << endl;
	continue;
      }

      //get reasonable limits on fit
      Double_t min = cellslice[i]->GetXaxis()->GetXmin();
      Double_t max = cellslice[i]->GetXaxis()->GetXmax();
      Int_t totalbins = cellslice[i]->GetNbinsX();
      Int_t binmax = cellslice[i]->GetMaximumBin();
      //Double_t maxcontent = cellslice[i]->GetMaximum();
      Double_t maxcontent = get_max_bin_content_ignore_overflow(cellslice[i]);

      if( maxcontent < noise_floor ) //catch, underflow/overflow may be off range
	continue;

      Double_t binmaxX = min+binmax*(max-min)/totalbins;
      Double_t llim = binmaxX-fwhm;
      Double_t ulim = binmaxX+fwhm;
      if( llim < min ){ //if the lower limit is below minimum, recalculate
	binmax = get_max_bin_on_range( cellslice[i], min, max );
	binmaxX = min+binmax*(max-min)/totalbins;
	llim = min;
	ulim = binmaxX + ( binmaxX - llim );
      }
    
      TF1 *fit1;
      cellslice[i]->Fit("gaus","Q","",llim,ulim);
      fit1 = cellslice[i]->GetFunction("gaus");
      meanval = fit1->GetParameter(1);
      errval = fit1->GetParameter(2);

      if( meanval > meanXmax || meanval < meanXmin ) //hard catch, if mean value out of known range use args
	continue;

      cell.push_back(cellval);
      mean.push_back(meanval);
      err.push_back(errval);
    }
  }

  void sliceHistoOnCanvas(TH2D *h2, Int_t Nslices, Double_t fwhm, Int_t min_ev, vector<Double_t> &cell, vector<Double_t> &mean, vector<Double_t> &err, TCanvas *c, Double_t meanXmin = -1e38, Double_t meanXmax = 1e38, Double_t noise_floor = 1e38) {
    
    TH1D *cellsliceC[Nslices];
    
    c->Divide(sqrt(Nslices)+1, sqrt(Nslices)+1); // Arrange the canvas in a grid close to square

    for (Int_t i = 0; i < Nslices; i++) {
      Double_t cellval = (Double_t)i + 0.5;
      Double_t meanval = 0.;
      Double_t errval = 0.;

      cellsliceC[i] = h2->ProjectionY(Form("cellsliceC_%d", i + 1), i + 1, i + 1);

      // Calculate the integral (total number of entries) in the range
      double sliceN = cellsliceC[i]->Integral(1, cellsliceC[i]->GetNbinsX());

      // Get maximum value for scaling and cuts
      Double_t histmaxa = cellsliceC[i]->GetMaximum(); 
      cellsliceC[i]->GetYaxis()->SetRangeUser(0,histmaxa*1.1);

      c->cd(i+1); // Go to the sub-pad
      cellsliceC[i]->Draw();

      // Create a legend
      TLegend *legR = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
      legR->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
      // Draw the legend
      legR->Draw();

      if (sliceN < min_ev) { // Continue if too sparse to get distribution
	//cout << "WARNING: Insufficient statistics on hist " << h2->GetName() << " cell " << cellval << " with " << sliceN << endl;
	cellsliceC[i]->SetLineColor(kRed);
	cellsliceC[i]->SetTitle(Form("i%d FAIL ev%0.2f",i,sliceN));
	c->Update();
	continue;
      }

      Double_t min = cellsliceC[i]->GetXaxis()->GetXmin();
      Double_t max = cellsliceC[i]->GetXaxis()->GetXmax();
      Int_t totalbins = cellsliceC[i]->GetNbinsX();
      Int_t binmax = cellsliceC[i]->GetMaximumBin();
      //Double_t maxcontent = cellsliceC[i]->GetMaximum();
      Double_t maxcontent = get_max_bin_content_ignore_overflow(cellsliceC[i]);

      if ( maxcontent < noise_floor ){ // Catch, underflow/overflow may be off range
	cellsliceC[i]->SetLineColor(kMagenta);
	cellsliceC[i]->SetTitle(Form("i%d FAIL noise_floor%0.2f",i,maxcontent));
	continue;
      }

      Double_t binmaxX = min + binmax * (max - min) / totalbins;
      Double_t llim = binmaxX - fwhm;
      Double_t ulim = binmaxX + fwhm;
      if (llim < min) {
	binmax = get_max_bin_on_range(cellsliceC[i], min, max);
	binmaxX = min + binmax * (max - min) / totalbins;
	llim = min;
	ulim = binmaxX + (binmaxX - llim);
      }

      cellsliceC[i]->Fit("gaus", "Q", "", llim, ulim);
      TF1 *fit1 = cellsliceC[i]->GetFunction("gaus");
      meanval = fit1->GetParameter(1);
      errval = fit1->GetParameter(2);
      cellsliceC[i]->SetTitle(Form("i%d m%0.2f s%0.2f max%0.2f",i,meanval,errval,maxcontent));

      if (meanval > meanXmax || meanval < meanXmin){ // Hard catch, if mean value out of known range use args
	cellsliceC[i]->SetLineColor(kOrange);
	continue;
      }
      //cell.push_back(cellval);
      //mean.push_back(meanval);
      //err.push_back(errval);
    }
    c->Update();
  }

  //slice histo function to include fits for all cells and write to c1 canvas. Hardcoded for HCal analysis.
  void sliceHCalIDHisto( TH2D *h2, Int_t Nslices, Double_t fwhm, Int_t min_ev, TCanvas *c1, TCanvas *c2, vector<Double_t> &cell, vector<Double_t> &mean, vector<Double_t> &err, Double_t meanXmin = -1e38, Double_t meanXmax = 1e38, Double_t noise_floor = 1e38 ){
  
    TH1D *hcalcellslice[Nslices];

    Int_t nBinsX = h2->GetNbinsX();

    // Assuming 288 bins in X
    if (nBinsX != hcal::maxHCalChan) {
      std::cerr << "Histogram does not have " << hcal::maxHCalChan << " bins in X" << std::endl;
        return;
    }

    c1->Divide(hcal::maxHCalCols,hcal::maxHCalRows/2); // 12 columns, 0-12 row (top half)
    c2->Divide(hcal::maxHCalCols,hcal::maxHCalRows/2); // 12 columns, 12-24 row (bottom half)

    for( Int_t i=0; i<Nslices; i++ ){
    
      Double_t cellval = (Double_t)i+0.5;
      Double_t meanval = 0.;
      Double_t errval = 0.;

      hcalcellslice[i] = h2->ProjectionY(Form("hcalcellslice_%d",i+1),i+1,i+1);

      if(i<144)
	c1->cd(i+1);
      else
	c2->cd(i-143);

      // Calculate the integral (total number of entries) in the range
      Double_t sliceN = hcalcellslice[i]->Integral(1, hcalcellslice[i]->GetNbinsX());

      // Get maximum value for scaling and cuts
      //Double_t histmax = hcalcellslice[i]->GetMaximum(); 
      Double_t histmax = get_max_bin_content_ignore_overflow(hcalcellslice[i]);

      hcalcellslice[i]->GetYaxis()->SetRangeUser(0,histmax*1.35);

      if( sliceN<min_ev ){ //continue if too sparse to get distribution
	hcalcellslice[i]->SetLineColor(kRed-5);
	hcalcellslice[i]->Draw();

	// Create a legend
	TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
	leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
	// Draw the legend
	leg->Draw();

	continue;
      }
      //get reasonable limits on fit
      Double_t min = hcalcellslice[i]->GetXaxis()->GetXmin();
      Double_t max = hcalcellslice[i]->GetXaxis()->GetXmax();
      Int_t totalbins = hcalcellslice[i]->GetNbinsX();
      Int_t binmax = hcalcellslice[i]->GetMaximumBin();

      if( histmax < noise_floor ){ //catch if all uniform noise
	hcalcellslice[i]->SetLineColor(kYellow-5);
	hcalcellslice[i]->Draw();

	// Create a legend
	TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
	leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
	// Draw the legend
	leg->Draw();

	continue;
      }

      Double_t binmaxX = min+binmax*(max-min)/totalbins;
      Double_t llim = binmaxX-fwhm;
      Double_t ulim = binmaxX+fwhm;
      if( llim < min ){ //if the lower limit is below minimum, recalculate
	binmax = get_max_bin_on_range( hcalcellslice[i], min, max );
	binmaxX = min+binmax*(max-min)/totalbins;
	llim = min;
	ulim = binmaxX + ( binmaxX - llim );
      }
    
      TF1 *fit1;
      hcalcellslice[i]->Fit("gaus","Q","",llim,ulim);
      fit1 = hcalcellslice[i]->GetFunction("gaus");
      meanval = fit1->GetParameter(1);
      errval = fit1->GetParameter(2);
      hcalcellslice[i]->SetTitle(Form("c:%d m:%0.2f s:%0.2f",i,meanval,errval));

      // Check if maximum bin is at the endpoints of the fit range
      Double_t tolerance = (ulim-llim)*0.025 ; // Define a small tolerance
      if (std::abs(meanval - llim) < tolerance || std::abs(meanval - ulim) < tolerance) {
	hcalcellslice[i]->SetLineColor(kMagenta);
      }

      if( meanval > meanXmax || meanval < meanXmin ){ //hard catch, if mean value out of known range use args
	hcalcellslice[i]->SetLineColor(kOrange-5);
	hcalcellslice[i]->Draw();

	// Create a legend
	TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
	leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
	// Draw the legend
	leg->Draw();

	continue;
      }

      hcalcellslice[i]->Draw();

      // Create a legend
      TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
      leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
      // Draw the legend
      leg->Draw();

      cell.push_back(cellval);
      mean.push_back(meanval);
      err.push_back(errval);
    }
  }

//slice histo function to include fits for all cells and write to c1 canvas. Hardcoded for HCal analysis.
  void sliceHCalIDHistoAlt( TH2D *h2, Int_t Nslices, Double_t fwhm, Int_t min_ev, TCanvas *c1, TCanvas *c2, vector<Double_t> &cell, vector<Double_t> &mean, vector<Double_t> &err, Double_t meanXmin = -1e38, Double_t meanXmax = 1e38, Double_t noise_floor = 1e38 ){
  
    TH1D *hcalcellslicealt[Nslices];

    Int_t nBinsX = h2->GetNbinsX();

    // Assuming 288 bins in X
    if (nBinsX != hcal::maxHCalChan) {
      std::cerr << "Histogram does not have " << hcal::maxHCalChan << " bins in X" << std::endl;
        return;
    }

    c1->Divide(hcal::maxHCalCols,hcal::maxHCalRows/2); // 12 columns, 0-12 row (top half)
    c2->Divide(hcal::maxHCalCols,hcal::maxHCalRows/2); // 12 columns, 12-24 row (bottom half)

    for( Int_t i=0; i<Nslices; i++ ){
    
      Double_t cellval = (Double_t)i+0.5;
      Double_t meanval = 0.;
      Double_t errval = 0.;

      hcalcellslicealt[i] = h2->ProjectionY(Form("hcalcellslicealt_%d",i+1),i+1,i+1);

      if(i<144)
	c1->cd(i+1);
      else
	c2->cd(i-143);

      // Calculate the integral (total number of entries) in the range
      Double_t sliceN = hcalcellslicealt[i]->Integral(1, hcalcellslicealt[i]->GetNbinsX());

      // Get maximum value for scaling and cuts
      Double_t histmax = hcalcellslicealt[i]->GetMaximum(); 
      hcalcellslicealt[i]->GetYaxis()->SetRangeUser(0,histmax*1.35);

      if( sliceN<min_ev ){ //continue if too sparse to get distribution
	hcalcellslicealt[i]->SetLineColor(kRed-5);
	hcalcellslicealt[i]->Draw();

	// Create a legend
	TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
	leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
	// Draw the legend
	leg->Draw();

	continue;
      }
      //get reasonable limits on fit
      Double_t min = hcalcellslicealt[i]->GetXaxis()->GetXmin();
      Double_t max = hcalcellslicealt[i]->GetXaxis()->GetXmax();
      Int_t totalbins = hcalcellslicealt[i]->GetNbinsX();
      Int_t binmax = hcalcellslicealt[i]->GetMaximumBin();

      if( histmax < noise_floor ){ //catch if all uniform noise
	hcalcellslicealt[i]->SetLineColor(kYellow-5);
	hcalcellslicealt[i]->Draw();

	// Create a legend
	TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
	leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
	// Draw the legend
	leg->Draw();

	continue;
      }

      Double_t binmaxX = min+binmax*(max-min)/totalbins;
      Double_t llim = binmaxX-fwhm;
      Double_t ulim = binmaxX+fwhm;
      if( llim < min ){ //if the lower limit is below minimum, recalculate
	binmax = get_max_bin_on_range( hcalcellslicealt[i], min, max );
	binmaxX = min+binmax*(max-min)/totalbins;
	llim = min;
	ulim = binmaxX + ( binmaxX - llim );
      }
    
      TF1 *fit1;
      hcalcellslicealt[i]->Fit("gaus","Q","",llim,ulim);
      fit1 = hcalcellslicealt[i]->GetFunction("gaus");
      meanval = fit1->GetParameter(1);
      errval = fit1->GetParameter(2);
      hcalcellslicealt[i]->SetTitle(Form("c:%d m:%0.2f s:%0.2f",i,meanval,errval));

      // Check if maximum bin is at the endpoints of the fit range
      Double_t tolerance = (ulim-llim)*0.025 ; // Define a small tolerance
      if (std::abs(meanval - llim) < tolerance || std::abs(meanval - ulim) < tolerance) {
	hcalcellslicealt[i]->SetLineColor(kMagenta);
      }

      if( meanval > meanXmax || meanval < meanXmin ){ //hard catch, if mean value out of known range use args
	hcalcellslicealt[i]->SetLineColor(kOrange-5);
	hcalcellslicealt[i]->Draw();

	// Create a legend
	TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
	leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
	// Draw the legend
	leg->Draw();

	continue;
      }

      hcalcellslicealt[i]->Draw();

      // Create a legend
      TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9); // x1, y1, x2, y2 coordinates (NDC)
      leg->AddEntry( (TObject*)0, Form("ev: %0.2f",sliceN), "");
      // Draw the legend
      leg->Draw();

      cell.push_back(cellval);
      mean.push_back(meanval);
      err.push_back(errval);
    }
  }

  //slice histo function configured to use course/fine gaussian fitting (fitGaussianAndGetFineParams())
  void sliceHistoFine( TH2D *h2, Int_t Nslices, Double_t fwhm, Int_t min_ev, vector<Double_t> &cell, vector<Double_t> &mean, vector<Double_t> &err ){
  
    TH1D *finecellslice[Nslices];

    for( Int_t i=0; i<Nslices; i++ ){
    
      Double_t cellval = (Double_t)i+0.5;
      Double_t meanval = 0.;
      Double_t errval = 0.;

      finecellslice[i] = h2->ProjectionY(Form("finecellslice_%d",i+1),i+1,i+1);
      Int_t sliceN = finecellslice[i]->GetEntries();
      if( sliceN<min_ev )
	continue;
    
      vector<Double_t> fitParams = fitGaussianAndGetFineParams(finecellslice[i],fwhm);

      meanval = fitParams[1];
      errval = fitParams[2];

      cell.push_back(cellval);
      mean.push_back(meanval);
      err.push_back(errval);
    }
  }

  //Fit functions
  //gaussian fit
  Double_t g_gfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
  }
 
  //gaussian fit with uniform offset
  Double_t g_gfit_bg(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t bg = par[3];
    return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.)) + bg;
  }
 
  //skewed gaussian fit
  Double_t g_sgfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];
    return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) );
  }

  //skewed gaussian fit with uniform offset
  Double_t g_sgfit_bg(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];
    Double_t bg = par[4];
    return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) ) + bg;
  }

  //expo fit
  Double_t g_expofit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t str = par[1];
    return exp(amp+str*x[0]);
  }

  //traditional timewalk fit
  Double_t g_tradtwfit(Double_t *x, Double_t *par){
    Double_t asymp = par[0];
    Double_t scale = par[1];
    Double_t expo = par[2];
    return asymp+scale/pow(x[0],expo);
  }

  //tuned exponential timewalk fit
  Double_t g_twfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t str = par[1];
    Double_t offset = par[2];
    return amp*exp(-str*x[0])+offset;
  } 

  //expo fit with offset
  Double_t g_scexpofit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t amp = par[1];
    Double_t str = par[2];
    return yint+exp(amp+str*x[0]);
  }

  //linear fit
  Double_t g_lfit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];;
    return yint+p1*x[0];
  }

  //3rd order poly fit
  Double_t g_p3fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3);
  }

  //4th order poly fit
  Double_t g_p4fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
  }

  //6th order poly fit
  Double_t g_p6fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6);
  }

} //::util
