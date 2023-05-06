void make_csv(){

	ofstream output_csv;
	output_csv.open("/home/users/ymehra/CMSSW_10_6_26/src/tttt/test/output/data_run1_2017D/2017D/DoubleEG_2017D.csv");

	TFile* file = new TFile("/home/users/ymehra/CMSSW_10_6_26/src/tttt/test/output/data_run1_2017D/2017D/DoubleEG_2017D.root","read");
	TTree* tree = (TTree*)file->Get("SkimTree");

	output_csv << "Event" << ',';
	output_csv << "has_OS" << ',';
	output_csv << "nGenMatched_leptons" << ',';
	output_csv << "nFlips" << ',';
	output_csv << "nTightCharge_leading" << ',';
	output_csv << "nTightCharge_trailing" << ',';
	output_csv << "mother_leading_gen" << ',';	
	output_csv << "mother_trailing_gen" << ',';
	output_csv << "has_SS_gen" << ',';
	output_csv << "has_fake" << ',';
	output_csv << "leading_eta" << ',';
	output_csv << "trailing_eta" << endl;

	ULong64_t  eventNumber = 0;
  std::vector<float>* leadingPdgId = nullptr;
  std::vector<float>* trailingPdgId = nullptr;
  std::vector<float>* leading_genmatchPdgId = nullptr;
  std::vector<float>* trailing_genmatchPdgId = nullptr;
  std::vector<float>* leading_tightCharge = nullptr;
  std::vector<float>* trailing_tightCharge = nullptr;
  std::vector<float>* leading_genmatch_momPdgId = nullptr;
  std::vector<float>* trailing_genmatch_momPdgId = nullptr;
	std::vector<float>* leading_eta = nullptr;
	std::vector<float>* trailing_eta = nullptr;
	
	tree->SetBranchAddress("EventNumber", &eventNumber);
  tree->SetBranchAddress("dileptons_pdgId_1", &leadingPdgId);
  tree->SetBranchAddress("dileptons_pdgId_2", &trailingPdgId);
  tree->SetBranchAddress("dileptons_match_pdgId_1", &leading_genmatchPdgId);
  tree->SetBranchAddress("dileptons_match_pdgId_2", &trailing_genmatchPdgId);
  tree->SetBranchAddress("dileptons_leading_tightcharge", &leading_tightCharge);
  tree->SetBranchAddress("dileptons_trailing_tightcharge", &trailing_tightCharge);
  tree->SetBranchAddress("dileptons_gen_leading_mother", &leading_genmatch_momPdgId);
  tree->SetBranchAddress("dileptons_gen_trailing_mother", &trailing_genmatch_momPdgId);
	tree->SetBranchAddress("dileptons_leading_eta", &leading_eta);
	tree->SetBranchAddress("dileptons_trailing_eta", &trailing_eta);	
	
	int entries = tree->GetEntries();

	for (int i=0; i<entries; i++){
		tree->GetEntry(i);
		output_csv << eventNumber << ',';
		
		int nGenMatched = 0, nFlips=0;		
				
		for (int j=0; j<leadingPdgId->size(); j++){
			
			if (((*leading_genmatchPdgId)[j] == -13) || ((*trailing_genmatchPdgId)[j] == -13)) cout << eventNumber << endl;


			bool hasOS = ((*leadingPdgId)[j]/(*trailingPdgId)[j]) == -1;
			
			if ((*leading_genmatchPdgId)[j] != 33)  nGenMatched++; 			

			if (((*leadingPdgId)[j]/(float)(*leading_genmatchPdgId)[j]) == -1) nFlips++; 	
				
			if ((*trailing_genmatchPdgId)[j] != 33)  nGenMatched++;			

			if (((*trailingPdgId)[j]/(float)(*trailing_genmatchPdgId)[j]) == -1) nFlips++; 	
			
			bool has_SS_gen = (*leading_genmatchPdgId)[j]*(*trailing_genmatchPdgId)[j] == 121;				
			//bool has_fake = ((*leading_genmatchPdgId)[j] > 11) || ((*trailing_genmatchPdgId)[j] > 11);			
			bool has_fake = (abs((*leading_genmatchPdgId)[j]) != 11) || (abs((*trailing_genmatchPdgId)[j]) != 11); 

			output_csv << hasOS << ',';
			output_csv << nGenMatched << ',';
			output_csv << nFlips << ',';
			output_csv << (*leading_tightCharge)[j] << ',';
			output_csv << (*trailing_tightCharge)[j] << ',';
			output_csv << (*leading_genmatch_momPdgId)[j] << ',';
			output_csv << (*trailing_genmatch_momPdgId)[j] << ',';
			output_csv << has_SS_gen << ',';
			output_csv << has_fake << ',';
			output_csv << (*leading_eta)[j] << ',';
			output_csv << (*trailing_eta)[j] << endl;
		
			}
	
	}

output_csv.close(); file->Close();

}
