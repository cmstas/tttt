void make_csv(){

	ofstream output_csv;
	output_csv.open("output_csv.csv");

	TFile* file = new TFile("DY_2l_M_50.root","read");
	TTree* tree = (TTree*)file->Get("SkimTree");

	output_csv << "Event" << ',';
	output_csv << "has_OS" << ',';
	output_csv << "nGenMatched_leptons" << ',';
	output_csv << "nFlips" << ',';
	output_csv << "nTightCharge_leading" << ',';
	output_csv << "nTightCharge_trailing" << ',';
	output_csv << "mother_leading_gen" << ',';	
	output_csv << "mother_trailing_gen" << endl;

	ULong64_t eventNumber = 0;
  std::vector<float>* leadingPdgId = nullptr;
  std::vector<float>* trailingPdgId = nullptr;
  std::vector<float>* leading_genmatchPdgId = nullptr;
  std::vector<float>* trailing_genmatchPdgId = nullptr;
  std::vector<float>* leading_tightCharge = nullptr;
  std::vector<float>* trailing_tightCharge = nullptr;
  std::vector<float>* leading_genmatch_momPdgId = nullptr;
  std::vector<float>* trailing_genmatch_momPdgId = nullptr;
	
	tree->SetBranchAddress("EventNumber", &eventNumber);
  tree->SetBranchAddress("dileptons_pdgId_1", &leadingPdgId);
  tree->SetBranchAddress("dileptons_pdgId_2", &trailingPdgId);
  tree->SetBranchAddress("dileptons_match_pdgId_1", &leading_genmatchPdgId);
  tree->SetBranchAddress("dileptons_match_pdgId_2", &trailing_genmatchPdgId);
  tree->SetBranchAddress("dileptons_leading_tightcharge", &leading_tightCharge);
  tree->SetBranchAddress("dileptons_trailing_tightcharge", &trailing_tightCharge);
  tree->SetBranchAddress("dileptons_gen_leading_mother", &leading_genmatch_momPdgId);
  tree->SetBranchAddress("dileptons_gen_trailing_mother", &trailing_genmatch_momPdgId);
		
	int entries = tree->GetEntries();

	for (int i=0; i<entries; i++){
		tree->GetEntry(i);
		output_csv << eventNumber << ',';
		
		int nGenMatched = 0, nFlips=0;		
				
		for (int j=0; j<leadingPdgId->size(); j++){

			bool hasOS = (((*leadingPdgId)[j]/(*trailingPdgId)[j]) == -1) ? 1 : 0;
			output_csv << hasOS << ',';
			
			if ((abs((*leadingPdgId)[j])/abs((*leading_genmatchPdgId)[j])) == 1) nGenMatched++;
			
			if (((*leadingPdgId)[j]/(*leading_genmatchPdgId)[j]) == -1) nFlips++; 	
				
			if ((abs((*trailingPdgId)[j])/abs((*trailing_genmatchPdgId)[j])) == 1) nGenMatched++;
			
			if (((*trailingPdgId)[j]/(*trailing_genmatchPdgId)[j]) == -1) nFlips++; 	
			
			output_csv << hasOS << ',';
			output_csv << nGenMatched << ',';
			output_csv << nFlips << ',';
			output_csv << (*leading_tightCharge)[j] << ',';
			output_csv << (*trailing_tightCharge)[j] << ',';
			output_csv << (*leading_genmatch_momPdgId)[j] << ',';
			output_csv << (*trailing_genmatch_momPdgId)[j] << endl;
		
			}
	
	}

output_csv.close(); file->Close();

}
