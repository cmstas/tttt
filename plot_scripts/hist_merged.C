void hist_merged(){
	
	TFile* file = new TFile("DY_2l_M_50.root","read");
	TFile* output = new TFile("output_merged.root","recreate");
	TTree* tree = (TTree*)file->Get("SkimTree");
	

	vector<float>* dileptons_leading_pt = nullptr;	
	vector<float>* dileptons_trailing_pt = nullptr;
	vector<float>* dileptons_leading_phi = nullptr;
	vector<float>* dileptons_trailing_phi = nullptr;
	vector<float>* dileptons_leading_eta = nullptr;
	vector<float>* dileptons_trailing_eta = nullptr;	
	vector<float>* dileptons_mass = nullptr;
	vector<float>* dileptons_pt_sum = nullptr;
	vector<float>* dileptons_pdgId_2 = nullptr;
	vector<float>* dileptons_match_pdgId_2 = nullptr;
	vector<float>* dileptons_pdgId_1 = nullptr;
	vector<float>* dileptons_match_pdgId_1 = nullptr;
	vector<int>* njets = nullptr;
	
	tree->SetBranchAddress("dileptons_lpt",&dileptons_leading_pt);
	tree->SetBranchAddress("dileptons_tpt",&dileptons_trailing_pt);
	tree->SetBranchAddress("dileptons_leading_phi",&dileptons_leading_phi);
	tree->SetBranchAddress("dileptons_trailing_phi",&dileptons_trailing_phi);
	tree->SetBranchAddress("dileptons_leading_eta",&dileptons_leading_eta);
	tree->SetBranchAddress("dileptons_trailing_eta",&dileptons_trailing_eta);
	tree->SetBranchAddress("dileptons_mass",&dileptons_mass);
	tree->SetBranchAddress("dileptons_pt_sum",&dileptons_pt_sum);
	tree->SetBranchAddress("dileptons_pdgId_2",&dileptons_pdgId_2);
	tree->SetBranchAddress("dileptons_match_pdgId_2",&dileptons_match_pdgId_2);
	tree->SetBranchAddress("dileptons_pdgId_1",&dileptons_pdgId_1);
	tree->SetBranchAddress("dileptons_match_pdgId_1",&dileptons_match_pdgId_1);
	tree->SetBranchAddress("dileptons_nJets",&njets);	

	TH1F* pt_leading_os = new TH1F("pt_leading_os","pt_leading_os",20,0,200);
	TH1F* eta_leading_os = new TH1F("eta_leading_os","eta_leading_os",20,0,3);
	TH1F* phi_leading_os = new TH1F("phi_leading_os","phi_leading_os",20,-5,5);

	TH1F* pt_trailing_os = new TH1F("pt_trailing_os","pt_trailing_os",20,0,200);
	TH1F* eta_trailing_os = new TH1F("eta_trailing_os","eta_trailing_os",20,0,3);
	TH1F* phi_trailing_os = new TH1F("phi_trailing_os","phi_trailing_os",20,-5,5);
	
	TH1F* boson_mass_os = new TH1F("zboson_mass_os","zboson_mass_os",40,0,200);
	TH1F* total_pt_os = new TH1F("total_pt_os","total_pt_os",50,0,200);
	TH2F* misidentification_hist_os = new TH2F("misidentification_plot_os","misidentification_plot_os",50,-15,15,70,-15,40);
	
	TH1F* njets_hist_os = new TH1F("nJets_os","nJets_os",20,0,10);	

	misidentification_hist_os->SetStats(0);
	misidentification_hist_os->SetFillColor(kGreen-9);

	TH1F* pt_leading_ss = new TH1F("pt_leading_ss","pt_leading_ss",20,0,200);
	TH1F* eta_leading_ss = new TH1F("eta_leading_ss","eta_leading_ss",20,0,3);
	TH1F* phi_leading_ss = new TH1F("phi_leading_ss","phi_leading_ss",20,-5,5);

	TH1F* pt_trailing_ss = new TH1F("pt_trailing_ss","pt_trailing_ss",20,0,200);
	TH1F* eta_trailing_ss = new TH1F("eta_trailing_ss","eta_trailing_ss",20,0,3);
	TH1F* phi_trailing_ss = new TH1F("phi_trailing_ss","phi_trailing_ss",20,-5,5);
	
	TH1F* boson_mass_ss = new TH1F("zboson_mass_ss","zboson_mass_ss",40,0,200);
	TH1F* total_pt_ss = new TH1F("total_pt_ss","total_pt_ss",50,0,200);
	TH2F* misidentification_hist_ss = new TH2F("misidentification_plot_ss","misidentification_plot_ss",50,-15,15,70,-15,40);
	
	TH1F* njets_hist_ss = new TH1F("nJets_ss","nJets_ss",20,0,10);	

	misidentification_hist_ss->SetStats(0);
	misidentification_hist_ss->SetFillColor(kGreen-9);
	int counter = 0;
	int entries = tree->GetEntries();

	for (int i=0; i<entries; i++){
		tree->GetEntry(i);
		for (int j=0; j<dileptons_mass->size(); j++){
				int product = (*dileptons_pdgId_1)[j] * (*dileptons_pdgId_2)[j];
				
				if (product < 0){
				boson_mass_os->Fill((*dileptons_mass)[j]);
				pt_leading_os->Fill((*dileptons_leading_pt)[j]);	
				pt_trailing_os->Fill((*dileptons_trailing_pt)[j]);
				eta_leading_os->Fill((*dileptons_leading_eta)[j]);
				eta_trailing_os->Fill((*dileptons_trailing_eta)[j]);
				phi_leading_os->Fill((*dileptons_leading_phi)[j]);
				phi_trailing_os->Fill((*dileptons_trailing_phi)[j]);
				total_pt_os->Fill((*dileptons_pt_sum)[j]);
				misidentification_hist_os->Fill((*dileptons_pdgId_1)[j],(*dileptons_match_pdgId_1)[j]);
				misidentification_hist_os->Fill((*dileptons_pdgId_2)[j],(*dileptons_match_pdgId_2)[j]);
				njets_hist_os->Fill((*njets)[j]);

				}
				else {
				counter++;
				boson_mass_ss->Fill((*dileptons_mass)[j]);
				pt_leading_ss->Fill((*dileptons_leading_pt)[j]);	
				pt_trailing_ss->Fill((*dileptons_trailing_pt)[j]);
				eta_leading_ss->Fill((*dileptons_leading_eta)[j]);
				eta_trailing_ss->Fill((*dileptons_trailing_eta)[j]);
				phi_leading_ss->Fill((*dileptons_leading_phi)[j]);
				phi_trailing_ss->Fill((*dileptons_trailing_phi)[j]);
				total_pt_ss->Fill((*dileptons_pt_sum)[j]);
				misidentification_hist_ss->Fill((*dileptons_pdgId_1)[j],(*dileptons_match_pdgId_1)[j]);
				misidentification_hist_ss->Fill((*dileptons_pdgId_2)[j],(*dileptons_match_pdgId_2)[j]);
				njets_hist_ss->Fill((*njets)[j]);
				}	

}			
}
	cout << counter << " yay" << endl;	
	output->Write();
	//pt_hist->Write();eta_hist->Write();phi_hist->Write();	
	//pt_leading->Write();eta_leading->Write();phi_leading->Write();	
	//pt_trailing->Write();eta_trailing->Write();phi_trailing->Write();
	//n_jets->Write();
	file->Close();output->Close();
	

}
