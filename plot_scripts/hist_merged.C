void normalize_histogram(TH1F* hist, Double_t norm=1.0) {
    Double_t integral = hist->Integral();
    if (integral > 0) {
        hist->Scale(norm / integral);
    }
}

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
	
	TH1F* boson_mass_os = new TH1F("zboson_mass_os","zboson_mass_os",50,50,110);
	TH1F* total_pt_os = new TH1F("total_pt_os","total_pt_os",50,0,200);
	TH2F* misidentification_hist_os = new TH2F("misidentification_plot_os","misidentification_plot_os",50,-15,15,70,-15,40);
	
	TH1F* njets_hist_os = new TH1F("nJets_os","nJets_os",20,0,10);	

	misidentification_hist_os->SetStats(0);
	//misidentification_hist_os->SetFillColor(kGreen-9);

	TH1F* pt_leading_ss = new TH1F("pt_leading_ss","pt_leading_ss",20,0,200);
	TH1F* eta_leading_ss = new TH1F("eta_leading_ss","eta_leading_ss",20,0,3);
	TH1F* phi_leading_ss = new TH1F("phi_leading_ss","phi_leading_ss",20,-5,5);

	TH1F* pt_trailing_ss = new TH1F("pt_trailing_ss","pt_trailing_ss",20,0,200);
	TH1F* eta_trailing_ss = new TH1F("eta_trailing_ss","eta_trailing_ss",20,0,3);
	TH1F* phi_trailing_ss = new TH1F("phi_trailing_ss","phi_trailing_ss",20,-5,5);
	
	TH1F* boson_mass_ss = new TH1F("zboson_mass_ss","zboson_mass_ss",50,50,110);
	TH1F* total_pt_ss = new TH1F("total_pt_ss","total_pt_ss",50,0,200);
	TH2F* misidentification_hist_ss = new TH2F("misidentification_plot_ss","misidentification_plot_ss",50,-15,15,70,-15,40);
	
	TH1F* njets_hist_ss = new TH1F("nJets_ss","nJets_ss",20,0,10);	

	misidentification_hist_ss->SetStats(0);
	//misidentification_hist_ss->SetFillColor(kGreen-9);

	TH1F* pt_correct_matched = new TH1F("SS_pt_correct_matched","SS_pt_correct_matched",6,15,300);
	TH1F* pt_opposite_matched = new TH1F("SS_pt_opposite_matched","SS_pt_opposite_matched",6,15,300);
	TH1F* pt_fake_matched = new TH1F("SS_pt_fake_matched","SS_pt_fake_matched",6,15,300);
	
	TH1F* eta_correct_matched = new TH1F("SS_eta_correct_matched","SS_eta_correct_matched",3,0,2.5);
	TH1F* eta_opposite_matched = new TH1F("SS_eta_opposite_matched","SS_eta_opposite_matched",3,0,2.5);
	TH1F* eta_fake_matched = new TH1F("SS_eta_fake_matched","SS_eta_fake_matched",3,0,2.5);

	TH2F* fraction_flips = new TH2F("fraction_flips","fraction_flips",6,15,300,3,0,2.5);
	TH2F* flip_rate_den = new TH2F("flip_rate_den","flip_rate_den",6,15,300,3,0,2.5);
//	TH2F* fraction_flips = new TH2F("fraction_flips","fraction_flips",50,0,200,20,0,3);

	TH2F* pt_eta_SS_matched = new TH2F("pt_eta_SS_matched","pt_eta_SS_matched",6,15,300,3,0,2.5);
 
	TH2F* pt_eta_OS_matched = new TH2F("pt_eta_OS_matched","pt_eta_OS_matched",6,15,300,3,0,2.5);

	TH2F* pt_eta_fake_matched = new TH2F("pt_eta_fake_matched","pt_eta_fake_matched",6,15,300,3,0,2.5);

	int entries = tree->GetEntries();
	int counter = 0;
	for (int i=0; i<entries; i++){
		tree->GetEntry(i);
		for (int j=0; j<dileptons_mass->size(); j++){

				int part1_id = (*dileptons_pdgId_1)[j], part2_id = (*dileptons_pdgId_2)[j];
				int match_part1_id = (*dileptons_match_pdgId_1)[j], match_part2_id = (*dileptons_match_pdgId_2)[j];

				bool correctly_matched_part1 = (part1_id * match_part1_id == 121), correctly_matched_part2 = (part2_id * match_part2_id == 121);
				bool opposite_matched_part1 = (part1_id * match_part1_id == -121), opposite_matched_part2 = (part2_id * match_part2_id == -121);
				bool not_matched_1 = (!correctly_matched_part1 && !opposite_matched_part1), not_matched_2 = (!correctly_matched_part2 && !opposite_matched_part2); 

				int product = part1_id * part2_id;

				if (opposite_matched_part1 || correctly_matched_part1){
										
					flip_rate_den->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));

					if (opposite_matched_part1)
						fraction_flips->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));

				}
				
							
				if (opposite_matched_part2 || correctly_matched_part2){
										
					flip_rate_den->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));

					if (opposite_matched_part2)
						fraction_flips->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));

				}

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

								if (correctly_matched_part1){
									pt_eta_SS_matched->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));
									pt_correct_matched->Fill((*dileptons_leading_pt)[j]);	
									eta_correct_matched->Fill(abs((*dileptons_leading_eta)[j]));		
								}
								
								else if (opposite_matched_part1){
									pt_eta_OS_matched->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));
									pt_opposite_matched->Fill((*dileptons_leading_pt)[j]);	
									eta_opposite_matched->Fill(abs((*dileptons_leading_eta)[j]));		
							
								}
								
								else{
									pt_eta_fake_matched->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));
									pt_fake_matched->Fill((*dileptons_leading_pt)[j]);	
									eta_fake_matched->Fill(abs((*dileptons_leading_eta)[j]));		
									
								}	


								if (correctly_matched_part2){

									pt_eta_SS_matched->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));

									pt_correct_matched->Fill((*dileptons_trailing_pt)[j]);	
									eta_correct_matched->Fill(abs((*dileptons_trailing_eta)[j]));		
								}
								
								else if (opposite_matched_part2){
									
									pt_eta_OS_matched->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));
									pt_opposite_matched->Fill((*dileptons_trailing_pt)[j]);	
									eta_opposite_matched->Fill(abs((*dileptons_trailing_eta)[j]));		
							
								}
								
								else{

									pt_eta_fake_matched->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));
									pt_fake_matched->Fill((*dileptons_trailing_pt)[j]);	
									eta_fake_matched->Fill(abs((*dileptons_trailing_eta)[j]));		
								}	


				boson_mass_ss->Fill((*dileptons_mass)[j]);
				pt_leading_ss->Fill((*dileptons_leading_pt)[j]);	
				pt_trailing_ss->Fill((*dileptons_trailing_pt)[j]);
				eta_leading_ss->Fill(abs((*dileptons_leading_eta)[j]));
				eta_trailing_ss->Fill(abs((*dileptons_trailing_eta)[j]));
				phi_leading_ss->Fill((*dileptons_leading_phi)[j]);
				phi_trailing_ss->Fill((*dileptons_trailing_phi)[j]);
				total_pt_ss->Fill((*dileptons_pt_sum)[j]);
				misidentification_hist_ss->Fill((*dileptons_pdgId_1)[j],(*dileptons_match_pdgId_1)[j]);
				misidentification_hist_ss->Fill((*dileptons_pdgId_2)[j],(*dileptons_match_pdgId_2)[j]);
				njets_hist_ss->Fill((*njets)[j]);
				}	

}			
}


	pt_correct_matched->SetLineColor(kRed);	
	
	pt_opposite_matched->SetLineColor(kBlue);	

	pt_fake_matched->SetLineColor(kGreen);	
	
	//normalize_histogram(pt_correct_matched,1.0);normalize_histogram(pt_opposite_matched,1.0);normalize_histogram(pt_fake_matched,1.0);
	
	TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
	
	legend->AddEntry(pt_correct_matched,"correctly matched","l");legend->AddEntry(pt_opposite_matched,"opposite matched","l");
	legend->AddEntry(pt_fake_matched,"fake matched","l");
	TCanvas* canvas = new TCanvas("pt_SS", "pt_SS", 800, 600);
	
	pt_correct_matched->Draw(); pt_opposite_matched->Draw("SAME"); pt_fake_matched->Draw("SAME"); legend->Draw("SAME");

	canvas->Write();

	canvas = new TCanvas("eta_SS", "eta_SS", 800, 600);

	eta_correct_matched->SetLineColor(kRed);eta_opposite_matched->SetLineColor(kBlue);eta_fake_matched->SetLineColor(kGreen);
	
	legend = new TLegend(0.7, 0.7, 0.9, 0.9);
		
	legend->AddEntry(eta_correct_matched,"correctly matched","l");legend->AddEntry(eta_opposite_matched,"opposite matched","l");
	legend->AddEntry(eta_fake_matched,"fake matched","l");
	eta_opposite_matched->Draw(); eta_correct_matched->Draw("SAME"); eta_fake_matched->Draw("SAME"); legend->Draw("SAME");
	
	canvas->Write();
	
	TH2F* flip_rate_num = (TH2F*)fraction_flips->Clone("flip_rate_num");
	fraction_flips->Divide(flip_rate_den);
	output->Write();
	file->Close();output->Close();
	
}
