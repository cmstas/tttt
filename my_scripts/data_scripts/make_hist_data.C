void normalize_histogram(TH1F* hist, Double_t norm=1.0) {
    Double_t integral = hist->Integral();
    if (integral > 0) {
        hist->Scale(norm / integral);
    }
}

void hist_merged(){
	
	TFile* file = new TFile("/home/users/ymehra/CMSSW_10_6_26/src/tttt/test/output/data_run1_2017D/2017D/DoubleEG_2017D.root","read");
	TFile* file_sim = new TFile("~/CMSSW_10_6_26/src/tttt/test/output/all_data/2018/all_data/first_run/output_merged.root");
	TH2F* fR = (TH2F*) file_sim->Get("flip_rate");
	TFile* output = new TFile("/home/users/ymehra/CMSSW_10_6_26/src/tttt/test/output/data_run1_2017D/2017D/output.root","recreate");
	
	TTree* tree = (TTree*)file->Get("SkimTree");
	
	// declaring variables in which branch data will be stored
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
	vector<float>* gen_leading_pt = nullptr;
	vector<float>* gen_trailing_pt = nullptr;
	vector<float>* gen_leading_eta = nullptr;
	vector<float>* gen_trailing_eta = nullptr;
	vector<float>* gen_leading_phi = nullptr;
	vector<float>* gen_trailing_phi = nullptr;
	vector<float>* dileptons_leading_iso = nullptr;
	vector<float>* dileptons_trailing_iso = nullptr;
	vector<float>* dileptons_dxy_leading = nullptr;
	vector<float>* dileptons_dxy_trailing = nullptr;
	vector<float>* dileptons_dz_leading = nullptr;
	vector<float>* dileptons_dz_trailing = nullptr;
	vector<float>* dileptons_mom_leading_gen = nullptr;
	vector<float>* dileptons_mom_trailing_gen = nullptr;

	// assigning a branch to variables above
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
	tree->SetBranchAddress("dileptons_gen_pt_leading",&gen_leading_pt);
	tree->SetBranchAddress("dileptons_gen_pt_trailing",&gen_trailing_pt);
	tree->SetBranchAddress("dileptons_gen_eta_leading",&gen_leading_eta);
	tree->SetBranchAddress("dileptons_gen_eta_trailing",&gen_trailing_eta);
	tree->SetBranchAddress("dileptons_gen_phi_leading",&gen_leading_phi);
	tree->SetBranchAddress("dileptons_gen_phi_trailing",&gen_trailing_phi);
	tree->SetBranchAddress("dileptons_leading_iso",&dileptons_leading_iso);
	tree->SetBranchAddress("dileptons_trailing_iso",&dileptons_trailing_iso);
	tree->SetBranchAddress("dileptons_dxy_leading",&dileptons_dxy_leading);	
	tree->SetBranchAddress("dileptons_dxy_trailing",&dileptons_dxy_trailing);	
	tree->SetBranchAddress("dileptons_dz_leading",&dileptons_dz_leading);	
	tree->SetBranchAddress("dileptons_dz_trailing",&dileptons_dz_trailing);	
	tree->SetBranchAddress("dileptons_gen_leading_mother",&dileptons_mom_leading_gen);
	tree->SetBranchAddress("dileptons_gen_trailing_mother",&dileptons_mom_trailing_gen);

	// declaring histograms
	TH1F* pt_leading_os = new TH1F("pt_leading_os","pt_leading_os",20,15,300);
	TH1F* eta_leading_os = new TH1F("eta_leading_os","eta_leading_os",20,0,2.501);
	TH1F* phi_leading_os = new TH1F("phi_leading_os","phi_leading_os",20,-5,5);

	TH1F* pt_trailing_os = new TH1F("pt_trailing_os","pt_trailing_os",20,0,200);
	TH1F* eta_trailing_os = new TH1F("eta_trailing_os","eta_trailing_os",20,0,3);
	TH1F* phi_trailing_os = new TH1F("phi_trailing_os","phi_trailing_os",20,-5,5);
	
	TH1F* boson_mass_os = new TH1F("zboson_mass_os","zboson_mass_os",50,50,110);
	TH1F* total_pt_os = new TH1F("total_pt_os","total_pt_os",50,0,200);
	TH2F* misidentification_hist_os = new TH2F("misidentification_plot_os","misidentification_plot_os",50,-15,15,70,-15,40);
	
	TH1F* njets_hist_os = new TH1F("nJets_os","nJets_os",20,0,10);	

	misidentification_hist_os->SetStats(0);

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

	TH1F* pt_correct_matched = new TH1F("SS_pt_correct_matched","SS_pt_correct_matched",6,15,300);
	TH1F* pt_opposite_matched = new TH1F("SS_pt_opposite_matched","SS_pt_opposite_matched",6,15,300);
	TH1F* pt_fake_matched = new TH1F("SS_pt_fake_matched","SS_pt_fake_matched",6,15,300);
	
	TH1F* eta_correct_matched = new TH1F("SS_eta_correct_matched","SS_eta_correct_matched",3,0,2.5);
	TH1F* eta_opposite_matched = new TH1F("SS_eta_opposite_matched","SS_eta_opposite_matched",3,0,2.5);
	TH1F* eta_fake_matched = new TH1F("SS_eta_fake_matched","SS_eta_fake_matched",3,0,2.5);

	Double_t pt_bins[] = {15, 40, 60, 80, 100, 200, 300}; Int_t nbins_pt = sizeof(pt_bins)/sizeof(Double_t) - 1;
	Double_t eta_bins[] = {0.0, 0.80, 1.479, 2.501}; Int_t nbins_eta = sizeof(eta_bins)/sizeof(Double_t) - 1;

	TH2F* fraction_flips = new TH2F("flip_rate_1","flip_rate_1",nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* flip_rate_den = new TH2F("flip_rate_den","flip_rate_den",nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* os_reco_pt_eta = new TH2F("os_reco", "os_reco", nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* expected_SS = new TH2F("ss_expected", "ss_expected", nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* true_SS = new TH2F("ss_true", "ss_true", nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* difference_SS = new TH2F("difference_SS", "difference", nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* factor_SS = new TH2F("factor_SS", "factor_SS", nbins_pt, pt_bins, nbins_eta, eta_bins); 
	TH2F* adjusted_fR = new TH2F("adjused_fR", "adjusted_fR", nbins_pt, pt_bins, nbins_eta, eta_bins);
	TH2F* new_expected = new TH2F("new_expected","new_expected", nbins_pt, pt_bins, nbins_eta, eta_bins);
	
	
	TH2F* pt_eta_SS_matched = new TH2F("pt_eta_SS_matched","pt_eta_SS_matched",6,15,300,3,0,2.5);
	TH2F* pt_eta_OS_matched = new TH2F("pt_eta_OS_matched","pt_eta_OS_matched",6,15,300,3,0,2.5);
	TH2F* pt_eta_fake_matched = new TH2F("pt_eta_fake_matched","pt_eta_fake_matched",6,15,300,3,0,2.5);

	TH1F* pt_gen_SS_matched = new TH1F("pt_gen_SS_matched","pt_gen_SS_matched",6,15,300); 
	TH1F* pt_gen_OS_matched = new TH1F("pt_gen_OS_matched","pt_gen_OS_matched",6,15,300); 
	TH1F* pt_gen_fake_matched = new TH1F("pt_gen_fake_matched","pt_gen_fake_matched",6,15,300); 
	
	TH1F* eta_gen_SS_matched = new TH1F("eta_gen_SS_matched","eta_gen_SS_matched",3,0,2.5); 
	TH1F* eta_gen_OS_matched = new TH1F("eta_gen_OS_matched","eta_gen_OS_matched",3,0,2.5); 
	TH1F* eta_gen_fake_matched = new TH1F("eta_gen_fake_matched","eta_gen_fake_matched",3,0,2.5); 

	TH1F* dR_SS_matched = new TH1F("dR_SS_matched","dR_SS_matched",100,0,0.2);
	TH1F* dR_OS_matched = new TH1F("dR_OS_matched","dR_OS_matched",100,0,0.2);

	TH1F* iso_SS = new TH1F("iso_SS","iso_SS",20,0,0.1);
	TH1F* iso_OS = new TH1F("iso_OS","iso_OS",20,0,0.1);
	TH1F* iso_fake = new TH1F("iso_fake","iso_fake",20,0,0.1);

	TH1F* dxy_SS = new TH1F("dxy_SS","dxy_SS",20,0,0.1);
	TH1F* dxy_OS = new TH1F("dxy_OS","dxy_OS",20,0,0.1);
	TH1F* dxy_fake = new TH1F("dxy_fake","dxy_fake",20,0,0.1);
	
	TH1F* dz_SS = new TH1F("dz_SS","dz_SS",20,0,0.1);
	TH1F* dz_OS = new TH1F("dz_OS","dz_OS",20,0,0.1);
	TH1F* dz_fake = new TH1F("dz_fake","dz_fake",20,0,0.1);

	TH1F* gen_mom_SS = new TH1F("gen_mother_SS","gen_mother_SS",60,30,-30);
	TH1F* gen_mom_OS = new TH1F("gen_mother_OS","gen_mother_OS",60,30,-30);
	TH1F* gen_mom_fake = new TH1F("gen_mother_fake","gen_mother_fake",60,-30,30);
	
	// getting entriess
	int entries = tree->GetEntries();

	for (int i=0; i<entries; i++){
		tree->GetEntry(i);
		for (int j=0; j<dileptons_mass->size(); j++){

				int part1_id = (*dileptons_pdgId_1)[j], part2_id = (*dileptons_pdgId_2)[j];
				int match_part1_id = (*dileptons_match_pdgId_1)[j], match_part2_id = (*dileptons_match_pdgId_2)[j];
			
				bool gen_SS_event = ((match_part1_id * match_part2_id) == 121);
				bool correctly_matched_part1 = (part1_id * match_part1_id == 121), correctly_matched_part2 = (part2_id * match_part2_id == 121);
				bool opposite_matched_part1 = (part1_id * match_part1_id == -121), opposite_matched_part2 = (part2_id * match_part2_id == -121);
				bool not_matched_1 = (!correctly_matched_part1 && !opposite_matched_part1), not_matched_2 = (!correctly_matched_part2 && !opposite_matched_part2);
				bool opposite_event = (opposite_matched_part1 && opposite_matched_part2);
				
				bool is_valid = ((!not_matched_1) && (!not_matched_2) && (!gen_SS_event));
 
				float dR_leading = sqrt(pow(((*dileptons_leading_phi)[j] - (*gen_leading_phi)[j]),2)	+ pow(((*dileptons_leading_eta)[j] - (*gen_leading_eta)[j]),2)); 			
				
				float dR_trailing = sqrt(pow(((*dileptons_trailing_phi)[j] - (*gen_trailing_phi)[j]),2)	+ pow(((*dileptons_trailing_eta)[j] - (*gen_trailing_eta)[j]),2)); 			
				int product = part1_id * part2_id;

				double pt_leading = (*dileptons_leading_pt)[j], pt_trailing = (*dileptons_trailing_pt)[j];
				
				if (pt_leading > 300) pt_leading = 299.9;
				if (pt_trailing > 300) pt_trailing = 299.9;

				if (is_valid){

					flip_rate_den->Fill(pt_leading, abs((*dileptons_leading_eta)[j]));
					if (opposite_matched_part1) fraction_flips->Fill(pt_leading ,abs((*dileptons_leading_eta)[j]));

					flip_rate_den->Fill(pt_trailing ,abs((*dileptons_trailing_eta)[j]));
					 if (opposite_matched_part2) fraction_flips->Fill(pt_trailing ,abs((*dileptons_trailing_eta)[j]));

				}
					
				if (product < 0){

								
								if (true){
									os_reco_pt_eta->Fill(pt_leading ,abs((*dileptons_leading_eta)[j]));
									os_reco_pt_eta->Fill(pt_trailing ,abs((*dileptons_trailing_eta)[j]));	
								}
								
							
								/*boson_mass_os->Fill((*dileptons_mass)[j]);
								pt_leading_os->Fill((*dileptons_leading_pt)[j]);	
								pt_trailing_os->Fill((*dileptons_trailing_pt)[j]);
								eta_leading_os->Fill((*dileptons_leading_eta)[j]);
								eta_trailing_os->Fill((*dileptons_trailing_eta)[j]);
								phi_leading_os->Fill((*dileptons_leading_phi)[j]);
								phi_trailing_os->Fill((*dileptons_trailing_phi)[j]);
								total_pt_os->Fill((*dileptons_pt_sum)[j]);
								misidentification_hist_os->Fill((*dileptons_pdgId_1)[j],(*dileptons_match_pdgId_1)[j]);
								misidentification_hist_os->Fill((*dileptons_pdgId_2)[j],(*dileptons_match_pdgId_2)[j]);
								//njets_hist_os->Fill((*njets)[j]);*/
																									
				}

				else {
							if (true){
										
								true_SS->Fill(pt_leading ,abs((*dileptons_leading_eta)[j]));
								true_SS->Fill(pt_trailing ,abs((*dileptons_trailing_eta)[j]));	
							}
								/*if (correctly_matched_part1){
								
									gen_mom_SS->Fill((*dileptons_mom_leading_gen)[j]);

									dz_SS->Fill((*dileptons_dz_leading)[j]);
									
									dxy_SS->Fill((*dileptons_dxy_leading)[j]);
									iso_SS->Fill((*dileptons_leading_iso)[j]);
									dR_SS_matched->Fill(dR_leading);	
									pt_gen_SS_matched->Fill((*gen_leading_pt)[j]); eta_gen_SS_matched->Fill(abs((*gen_leading_eta)[j]));			

									pt_eta_SS_matched->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));
									pt_correct_matched->Fill((*dileptons_leading_pt)[j]);	
									eta_correct_matched->Fill(abs((*dileptons_leading_eta)[j]));		
								}
								
								else if (opposite_matched_part1){
									gen_mom_OS->Fill((*dileptons_mom_leading_gen)[j]);
									
									dz_OS->Fill((*dileptons_dz_leading)[j]);	
									iso_OS->Fill((*dileptons_leading_iso)[j]);
									dxy_OS->Fill((*dileptons_dxy_leading)[j]);
									dR_OS_matched->Fill(dR_leading);
									pt_gen_OS_matched->Fill((*gen_leading_pt)[j]); eta_gen_OS_matched->Fill(abs((*gen_leading_eta)[j]));			
					
									pt_eta_OS_matched->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));
									pt_opposite_matched->Fill((*dileptons_leading_pt)[j]);	
									eta_opposite_matched->Fill(abs((*dileptons_leading_eta)[j]));		
							
								}
								
								else{

									gen_mom_fake->Fill((*dileptons_mom_leading_gen)[j]);

								//	pt_gen_fake_matched->Fill((*gen_leading_pt)[j]); eta_gen_fake_matched->Fill(abs((*gen_leading_eta)[j]));
									dz_fake->Fill((*dileptons_dz_leading)[j]);	
									dxy_fake->Fill((*dileptons_dxy_leading)[j]);
									iso_fake->Fill((*dileptons_leading_iso)[j]);
									pt_eta_fake_matched->Fill((*dileptons_leading_pt)[j],abs((*dileptons_leading_eta)[j]));
									pt_fake_matched->Fill((*dileptons_leading_pt)[j]);	
									eta_fake_matched->Fill(abs((*dileptons_leading_eta)[j]));		
									
								}	


								if (correctly_matched_part2){
									
									gen_mom_SS->Fill((*dileptons_mom_trailing_gen)[j]);
									
									dz_SS->Fill((*dileptons_dz_trailing)[j]);
									dxy_SS->Fill((*dileptons_dxy_trailing)[j]);
									dR_SS_matched->Fill(dR_trailing);  
									iso_SS->Fill((*dileptons_trailing_iso)[j]);

									pt_gen_SS_matched->Fill((*gen_trailing_pt)[j]); eta_gen_SS_matched->Fill(abs((*gen_trailing_eta)[j]));
									pt_eta_SS_matched->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));

									pt_correct_matched->Fill((*dileptons_trailing_pt)[j]);	
									eta_correct_matched->Fill(abs((*dileptons_trailing_eta)[j]));		
								}
								
								else if (opposite_matched_part2){

									gen_mom_OS->Fill((*dileptons_mom_trailing_gen)[j]);
								
									dz_OS->Fill((*dileptons_dz_trailing)[j]);
									dxy_OS->Fill((*dileptons_dxy_trailing)[j]);
									iso_OS->Fill((*dileptons_trailing_iso)[j]);
									dR_OS_matched->Fill(dR_trailing);	
									pt_gen_OS_matched->Fill((*gen_trailing_pt)[j]); eta_gen_OS_matched->Fill(abs((*gen_trailing_eta)[j]));

									pt_eta_OS_matched->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));
									pt_opposite_matched->Fill((*dileptons_trailing_pt)[j]);	
									eta_opposite_matched->Fill(abs((*dileptons_trailing_eta)[j]));		
							
								}
								
								else{
									gen_mom_fake->Fill((*dileptons_mom_trailing_gen)[j]);

									dz_fake->Fill((*dileptons_dz_trailing)[j]);
									dxy_fake->Fill((*dileptons_dxy_trailing)[j]);
								//	pt_gen_fake_matched->Fill((*gen_trailing_pt)[j]); eta_gen_fake_matched->Fill(abs((*gen_trailing_eta)[j]));
									pt_eta_fake_matched->Fill((*dileptons_trailing_pt)[j],abs((*dileptons_trailing_eta)[j]));
									iso_fake->Fill((*dileptons_trailing_iso)[j]);
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
				njets_hist_ss->Fill((*njets)[j]);*/
				}	

	}			
}

	// pt_SS_histograms
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
	
	// end of pt_SS_histograms
	
	// eta_SS_histograms
	canvas = new TCanvas("eta_SS", "eta_SS", 800, 600);

	eta_correct_matched->SetLineColor(kRed);eta_opposite_matched->SetLineColor(kBlue);eta_fake_matched->SetLineColor(kGreen);
	
	legend = new TLegend(0.7, 0.7, 0.9, 0.9);
		
	legend->AddEntry(eta_correct_matched,"correctly matched","l");legend->AddEntry(eta_opposite_matched,"opposite matched","l");
	legend->AddEntry(eta_fake_matched,"fake matched","l");
	eta_opposite_matched->Draw(); eta_correct_matched->Draw("SAME"); eta_fake_matched->Draw("SAME"); legend->Draw("SAME");
	
	canvas->Write();
	canvas->Clear();
	// end of eta_SS_histograms
	
	// pt_gen_match_histograms
	
	canvas = new TCanvas("pt_gen", "pt_gen", 800, 600);
		
	pt_gen_SS_matched->SetLineColor(kRed); pt_gen_OS_matched->SetLineColor(kBlue);// pt_gen_fake_matched->SetLineColor(kGreen);

	legend = new TLegend(0.7,0.7,0.9,0.9);
		
	legend->AddEntry(pt_gen_SS_matched,"correctly matched","l");legend->AddEntry(pt_gen_OS_matched,"opposite matched","l");
	//legend->AddEntry(pt_gen_fake_matched,"fake matched","l");

	pt_gen_SS_matched->Draw();pt_gen_OS_matched->Draw("SAME");/*pt_gen_fake_matched->Draw("SAME");*/legend->Draw("SAME");

	canvas->Write();
	canvas->Clear();
	// end of pt_gen_match_histograms
	
	// eta_gen_match_histograms
	
	canvas = new TCanvas("eta_gen", "eta_gen", 800, 600);
		
	eta_gen_SS_matched->SetLineColor(kRed); eta_gen_OS_matched->SetLineColor(kBlue); //eta_gen_fake_matched->SetLineColor(kGreen);

	legend = new TLegend(0.7,0.7,0.9,0.9);
		
	legend->AddEntry(eta_gen_SS_matched,"correctly matched","l");legend->AddEntry(eta_gen_OS_matched,"opposite matched","l");
	//legend->AddEntry(eta_gen_fake_matched,"fake matched","l");

	eta_gen_OS_matched->Draw();eta_gen_SS_matched->Draw("SAME");/*eta_gen_fake_matched->Draw("SAME");*/legend->Draw("SAME");

	canvas->Write();
	canvas->Clear();
  // end of eta_gen_match_histograms

	// dR to genmatched plot
	
	dR_SS_matched->SetLineColor(kRed); dR_OS_matched->SetLineColor(kBlue);
	
	canvas = new TCanvas("dR_to_gen", "dR_to_gen", 800, 600);
	legend = new TLegend(0.7,0.7,0.9,0.9);

	legend->AddEntry(dR_SS_matched,"correctly matched","l");legend->AddEntry(dR_OS_matched,"opposite matched","l");
	dR_SS_matched->Draw(); dR_OS_matched->Draw("SAME"); legend->Draw("SAME");
	canvas->Write();
	canvas->Clear();	

	// end of dR plot
	

	// rel iso plot
	
	iso_SS->SetLineColor(kRed); iso_OS->SetLineColor(kBlue); iso_fake->SetLineColor(kGreen);

	legend = new TLegend(0.7,0.7,0.9,0.9);

	canvas = new TCanvas("relative_iso","relative_iso",800,400);
	legend->AddEntry(iso_SS,"iso_SS","l");legend->AddEntry(iso_OS,"iso_OS","l");legend->AddEntry(iso_fake,"iso_fake","l");

	iso_SS->Draw(); iso_OS->Draw("SAME"); iso_fake->Draw("SAME"); legend->Draw("SAME");

	canvas->Write();
	// end of iso plot

	// dxy plot
	
	dxy_SS->SetLineColor(kRed); dxy_OS->SetLineColor(kBlue); dxy_fake->SetLineColor(kGreen);

	legend = new TLegend(0.7,0.7,0.9,0.9);

	canvas = new TCanvas("dxy_parameter","dxy_parameter",800,400);
	legend->AddEntry(dxy_SS,"dxy_SS","l");legend->AddEntry(dxy_OS,"dxy_OS","l");legend->AddEntry(dxy_fake,"dxy_fake","l");

	dxy_SS->Draw(); dxy_OS->Draw("SAME"); dxy_fake->Draw("SAME"); legend->Draw("SAME");

	canvas->Write();

	// end of dxy plot

	// dz plot
	
	dz_SS->SetLineColor(kRed); dz_OS->SetLineColor(kBlue); dz_fake->SetLineColor(kGreen);

	legend = new TLegend(0.7,0.7,0.9,0.9);

	canvas = new TCanvas("dz_parameter","dz_parameter",800,400);
	legend->AddEntry(dz_SS,"dz_SS","l");legend->AddEntry(dz_OS,"dz_OS","l");legend->AddEntry(dz_fake,"dz_fake","l");

	dz_SS->Draw(); dz_OS->Draw("SAME"); dz_fake->Draw("SAME"); dz_fake->Draw("SAME"); legend->Draw("SAME");

	canvas->Write();
	
	// end of dz plot
	
	// mother gen plot
	
	gen_mom_SS->SetLineColor(kRed); gen_mom_OS->SetLineColor(kBlue); gen_mom_fake->SetLineColor(kGreen);	
	
	legend = new TLegend(0.7,0.7,0.9,0.9);	

	canvas = new TCanvas("gen_mom_pdgId","gen_mom_pdgId",800,400);

	legend->AddEntry(gen_mom_SS,"gen_mom_SS","l");legend->AddEntry(gen_mom_OS,"gen_mom_OS","l");legend->AddEntry(gen_mom_fake,"gen_mom_fake","l");

	gen_mom_SS->Draw(); gen_mom_OS->Draw("SAME"); gen_mom_fake->Draw("SAME"); legend->Draw("SAME");

	canvas->Write();

	// end of gen mother plot
	

	TH2F* flip_rate_num = (TH2F*)fraction_flips->Clone("flip_rate_num");
	fraction_flips->Divide(flip_rate_den);
	
	// start of expected_SS plot

	for (int i=1; i<=os_reco_pt_eta->GetNbinsX(); i++){
		for (int j=1; j<=os_reco_pt_eta->GetNbinsY(); j++){
		
			double current = os_reco_pt_eta->GetBinContent(i, j);
			double number = fR->GetBinContent(i, j);
			expected_SS->SetBinContent(i, j, current*(number/(1-number)));		

			}

	}

	// end of expected_SS plot

	// start of difference plot

	for (int i=1; i<=difference_SS->GetNbinsX(); i++){
		for (int j=1; j<=difference_SS->GetNbinsY(); j++){

			double difference = abs(true_SS->GetBinContent(i, j) - expected_SS->GetBinContent(i, j));
			difference_SS->SetBinContent(i, j, difference);		
		}

	}

	// end of difference plot
	
	// new plots
	
	TH2F* scaledHist = (TH2F*) true_SS->Clone("clone_true_SS");
	scaledHist->Scale(0.5);	
	factor_SS->Divide(scaledHist, expected_SS);
	adjusted_fR->Multiply(fR,factor_SS);
	new_expected->Multiply(adjusted_fR, os_reco_pt_eta );				
	
	// end of new plots	

	output->Write();

	output->Close();
	file_sim->Close();
	file->Close();
}
