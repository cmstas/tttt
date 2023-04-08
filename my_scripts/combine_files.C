void combine_files() {
  // Create a new TChain object
  TChain *chain = new TChain("SkimTree");

  // Loop over all the files and add them to the chain
  for (int i=0; i<283; i++) {	// file 0 through file 283
    TString filename = Form("/ceph/cms/store/group/tttt/Worker/ymehra/output/230209_DY_2l_M_50/2018/DY_2l_M_50_%d_of_283.root", i); 
    chain->Add(filename);
  }

  // Create a new output file and clone the TTree from the first file
  TFile *outfile = new TFile("combined_DY_2l_M_50.root", "RECREATE");
  TTree *outtree = chain->CloneTree(0);

  // Loop over all entries in the chain and fill the output tree
  for (int i=0; i<chain->GetEntries(); i++) {
    chain->GetEntry(i);
    outtree->Fill();
  }

  // Write the output file and clean up
  outfile->Write();
  outfile->Close();
  delete chain;
}
