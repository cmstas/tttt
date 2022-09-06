#include <string>
#include <iostream>
#include <iomanip>
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "IvyFramework/IvyDataTools/interface/HelperFunctions.h"
#include "IvyFramework/IvyDataTools/interface/IvyStreamHelpers.hh"
#include "IvyFramework/IvyDataTools/interface/BaseTree.h"
#include "IvyFramework/IvyDataTools/interface/SimpleEntry.h"
#include "GlobalCollectionNames.h"
#include "RunLumiEventBlock.h"
#include "SamplesCore.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "IsotrackSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "ParticleDisambiguator.h"
#include "MuonHandler.h"
#include "ElectronHandler.h"
#include "JetMETHandler.h"
#include "EventFilterHandler.h"
#include "SimEventHandler.h"
#include "GenInfoHandler.h"
#include "IsotrackHandler.h"
#include "SamplesCore.h"
#include "FourTopTriggerHelpers.h"
#include "DileptonHandler.h"
#include "SplitFileAndAddForTransfer.h"


using namespace std;
using namespace HelperFunctions;
using namespace IvyStreamHelpers;


struct SelectionTracker{
  std::vector<TString> ordered_reqs;
  std::unordered_map<TString, std::pair<double, double>> req_sumws_pair_map;

  void accumulate(TString const& strsel, double const& wgt);
  void print() const;
};
void SelectionTracker::accumulate(TString const& strsel, double const& wgt){
  if (!HelperFunctions::checkListVariable(ordered_reqs, strsel)){
    req_sumws_pair_map[strsel] = std::pair<double, double>(0, 0);
    ordered_reqs.push_back(strsel);
  }
  auto it_req_sumws_pair = req_sumws_pair_map.find(strsel);
  it_req_sumws_pair->second.first += wgt;
  it_req_sumws_pair->second.second += wgt*wgt;
}
void SelectionTracker::print() const{
  IVYout << "Selection summary:" << endl;
  for (auto const& strsel:ordered_reqs){
    auto it_req_sumws_pair = req_sumws_pair_map.find(strsel);
    IVYout << "\t- " << strsel << ": " << setprecision(15) << it_req_sumws_pair->second.first << " +- " << std::sqrt(it_req_sumws_pair->second.second) << endl;
  }
}

int ScanChain(std::string const& strdate, std::string const& dset, std::string const& proc, double const& xsec, SimpleEntry const& extra_arguments){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  TDirectory* curdir = gDirectory;

  // This is the output directory.
  // Output should always be recorded as if you are running the job locally.
  // We will inform the Condor job later on that some files would need transfer if we are running on Condor.
  TString coutput_main = ANALYSISPKGPATH + "test/output/Analysis_CutBased/" + strdate.data() + "/" + SampleHelpers::getDataPeriod();
  HostHelpers::ExpandEnvironmentVariables(coutput_main);
  gSystem->mkdir(coutput_main, true);
  TString stroutput = coutput_main + "/" + proc.data() + ".root"; // This is the output file.
  TString stroutput_log = coutput_main + "/log_" + proc.data() + ".out"; // This is the output file.
  IVYout.open(stroutput_log.Data());

  constexpr bool cleanJetsFromFakeableObjects = true;
  ParticleSelectionHelpers::setUseFakeableIdForJetCleaning(cleanJetsFromFakeableObjects);

  float const absEtaThr_ak4jets = (SampleHelpers::getDataYear()<=2016 ? AK4JetSelectionHelpers::etaThr_btag_Phase0Tracker : AK4JetSelectionHelpers::etaThr_btag_Phase1Tracker);

  // Turn on synchronization exercise options
  std::string input_files;
  extra_arguments.getNamedVal("input_files", input_files);
  bool runSyncExercise = false;
  extra_arguments.getNamedVal("run_sync", runSyncExercise);
  bool writeSyncObjects = false;
  bool forceApplyPreselection = false;
  if (runSyncExercise){
    extra_arguments.getNamedVal("write_sync_objects", writeSyncObjects);
    extra_arguments.getNamedVal("force_sync_preselection", forceApplyPreselection);
  }

  // Shorthand option for the Run 2 UL analysis proposal
  bool use_shorthand_Run2_UL_proposal_config;
  extra_arguments.getNamedVal("shorthand_Run2_UL_proposal_config", use_shorthand_Run2_UL_proposal_config);

  // Options to set alternative muon and electron IDs, or b-tagging WP
  std::string muon_id_name;
  std::string electron_id_name;
  std::string btag_WPname;
  if (use_shorthand_Run2_UL_proposal_config){
    muon_id_name = electron_id_name = "TopMVA_Run2";
    btag_WPname = "loose";
  }
  else{
    extra_arguments.getNamedVal("muon_id", muon_id_name);
    extra_arguments.getNamedVal("electron_id", electron_id_name);
    extra_arguments.getNamedVal("btag", btag_WPname);
  }

  if (muon_id_name!=""){
    IVYout << "Switching to muon id " << muon_id_name << "..." << endl;
    MuonSelectionHelpers::setSelectionTypeByName(muon_id_name);
  }
  else IVYout << "Using default muon id = " << MuonSelectionHelpers::selection_type << "..." << endl;

  if (electron_id_name!=""){
    IVYout << "Switching to electron id " << electron_id_name << "..." << endl;
    ElectronSelectionHelpers::setSelectionTypeByName(electron_id_name);
  }
  else IVYout << "Using default electron id = " << ElectronSelectionHelpers::selection_type << "..." << endl;

  AK4JetSelectionHelpers::SelectionBits bit_preselection_btag = AK4JetSelectionHelpers::kPreselectionTight_BTagged_Medium;
  if (btag_WPname!=""){
    std::string btag_WPname_lower;
    HelperFunctions::lowercase(btag_WPname, btag_WPname_lower);
    IVYout << "Switching to b-tagging WP " << btag_WPname_lower << "..." << endl;
    if (btag_WPname_lower=="loose") bit_preselection_btag = AK4JetSelectionHelpers::kPreselectionTight_BTagged_Loose;
    else if (btag_WPname_lower=="medium") bit_preselection_btag = AK4JetSelectionHelpers::kPreselectionTight_BTagged_Medium;
    else if (btag_WPname_lower=="tight") bit_preselection_btag = AK4JetSelectionHelpers::kPreselectionTight_BTagged_Tight;
    else{
      IVYerr << "btag=" << btag_WPname << " is not implemented." << endl;
      assert(0);
    }
  }
  else IVYout << "Using default b-tagging WP = " << static_cast<int>(bit_preselection_btag)-static_cast<int>(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Loose) << "..." << endl;

  // Flag to control whether any preselection other than nleps>=2 to be applied
  bool const applyPreselection = !runSyncExercise || forceApplyPreselection;

  // Trigger configuration
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton{
    TriggerHelpers::kDoubleMu,
    TriggerHelpers::kDoubleEle,
    TriggerHelpers::kMuEle
  };
  // These PFHT triggers were used in the 2016 analysis. They are a bit more efficient.
  if (SampleHelpers::getDataYear()==2016) requiredTriggers_Dilepton = std::vector<TriggerHelpers::TriggerType>{
      TriggerHelpers::kDoubleMu_PFHT,
      TriggerHelpers::kDoubleEle_PFHT,
      TriggerHelpers::kMuEle_PFHT
  };
  std::vector<std::string> const hltnames_Dilepton = TriggerHelpers::getHLTMenus(requiredTriggers_Dilepton);

  // Related to triggers is how we apply loose and fakeable IDs in electrons.
  // This setting should vary for 2016 when analyzing fake rates instead of the signal or SR-like control regions.
  if (ElectronSelectionHelpers::selection_type == ElectronSelectionHelpers::kCutbased_Run2) ElectronSelectionHelpers::setApplyMVALooseFakeableNoIsoWPs((SampleHelpers::getDataYear()==2016));

  // Declare handlers
  GenInfoHandler genInfoHandler;
  SimEventHandler simEventHandler;
  EventFilterHandler eventFilter(requiredTriggers_Dilepton);
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;

  // These are called handlers, but they are more like helpers.
  DileptonHandler dileptonHandler;
  ParticleDisambiguator particleDisambiguator;

  // Some advanced event filters
  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(true);
  eventFilter.setCheckHLTPathRunRanges(true);

  curdir->cd();

  // Acquire input tree/chain
  TString strinput = SampleHelpers::getInputDirectory() + "/" + SampleHelpers::getDataPeriod() + "/" + proc.data();
  TString cinput = (input_files=="" ? strinput + "/*.root" : strinput + "/" + input_files.data());
  IVYout << "Accessing input files " << cinput << "..." << endl;
  BaseTree* tin = new BaseTree(cinput, "Events", "", "");
  tin->sampleIdentifier = SampleHelpers::getSampleIdentifier(dset);
  bool const isData = SampleHelpers::checkSampleIsData(tin->sampleIdentifier);
  if (!isData && xsec<0.){
    IVYerr << "xsec = " << xsec << " is not valid." << endl;
    assert(0);
  }

  double sum_wgts = (isData ? 1 : 0);
  if (!isData){
    for (auto const& fname:SampleHelpers::lsdir(strinput.Data())){
      if (input_files!="" && fname!=input_files.data()) continue;
      if (fname.EndsWith(".root")){
        TFile* ftmp = TFile::Open(strinput + "/" + fname, "read");
        TH2D* hCounters = (TH2D*) ftmp->Get("Counters");
        sum_wgts += hCounters->GetBinContent(1, 1);
        ftmp->Close();
      }
    }
  }
  if (sum_wgts==0.){
    IVYerr << "Sum of pre-recorded weights cannot be zero." << endl;
    assert(0);
  }

  curdir->cd();

  // Calculate the overall normalization scale on the events.
  // Includes xsec (in fb), lumi (in fb-1), and 1/sum of weights in all of the MC.
  // Data normalizaion factor is always 1.
  double const lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  double norm_scale = (isData ? 1. : xsec * xsecScale * lumi)/sum_wgts;

  IVYout << "Valid data periods for " << SampleHelpers::getDataPeriod() << ": " << SampleHelpers::getValidDataPeriods() << endl;
  IVYout << "Integrated luminosity: " << lumi << endl;
  IVYout << "Acquired a sum of weights of " << sum_wgts << ". Overall normalization will be " << norm_scale << "." << endl;

  curdir->cd();

  // Wrap the ivies around the input tree:
  // Booking is basically SetBranchStatus+SetBranchAddress. You can book for as many trees as you would like.
  // In some cases, bookBranches also informs the ivy dynamically that it is supposed to consume certain entries.
  // For entries common to all years or any data or MC, the consumption information is handled in the ivy constructor already.
  // None of these mean the ivy establishes its access to the input tree yet.
  // Wrapping a tree informs the ivy that it is supposed to consume the booked entries from that particular tree.
  // Without wrapping, you are not really accessing the entries from the input tree to construct the physics objects;
  // all you would get are 0 electrons, 0 jets, everything failing event filters etc.
  genInfoHandler.bookBranches(tin);
  genInfoHandler.wrapTree(tin);

  simEventHandler.bookBranches(tin);
  simEventHandler.wrapTree(tin);

  eventFilter.bookBranches(tin);
  eventFilter.wrapTree(tin);

  muonHandler.bookBranches(tin);
  muonHandler.wrapTree(tin);

  electronHandler.bookBranches(tin);
  electronHandler.wrapTree(tin);

  jetHandler.bookBranches(tin);
  jetHandler.wrapTree(tin);

  isotrackHandler.bookBranches(tin);
  isotrackHandler.wrapTree(tin);

  EventNumber_t* ptr_EventNumber = nullptr;
  tin->bookBranch<EventNumber_t>("event", 0);
  tin->getValRef("event", ptr_EventNumber);

  curdir->cd();

  TFile* foutput_sync_objects = nullptr;
  BaseTree* tout_sync_objects = nullptr;
  SimpleEntry rcd_sync_objects;
  foutput_sync_objects = TFile::Open(coutput_main + "/" + Form("sync_objects_%s.root", proc.data()), "recreate");
  tout_sync_objects = new BaseTree("SyncObjects");
  tout_sync_objects->setAutoSave(0);
  curdir->cd();
  

  // Keep track of sums of weights
  SelectionTracker seltracker;

  bool firstOutputEvent = true;
  unsigned int n_traversed = 0;
  unsigned int n_recorded = 0;
  int nEntries = tin->getNEvents();
  IVYout << "Looping over " << nEntries << " events..." << endl;
  for (int ev=0; ev<nEntries; ev++){
    if (SampleHelpers::doSignalInterrupt==1) break;

    tin->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);
    n_traversed++;

    genInfoHandler.constructGenInfo();
    auto const& genInfo = genInfoHandler.getGenInfo();

    simEventHandler.constructSimEvent();

    double wgt = 1;
    if (!isData){
      double genwgt = 1;
      genwgt = genInfo->getGenWeight(SystematicsHelpers::sNominal);

      double puwgt = 1;
      puwgt = simEventHandler.getPileUpWeight(SystematicsHelpers::sNominal);

      wgt = genwgt * puwgt;

      // Add L1 prefiring weight for 2016 and 2017
      wgt *= simEventHandler.getL1PrefiringWeight(SystematicsHelpers::sNominal);
    }

    muonHandler.constructMuons();
    electronHandler.constructElectrons();
    jetHandler.constructJetMET(&simEventHandler);

    // !!!IMPORTANT!!!
    // NEVER USE LEPTONS AND JETS IN AN ANALYSIS BEFORE DISAMBIGUATING THEM!
    // Muon and electron handlers do not apply any selection, so the selection bits are all 0.
    // In order to compute pTratio and pTrel, you need jets!
    // ParticleDisambiguator does the matching, and assigns the overlapping jets (or closest ones) as 'mothers' of the leptons.
    // Once mothers are assigned, ParticleObject::ptratio and ptrel functions work as intended,
    // and you can apply the additional selections on these variables this way.
    // ParticleDisambiguator then cleans all geometrically overlapping jets by resetting their selection bits, which makes them unusable.
    particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, nullptr, &jetHandler);

    bool const printObjInfo = runSyncExercise
      &&
      HelperFunctions::checkListVariable(std::vector<int>{ 6994, 7046, 11794 }, ev);
      //HelperFunctions::checkListVariable(std::vector<int>{ 1233, 1475, 1546, 1696, 2011, 2103, 2801, 2922, 3378, 3407, 3575, 3645, 5021, 5127, 6994, 7000, 7046, 7341, 7351, 8050, 9931, 10063, 10390, 10423, 10623, 10691, 10791, 10796, 11127, 11141, 11279, 11794, 12231, 12996, 13115, 13294, 13550, 14002, 14319, 15062, 15754, 16153, 16166, 16316, 16896, 16911, 17164 }, ev);
      //HelperFunctions::checkListVariable(std::vector<int>{663, 1469, 3087, 3281}, ev);
      //HelperFunctions::checkListVariable(std::vector<int>{204, 353, 438, 1419}, ev);
      //HelperFunctions::checkListVariable(std::vector<int>{3, 15, 30, 31, 32, 41, 153, 154, 162, 197, 215, 284, 317, 360, 572, 615, 747, 1019, 1119, 1129}, ev);

    // Muon sync. write variables
#define SYNC_MUONS_BRANCH_VECTOR_COMMANDS \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_loose) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_fakeable) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_tight) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, pt) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, eta) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, phi) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, mass) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, ptrel_final) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, ptratio_final) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, bscore) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, muons, extMVAscore) \
    MUON_EXTRA_VARIABLES
    // Electron sync. write variables
#define SYNC_ELECTRONS_BRANCH_VECTOR_COMMANDS \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_loose) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_fakeable) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_tight) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, pt) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, eta) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, etaSC) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, phi) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, mass) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, ptrel_final) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, ptratio_final) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, mvaFall17V2noIso_raw) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, bscore) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, electrons, extMVAscore) \
    ELECTRON_EXTRA_VARIABLES
    // ak4jet sync. write variables
#define SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, ak4jets, is_tight) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, ak4jets, is_btagged) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, ak4jets, is_clean) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, pt) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, eta) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, phi) \
    SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, mass) \
    AK4JET_EXTRA_INPUT_VARIABLES
    // All sync. write objects
#define SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS \
   SYNC_MUONS_BRANCH_VECTOR_COMMANDS \
   SYNC_ELECTRONS_BRANCH_VECTOR_COMMANDS \
   SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) std::vector<TYPE> COLLNAME##_##NAME;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, muons, NAME)
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, electrons, NAME)
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, ak4jets, NAME)
    SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef ELECTRON_VARIABLE
#undef MUON_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND

    if (printObjInfo) IVYout << "Lepton info for event " << ev << ":" << endl;

    auto const& muons = muonHandler.getProducts();
    std::vector<MuonObject*> muons_selected;
    std::vector<MuonObject*> muons_tight;
    std::vector<MuonObject*> muons_fakeable;
    std::vector<MuonObject*> muons_loose;
    for (auto const& part:muons){
      float pt = part->pt();
      float eta = part->eta();
      float phi = part->phi();
      float mass = part->mass();

      bool is_tight = false;
      bool is_fakeable = false;
      bool is_loose = false;

      if (ParticleSelectionHelpers::isTightParticle(part)){
        muons_tight.push_back(part);
        is_loose = is_fakeable = is_tight = true;
      }
      else if (ParticleSelectionHelpers::isFakeableParticle(part)){
        muons_fakeable.push_back(part);
        is_loose = is_fakeable = true;
      }
      else if (ParticleSelectionHelpers::isLooseParticle(part)){
        muons_loose.push_back(part);
        is_loose = true;
      }

      float ptrel_final = part->ptrel();
      float ptratio_final = part->ptratio();

      float extMVAscore=-99;
      bool has_extMVAscore = part->getExternalMVAScore(MuonSelectionHelpers::selection_type, extMVAscore);

      float bscore = 0;
      AK4JetObject* mother = nullptr;
      for (auto const& mom:part->getMothers()){
        mother = dynamic_cast<AK4JetObject*>(mom);
        if (mother) break;
      }
      if (mother) bscore = mother->extras.btagDeepFlavB;

      if (writeSyncObjects){
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.push_back(NAME);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) muons_##NAME.push_back(part->extras.NAME);
        SYNC_MUONS_BRANCH_VECTOR_COMMANDS;
#undef MUON_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND
      }
    }
    HelperFunctions::appendVector(muons_selected, muons_tight);
    HelperFunctions::appendVector(muons_selected, muons_fakeable);
    HelperFunctions::appendVector(muons_selected, muons_loose);

    auto const& electrons = electronHandler.getProducts();
    std::vector<ElectronObject*> electrons_selected;
    std::vector<ElectronObject*> electrons_tight;
    std::vector<ElectronObject*> electrons_fakeable;
    std::vector<ElectronObject*> electrons_loose;
    for (auto const& part:electrons){
      float pt = part->pt();
      float eta = part->eta();
      float etaSC = part->etaSC();
      float phi = part->phi();
      float mass = part->mass();

      bool is_tight = false;
      bool is_fakeable = false;
      bool is_loose = false;

      if (ParticleSelectionHelpers::isTightParticle(part)){
        electrons_tight.push_back(part);
        is_loose = is_fakeable = is_tight = true;
      }
      else if (ParticleSelectionHelpers::isFakeableParticle(part)){
        electrons_fakeable.push_back(part);
        is_loose = is_fakeable = true;
      }
      else if (ParticleSelectionHelpers::isLooseParticle(part)){
        electrons_loose.push_back(part);
        is_loose = true;
      }

      float mvaFall17V2noIso_raw = 0.5 * std::log((1. + part->extras.mvaFall17V2noIso)/(1. - part->extras.mvaFall17V2noIso));
      float ptrel_final = part->ptrel();
      float ptratio_final = part->ptratio();

      float extMVAscore=-99;
      bool has_extMVAscore = part->getExternalMVAScore(ElectronSelectionHelpers::selection_type, extMVAscore);

      float bscore = 0;
      AK4JetObject* mother = nullptr;
      for (auto const& mom:part->getMothers()){
        mother = dynamic_cast<AK4JetObject*>(mom);
        if (mother) break;
      }
      if (mother) bscore = mother->extras.btagDeepFlavB;

      if (writeSyncObjects){
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.push_back(NAME);
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) electrons_##NAME.push_back(part->extras.NAME);
        SYNC_ELECTRONS_BRANCH_VECTOR_COMMANDS;
#undef ELECTRON_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND
      }
    }
    HelperFunctions::appendVector(electrons_selected, electrons_tight);
    HelperFunctions::appendVector(electrons_selected, electrons_fakeable);
    HelperFunctions::appendVector(electrons_selected, electrons_loose);

    unsigned int const nleptons_tight = muons_tight.size() + electrons_tight.size();
    unsigned int const nleptons_fakeable = muons_fakeable.size() + electrons_fakeable.size();
    unsigned int const nleptons_loose = muons_loose.size() + electrons_loose.size();
    unsigned int const nleptons_selected = nleptons_tight + nleptons_fakeable + nleptons_loose;

    std::vector<ParticleObject*> leptons_selected; leptons_selected.reserve(nleptons_selected);
    for (auto const& part:muons_selected) leptons_selected.push_back(dynamic_cast<ParticleObject*>(part));
    for (auto const& part:electrons_selected) leptons_selected.push_back(dynamic_cast<ParticleObject*>(part));
    ParticleObjectHelpers::sortByGreaterPt(leptons_selected);

    std::vector<ParticleObject*> leptons_tight; leptons_tight.reserve(nleptons_tight);
    for (auto const& part:muons_tight) leptons_tight.push_back(dynamic_cast<ParticleObject*>(part));
    for (auto const& part:electrons_tight) leptons_tight.push_back(dynamic_cast<ParticleObject*>(part));
    ParticleObjectHelpers::sortByGreaterPt(leptons_tight);

    double ak4jets_pt40_HT=0;
    auto const& ak4jets = jetHandler.getAK4Jets();
    std::vector<AK4JetObject*> ak4jets_tight_pt40;
    std::vector<AK4JetObject*> ak4jets_tight_pt25_btagged;
    for (auto const& jet:ak4jets){
      float pt = jet->pt();
      float eta = jet->eta();
      float phi = jet->phi();
      float mass = jet->mass();

      bool is_tight = ParticleSelectionHelpers::isTightJet(jet);
      bool is_btagged = jet->testSelectionBit(bit_preselection_btag);
      constexpr bool is_clean = true;

      if (writeSyncObjects){
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.push_back(NAME);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) ak4jets_##NAME.push_back(jet->extras.NAME);
        SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND
      }

      if (is_tight && pt>=25. && std::abs(eta)<absEtaThr_ak4jets){
        if (is_btagged) ak4jets_tight_pt25_btagged.push_back(jet);
        if (pt>=40.){
          ak4jets_tight_pt40.push_back(jet);
          ak4jets_pt40_HT += pt;
        }
      }
    }

    unsigned int const nak4jets_tight_pt40 = ak4jets_tight_pt40.size();
    unsigned int const nak4jets_tight_pt25_btagged = ak4jets_tight_pt25_btagged.size();

    auto const& eventmet = jetHandler.getPFMET();


    // BEGIN PRESELECTION
    seltracker.accumulate("Full sample", wgt);


    if (nleptons_fakeable!=1) continue; // Skims are required to apply this selection, so no additional test on applyPreselection.


    // Put event filters to the last because data has unique event tracking enabled.
    eventFilter.constructFilters(&simEventHandler);
    //if (!eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets)) continue; // Test for 2018 partial HEM failure
    //if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons, nullptr, nullptr)) continue; // Test for 2018 partial HEM failure
    //seltracker.accumulate("Pass HEM veto", wgt);
    if (!eventFilter.passMETFilters()) continue; // Test for MET filters
    seltracker.accumulate("Pass MET filters", wgt);
    if (!eventFilter.isUniqueDataEvent()) continue; // Test if the data event is unique (i.e., dorky). Does not do anything in the MC.
    seltracker.accumulate("Pass unique event check", wgt);

    // TODO: get code from ulascan
    // does the trigger selection and records the flag for the control trigger


    /*************************************************/
    /* NO MORE CALLS TO SELECTION BEYOND THIS POINT! */
    /*************************************************/
    // Write object sync. info.
    if (writeSyncObjects){
      rcd_sync_objects.setNamedVal("EventNumber", *ptr_EventNumber);
      rcd_sync_objects.setNamedVal("PFMET_pt_final", eventmet->pt());
      rcd_sync_objects.setNamedVal("PFMET_phi_final", eventmet->phi());

#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) rcd_sync_objects.setNamedVal(Form("%s_%s", #COLLNAME, #NAME), COLLNAME##_##NAME);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, muons, NAME)
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, electrons, NAME)
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, ak4jets, NAME)
      SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef ELECTRON_VARIABLE
#undef MUON_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND

      if (firstOutputEvent){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_sync_objects.named##name_t##s.begin(); itb!=rcd_sync_objects.named##name_t##s.end(); itb++) tout_sync_objects->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_sync_objects.namedV##name_t##s.begin(); itb!=rcd_sync_objects.namedV##name_t##s.end(); itb++) tout_sync_objects->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_sync_objects.namedVV##name_t##s.begin(); itb!=rcd_sync_objects.namedVV##name_t##s.end(); itb++) tout_sync_objects->putBranch(itb->first, &(itb->second));
        SIMPLE_DATA_OUTPUT_DIRECTIVES;
        VECTOR_DATA_OUTPUT_DIRECTIVES;
        DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES;
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
      }

      // Record whatever is in rcd_sync_objects into the tree.
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_sync_objects.named##name_t##s.begin(); itb!=rcd_sync_objects.named##name_t##s.end(); itb++) tout_sync_objects->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_sync_objects.namedV##name_t##s.begin(); itb!=rcd_sync_objects.namedV##name_t##s.end(); itb++) tout_sync_objects->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_sync_objects.namedVV##name_t##s.begin(); itb!=rcd_sync_objects.namedVV##name_t##s.end(); itb++) tout_sync_objects->setVal(itb->first, &(itb->second));
      SIMPLE_DATA_OUTPUT_DIRECTIVES;
      VECTOR_DATA_OUTPUT_DIRECTIVES;
      DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES;
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

      tout_sync_objects->fill();
    }

    n_recorded++;

    if (firstOutputEvent) firstOutputEvent = false;
  }

  IVYout << "Number of events recorded: " << n_recorded << " / " << n_traversed << " / " << nEntries << endl;
  seltracker.print();

  tout_sync_objects->writeToFile(foutput_sync_objects);
  delete tout_sync_objects;

  foutput_sync_objects->Close();

  if (runSyncExercise){
    SampleHelpers::splitFileAndAddForTransfer(stroutput_sync);
  }

  curdir->cd();
  delete tin;

  IVYout.close();

  // Split large files, and add them to the transfer queue from Condor to the target site
  // Does nothing if you are running the program locally because your output is already in the desired location.
  SampleHelpers::splitFileAndAddForTransfer(stroutput);
  return 0;
}

int main(int argc, char** argv){
  constexpr int iarg_offset=1; // argv[0]==[Executable name]

  bool print_help=false, has_help=false;
  std::string str_dset;
  std::string str_proc;
  std::string str_period;
  std::string str_tag;
  std::string str_outtag;
  SimpleEntry extra_arguments;
  double xsec = -1;
  for (int iarg=iarg_offset; iarg<argc; iarg++){
    std::string strarg = argv[iarg];
    std::string wish, value;
    splitOption(strarg, wish, value, '=');

    if (wish.empty()){
      if (value=="help"){ print_help=has_help=true; }
      else if (value=="shorthand_Run2_UL_proposal_config") extra_arguments.setNamedVal<bool>(value, true);
      else{
        IVYerr << "ERROR: Unknown argument " << value << endl;
        print_help=true;
      }
    }
    else if (wish=="dataset") str_dset = value;
    else if (wish=="short_name") str_proc = value;
    else if (wish=="period") str_period = value;
    else if (wish=="input_tag") str_tag = value;
    else if (wish=="output_tag") str_outtag = value;
    else if (wish=="input_files"){
      if (value.find("/")!=std::string::npos){
        IVYerr << "ERROR: Input file specification cannot contain directory structure." << endl;
        print_help=true;
      }
      extra_arguments.setNamedVal(wish, value);
    }
    else if (wish=="run_sync" || wish=="write_sync_objects" || wish=="force_sync_preselection"){
      bool tmpval;
      HelperFunctions::castStringToValue(value, tmpval);
      extra_arguments.setNamedVal(wish, tmpval);
    }
    else if (wish=="muon_id" || wish=="electron_id" || wish=="btag") extra_arguments.setNamedVal(wish, value);
    else if (wish=="xsec"){
      if (xsec<0.) xsec = 1;
      xsec *= std::stod(value);
    }
    else if (wish=="BR"){
      if (xsec<0.) xsec = 1;
      xsec *= std::stod(value);
    }
    else{
      IVYerr << "ERROR: Unknown argument " << wish << "=" << value << endl;
      print_help=true;
    }
  }

  if (!print_help && (str_proc=="" || str_dset=="" || str_period=="" || str_tag=="" || str_outtag=="")){
    IVYerr << "ERROR: Not all mandatory inputs are present." << endl;
    print_help=true;
  }

  if (print_help){
    IVYout << "skim_UL options:\n\n";
    IVYout << "- help: Prints this help message.\n";
    IVYout << "- dataset: Data set name. Mandatory.\n";
    IVYout << "- short_name: Process short name. Mandatory.\n";
    IVYout << "- period: Data period. Mandatory.\n";
    IVYout << "- input_tag: Version of the input skims. Mandatory.\n";
    IVYout << "- output_tag: Version of the output. Mandatory.\n";
    IVYout << "- xsec: Cross section value. Mandatory in the MC.\n";
    IVYout << "- BR: BR value. Mandatory in the MC.\n";
    IVYout << "- input_files: Input files to run. Optional. Default is to run on all files.\n";
    IVYout << "- run_sync: Turn on synchronization output. Optional. Default is to run without synchronization output.\n";
    IVYout << "- write_sync_objects: Create a file that contains the info. for all leptons and jets, and event identifiers. Ignored if run_sync=false. Optional. Default is to not produce such a file.\n";
    IVYout << "- force_sync_preselection: When sync. mode is on, also apply SR/CR preselection.  Ignored if run_sync=false. Optional. Default is to run without preselection.\n";
    IVYout
      << "- shorthand_Run2_UL_proposal_config: Shorthand flag for the switches for the Run 2 UL analysis proposal:\n"
      << "  * muon_id='TopMVA_Run2'\n"
      << "  * electron_id='TopMVA_Run2'\n"
      << "  * btag='loose'\n"
      << "  The use of this shorthand will ignore the user-defined setting of these options above.\n";
    IVYout
      << "- muon_id: Can be 'Cutbased_Run2', 'TopMVA_Run2', or 'TopMVAv2_Run2'.\n"
      << "  Default is whatever is in MuonSelectionHelpers (currently 'Cutbased_Run2') if no value is given.\n";
    IVYout
      << "- electron_id: Can be 'Cutbased_Run2', 'TopMVA_Run2', or 'TopMVAv2_Run2'.\n"
      << "  Default is whatever is in ElectronSelectionHelpers (currently 'Cutbased_Run2') if no value is given.\n";
    IVYout << "- btag: Name of the b-tagging WP. Can be 'medium' or 'loose' (case-insensitive). Default='medium'.\n";

    IVYout << endl;
    return (has_help ? 0 : 1);
  }

  SampleHelpers::configure(str_period, Form("skims:%s", str_tag.data()), HostHelpers::kUCSDT2);

  return ScanChain(str_outtag, str_dset, str_proc, xsec, extra_arguments);
}