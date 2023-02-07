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
#include "BtagScaleFactorHandler.h"
#include "SamplesCore.h"
#include "FourTopTriggerHelpers.h"
#include "DileptonHandler.h"
#include "InputChunkSplitter.h"
#include "SplitFileAndAddForTransfer.h"


using namespace std;
using namespace HelperFunctions;
using namespace IvyStreamHelpers;


struct SelectionTracker{
  std::vector<TString> ordered_reqs;
  std::unordered_map<TString, std::pair<double, double>> req_sumws_pair_map;

  void accumulate(TString const& strsel, double const& wgt, bool printPass=false);
  void print() const;
};
void SelectionTracker::accumulate(TString const& strsel, double const& wgt, bool printPass){
  if (!HelperFunctions::checkListVariable(ordered_reqs, strsel)){
    req_sumws_pair_map[strsel] = std::pair<double, double>(0, 0);
    ordered_reqs.push_back(strsel);
  }
  auto it_req_sumws_pair = req_sumws_pair_map.find(strsel);
  it_req_sumws_pair->second.first += wgt;
  it_req_sumws_pair->second.second += wgt*wgt;
  if (printPass && wgt!=0.) IVYout << strsel << endl;
}
void SelectionTracker::print() const{
  IVYout << "Selection summary:" << endl;
  for (auto const& strsel:ordered_reqs){
    auto it_req_sumws_pair = req_sumws_pair_map.find(strsel);
    IVYout << "\t- " << strsel << ": " << setprecision(15) << it_req_sumws_pair->second.first << " +- " << std::sqrt(it_req_sumws_pair->second.second) << endl;
  }
}

int ScanChain(std::string const& strdate, std::string const& dset, std::string const& proc, double const& xsec, int const& ichunk, int const& nchunks, SimpleEntry const& extra_arguments){
  bool const isCondorRun = SampleHelpers::checkRunOnCondor();
  if (!isCondorRun) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  TDirectory* curdir = gDirectory;

  // Systematics hypothesis
  // FIXME: Should be generalized to include all systematics at once.
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  // Data period quantities
  auto const& theDataPeriod = SampleHelpers::getDataPeriod();
  auto const& theDataYear = SampleHelpers::getDataYear();

  // Configure analysis-specific stuff
  bool no_FOs_phys_checks = false;
  extra_arguments.getNamedVal("no_FOs_phys_checks", no_FOs_phys_checks);
  bool const useFakeableIdForPhysicsChecks = !no_FOs_phys_checks;
  ParticleSelectionHelpers::setUseFakeableIdForPhysicsChecks(useFakeableIdForPhysicsChecks);
  IVYout << "Configured to use " << (useFakeableIdForPhysicsChecks ? "fakeable" : "only tight") << " leptons in physics objects checks..." << endl;

  bool only_tight_dileptons = false;
  extra_arguments.getNamedVal("only_tight_dileptons", only_tight_dileptons);
  IVYout << "Configured to use " << (only_tight_dileptons ? "only tight" : "loose, fakeable, or tight") << " leptons in pairing..." << endl;

  bool ignore_PU_rewgt = false;
  extra_arguments.getNamedVal("ignore_PU_rewgt", ignore_PU_rewgt);
  if (ignore_PU_rewgt) IVYout << "The PU weight will be ignored from the 'event_wgt' variable and will be recorded separately instead." << endl;

  bool ignore_L1Prefiring_rewgt = false;
  extra_arguments.getNamedVal("ignore_L1Prefiring_rewgt", ignore_L1Prefiring_rewgt);
  if (ignore_L1Prefiring_rewgt) IVYout << "The L1 prefiring weight will be ignored from the 'event_wgt' variable and will be recorded separately instead." << endl;

  constexpr bool useIsotrackVeto = false;

  float const absEtaThr_ak4jets = (theDataYear<=2016 ? AK4JetSelectionHelpers::etaThr_btag_Phase0Tracker : AK4JetSelectionHelpers::etaThr_btag_Phase1Tracker);

  double const lumi = SampleHelpers::getIntegratedLuminosity(theDataPeriod);
  IVYout << "Valid data periods for " << theDataPeriod << ": " << SampleHelpers::getValidDataPeriods() << endl;
  IVYout << "Integrated luminosity: " << lumi << endl;

  // This is the output directory.
  // Output should always be recorded as if you are running the job locally.
  // If we are running on Condor, ee will inform the batch job later on that some files would need transfer.
  TString coutput_main = TString("output/Analysis_CutBased/") + strdate.data() + "/" + theDataPeriod;
  if (!isCondorRun) coutput_main = ANALYSISPKGPATH + "/test/" + coutput_main;
  HostHelpers::ExpandEnvironmentVariables(coutput_main);
  gSystem->mkdir(coutput_main, true);

  std::vector<std::pair<std::string, std::string>> dset_proc_pairs;
  {
    std::vector<std::string> dsets, procs;
    HelperFunctions::splitOptionRecursive(dset, dsets, ',', false);
    HelperFunctions::splitOptionRecursive(proc, procs, ',', false);
    HelperFunctions::zipVectors(dsets, procs, dset_proc_pairs);
  }
  bool const has_multiple_dsets = dset_proc_pairs.size()>1;

  // Configure output file name
  std::string output_file;
  extra_arguments.getNamedVal("output_file", output_file);
  if (output_file==""){
    if (has_multiple_dsets){
      IVYerr << "Multiple data sets are passed. An output file identifier must be specified through the option 'output_file'." << endl;
      assert(0);
    }
    output_file = proc;
  }
  if (nchunks>0) output_file = Form("%s_%i_of_%i", output_file.data(), ichunk, nchunks);

  // Configure full output file names with path
  TString stroutput = coutput_main + "/" + output_file.data() + ".root"; // This is the output ROOT file.
  TString stroutput_log = coutput_main + "/log_" + output_file.data() + ".out"; // This is the output log file.
  TString stroutput_err = coutput_main + "/log_" + output_file.data() + ".err"; // This is the error log file.
  IVYout.open(stroutput_log.Data());
  IVYerr.open(stroutput_err.Data());

  // Produce output trees instead of histograms
  // Note: With this switch is on, keeping track of the sums of weights becomes meaningless.
  bool produce_trees = false;
  extra_arguments.getNamedVal("produce_trees", produce_trees);
  if (produce_trees) IVYout << "Producing output trees instead of histograms..." << endl;

  // Turn on synchronization exercise options
  std::string input_files;
  extra_arguments.getNamedVal("input_files", input_files);
  bool runSyncExercise = false;
  if (!has_multiple_dsets) extra_arguments.getNamedVal("run_sync", runSyncExercise);
  bool writeSyncObjects = false;
  bool forceApplyPreselection = false;
  std::string sync_evtno_filter_input;
  if (runSyncExercise){
    extra_arguments.getNamedVal("write_sync_objects", writeSyncObjects);
    extra_arguments.getNamedVal("force_sync_preselection", forceApplyPreselection);
    extra_arguments.getNamedVal("sync_evtno_filter_input", sync_evtno_filter_input);
    IVYout << "Sync options:" << endl;
    IVYout << "\t- write_sync_objects: " << writeSyncObjects << endl;
    IVYout << "\t- force_sync_preselection: " << forceApplyPreselection << endl;
    IVYout << "\t- Event number filter file: " << sync_evtno_filter_input << endl;
  }

  if (produce_trees && runSyncExercise){
    IVYerr << "Tree output and sync exercise should not be run at the same time." << endl;
    assert(0);
  }

  // Create the sync files
  ofstream foutput_sync;
  TString stroutput_sync;
  if (runSyncExercise){
    stroutput_sync = coutput_main + "/" + Form("sync_%s.csv", output_file.data());
    foutput_sync.open(stroutput_sync.Data(), std::ios_base::out);
    foutput_sync << "Event#,#inFile,MET,MET_phi,Nlooseleptons,Nfakeableleptons,Ntightleptons,Goodsspair?,HT,Njets,Nbjets,SR/CR" << endl;
  }

  TFile* foutput_sync_objects = nullptr;
  BaseTree* tout_sync_objects = nullptr;
  SimpleEntry rcd_sync_objects;
  if (writeSyncObjects){
    foutput_sync_objects = TFile::Open(coutput_main + "/" + Form("sync_objects_%s.root", output_file.data()), "recreate");
    tout_sync_objects = new BaseTree("SyncObjects");
    tout_sync_objects->setAutoSave(0);
    curdir->cd();
  }

  std::vector<EventNumber_t> eventnumber_filter;
  if (sync_evtno_filter_input!=""){
    if (sync_evtno_filter_input.find(".root")!=std::string::npos){
      std::string sync_evtno_filter_input_fname = sync_evtno_filter_input;
      std::string sync_evtno_filter_input_treename = "SyncObjects";
      if (sync_evtno_filter_input.find(":")!=std::string::npos) HelperFunctions::splitOption(sync_evtno_filter_input, sync_evtno_filter_input_fname, sync_evtno_filter_input_treename, ':');
      BaseTree tin_sync_evtno_input(sync_evtno_filter_input_fname.data(), { sync_evtno_filter_input_treename.data() }, "");
      tin_sync_evtno_input.bookBranch<BaseTree::BranchType_unknown_t>("EventNumber");
      BaseTree::BranchType sync_evtno_input_type = BaseTree::BranchType_unknown_t;
      if (!tin_sync_evtno_input.branchExists("EventNumber", &sync_evtno_input_type)){
        IVYerr << "Failed to book 'EventNumber' in the input sync. tree." << endl;
        assert(0);
      }
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) \
      type* ptr_sync_evtno_input_##name = nullptr; \
      if (sync_evtno_input_type==BaseTree::BranchType_##name##_t){ \
        IVYout << "Found a sync 'EventNumber' branch of type name '" << #name << "'." << endl; \
        tin_sync_evtno_input.getValRef("EventNumber", ptr_sync_evtno_input_##name); \
      }
      FUNDAMENTAL_DATA_INPUT_DIRECTIVES;
#undef SIMPLE_DATA_INPUT_DIRECTIVE
      for (int ev=0; ev<tin_sync_evtno_input.getNEvents(); ev++){
        tin_sync_evtno_input.getEvent(ev);
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) if (sync_evtno_input_type==BaseTree::BranchType_##name##_t) eventnumber_filter.push_back(*ptr_sync_evtno_input_##name);
        FUNDAMENTAL_DATA_INPUT_DIRECTIVES;
#undef SIMPLE_DATA_INPUT_DIRECTIVE
      }
    }
    else{
      ifstream fin_evtnos(sync_evtno_filter_input.data());
      while (!fin_evtnos.eof()){
        std::string strline;
        std::getline(fin_evtnos, strline);
        HelperFunctions::lrstrip(strline, " \n\t{}");
        if (strline.empty()) continue;
        std::vector<std::string> tmpevtnos;
        HelperFunctions::splitOptionRecursive(strline, tmpevtnos, ',', false);
        for (auto const& evtno:tmpevtnos) eventnumber_filter.push_back(std::stoi(evtno));
      }
    }
  }
  bool const has_eventnumber_filter = !eventnumber_filter.empty();
  if (has_eventnumber_filter) IVYout << "Applying filtering over " << eventnumber_filter.size() << " event numbers..." << endl;

  // Shorthand option for the Run 2 UL analysis proposal
  bool use_shorthand_Run2_UL_proposal_config = false;
  extra_arguments.getNamedVal("shorthand_Run2_UL_proposal_config", use_shorthand_Run2_UL_proposal_config);

  // Options to set alternative muon and electron IDs, or b-tagging WP
  std::string muon_id_name;
  std::string electron_id_name;
  std::string btag_WPname;
  double minpt_jets = 40;
  double minpt_bjets = 25;
  double minpt_l3 = 20;
  if (use_shorthand_Run2_UL_proposal_config){
    muon_id_name = electron_id_name = "TopMVA_Run2";
    btag_WPname = "loose";
    minpt_jets = minpt_bjets = 25;
    minpt_l3 = 10;
  }
  else{
    extra_arguments.getNamedVal("muon_id", muon_id_name);
    extra_arguments.getNamedVal("electron_id", electron_id_name);
    extra_arguments.getNamedVal("btag", btag_WPname);
    extra_arguments.getNamedVal("minpt_jets", minpt_jets);
    extra_arguments.getNamedVal("minpt_bjets", minpt_bjets);
    extra_arguments.getNamedVal("minpt_l3", minpt_l3);
  }

  double minpt_miss = 50;
  extra_arguments.getNamedVal("minpt_miss", minpt_miss);
  double minHT_jets = 300;
  extra_arguments.getNamedVal("minHT_jets", minHT_jets);

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

  double const minpt_anyjet = std::min(minpt_jets, minpt_bjets);
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
  if (theDataYear==2016){
    requiredTriggers_Dilepton = std::vector<TriggerHelpers::TriggerType>{
      TriggerHelpers::kDoubleMu_PFHT,
      TriggerHelpers::kDoubleEle_PFHT,
      TriggerHelpers::kMuEle_PFHT
    };
    // Related to triggers is how we apply loose and fakeable IDs in electrons.
    // This setting should vary for 2016 when analyzing fake rates instead of the signal or SR-like control regions.
    // If trigger choices change, this setting may not be relevant either.
    if (ElectronSelectionHelpers::selection_type == ElectronSelectionHelpers::kCutbased_Run2) ElectronSelectionHelpers::setApplyMVALooseFakeableNoIsoWPs(true);
  }
  std::vector<std::string> const hltnames_Dilepton = TriggerHelpers::getHLTMenus(requiredTriggers_Dilepton);
  auto triggerPropsCheckList_Dilepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton);

  // Declare handlers
  GenInfoHandler genInfoHandler;
  SimEventHandler simEventHandler;
  EventFilterHandler eventFilter(requiredTriggers_Dilepton);
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;

  // SF handlers
  BtagScaleFactorHandler btagSFHandler;

  // These are called handlers, but they are more like helpers.
  DileptonHandler dileptonHandler;
  ParticleDisambiguator particleDisambiguator;

  // Some advanced event filters
  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(has_multiple_dsets);
  eventFilter.setCheckHLTPathRunRanges(true);

  // Allow extraction of optional lepton information
  muonHandler.setExtractOptionalInfo(writeSyncObjects);
  electronHandler.setExtractOptionalInfo(writeSyncObjects);

  curdir->cd();

  // We conclude the setup of event processing specifications and move on to I/O configuration next.

  // Open the output ROOT file
  TFile* foutput = TFile::Open(stroutput, "recreate"); foutput->cd();
  BaseTree* tout = nullptr;
  SimpleEntry rcd_output;
  if (produce_trees){
    tout = new BaseTree("SkimTree");
    tout->setAutoSave(0);
  }

  std::vector<TString> const gencats{ "AllLeptonsGenMatched_SameCharge", "hasGenPromptLepton_ChargeFlip", "hasGenMatchedPromptPhotonProduct", "hasNonPromptLepton" };
  std::vector<TH2D*> hCats;
  if (!produce_trees){
    hCats.reserve(gencats.size());
    for (auto const& gencat:gencats){
      TH2D* hCat = new TH2D(Form("hCat_%s", gencat.Data()), "", 15, 0, 15, 2, 0, 2); hCat->Sumw2();
      for (int ix=1; ix<=15; ix++){
        if (ix<15) hCat->GetXaxis()->SetBinLabel(ix, Form("SR%i", ix));
        else hCat->GetXaxis()->SetBinLabel(ix, "CRW");
      }
      hCat->GetYaxis()->SetBinLabel(1, "SR");
      hCat->GetYaxis()->SetBinLabel(2, "CRZ");
      hCats.push_back(hCat);
    }
  }
  curdir->cd();

  // Acquire input tree/chains
  TString strinputdpdir = theDataPeriod;
  if (SampleHelpers::testDataPeriodIsLikeData(theDataPeriod)){
    if (theDataYear==2016){
      if (SampleHelpers::isAPV2016Affected(theDataPeriod)) strinputdpdir = Form("%i_APV", theDataYear);
      else strinputdpdir = Form("%i_NonAPV", theDataYear);
    }
    else strinputdpdir = Form("%i", theDataYear);
  }

  signed char is_sim_data_flag = -1; // =0 for sim, =1 for data
  int nevents_total = 0;
  unsigned int nevents_total_traversed = 0;
  std::vector<BaseTree*> tinlist; tinlist.reserve(dset_proc_pairs.size());
  std::unordered_map<BaseTree*, double> tin_normScale_map;
  std::unordered_map<BaseTree*, double> tin_normScale_noPU_map;
  for (auto const& dset_proc_pair:dset_proc_pairs){
    TString strinput = SampleHelpers::getInputDirectory() + "/" + strinputdpdir + "/" + dset_proc_pair.second.data();
    TString cinput = (input_files=="" ? strinput + "/*.root" : strinput + "/" + input_files.data());
    IVYout << "Accessing input files " << cinput << "..." << endl;
    TString const sid = SampleHelpers::getSampleIdentifier(dset_proc_pair.first);
    bool const isData = SampleHelpers::checkSampleIsData(sid);
    BaseTree* tin = new BaseTree(cinput, "Events", "", (isData ? "" : "Counters"));
    if (!tin->isValid()){
      IVYout << "An error occured while acquiring the input from " << cinput << ". Aborting..." << endl;
      assert(0);
    }
    tin->sampleIdentifier = sid;
    if (!isData){
      if (xsec<0.){
        IVYerr << "xsec = " << xsec << " is not valid." << endl;
        assert(0);
      }
      if (has_multiple_dsets){
        IVYerr << "Should not process multiple data sets in sim. mode. Aborting..." << endl;
        assert(0);
      }
    }
    if (is_sim_data_flag==-1) is_sim_data_flag = (isData ? 1 : 0);
    else if (is_sim_data_flag != (isData ? 1 : 0)){
      IVYerr << "Should not process data and simulation at the same time." << endl;
      assert(0);
    }

    double sum_wgts = (isData ? 1 : 0);
    double sum_wgts_noPU = (isData ? 1 : 0);
    if (!isData){
      int ix = 1;
      switch (theGlobalSyst){
      case SystematicsHelpers::ePUDn:
        ix = 2;
        break;
      case SystematicsHelpers::ePUUp:
        ix = 3;
        break;
      default:
        break;
      }
      TH2D* hCounters = (TH2D*) tin->getCountersHistogram();
      sum_wgts = hCounters->GetBinContent(ix, 1);
      sum_wgts_noPU = hCounters->GetBinContent(0, 0);
    }
    if (sum_wgts==0.){
      IVYerr << "Sum of pre-recorded weights cannot be zero." << endl;
      assert(0);
    }

    // Add the tree to the list of trees to process
    tinlist.push_back(tin);

    // Calculate the overall normalization scale on the events.
    // Includes xsec (in fb), lumi (in fb-1), and 1/sum of weights in all of the MC.
    // Data normalizaion factor is always 1.
    double norm_scale = (isData ? 1. : xsec * xsecScale * lumi)/sum_wgts;
    tin_normScale_map[tin] = norm_scale;
    double norm_scale_noPU = (isData ? 1. : xsec * xsecScale * lumi)/sum_wgts_noPU;
    tin_normScale_noPU_map[tin] = norm_scale_noPU;
    IVYout << "Acquired a sum of weights of " << sum_wgts << ". Overall normalization will be " << norm_scale << "." << endl;
    IVYout << "\t- Without PU reweighting, the sum of weights would have been " << sum_wgts_noPU << " instead. Overall normalization would have become " << norm_scale_noPU << "." << endl;

    nevents_total += tin->getNEvents();

    curdir->cd();

    // Book the necessary branches
    genInfoHandler.bookBranches(tin);
    simEventHandler.bookBranches(tin);
    eventFilter.bookBranches(tin);
    muonHandler.bookBranches(tin);
    electronHandler.bookBranches(tin);
    jetHandler.bookBranches(tin);
    if (useIsotrackVeto) isotrackHandler.bookBranches(tin);

    // Book a few additional branches
    tin->bookBranch<EventNumber_t>("event", 0);
    if ((runSyncExercise || produce_trees) && !isData){
      tin->bookBranch<float>("GenMET_pt", 0);
      tin->bookBranch<float>("GenMET_phi", 0);
    }
    if (isData){
      tin->bookBranch<RunNumber_t>("run", 0);
      tin->bookBranch<LuminosityBlock_t>("luminosityBlock", 0);
    }
  }

  curdir->cd();

  // Prepare to loop!
  // Keep track of sums of weights
  SelectionTracker seltracker;

  // Keep track of the traversed events
  bool firstOutputEvent = true;
  bool firstSyncObjectsOutputEvent = true;
  int eventIndex_begin = -1;
  int eventIndex_end = -1;
  int eventIndex_tracker = 0;
  splitInputEventsIntoChunks((is_sim_data_flag==1), nevents_total, ichunk, nchunks, eventIndex_begin, eventIndex_end);

  for (auto const& tin:tinlist){
    if (SampleHelpers::doSignalInterrupt==1) break;

    auto const& norm_scale = tin_normScale_map.find(tin)->second;
    auto const& norm_scale_noPU = tin_normScale_noPU_map.find(tin)->second;
    bool const isData = (is_sim_data_flag==1);

    // Wrap the ivies around the input tree:
    // Booking is basically SetBranchStatus+SetBranchAddress. You can book for as many trees as you would like.
    // In some cases, bookBranches also informs the ivy dynamically that it is supposed to consume certain entries.
    // For entries common to all years or any data or MC, the consumption information is handled in the ivy constructor already.
    // None of these mean the ivy establishes its access to the input tree yet.
    // Wrapping a tree informs the ivy that it is supposed to consume the booked entries from that particular tree.
    // Without wrapping, you are not really accessing the entries from the input tree to construct the physics objects;
    // all you would get are 0 electrons, 0 jets, everything failing event filters etc.

    genInfoHandler.wrapTree(tin);
    simEventHandler.wrapTree(tin);
    eventFilter.wrapTree(tin);
    muonHandler.wrapTree(tin);
    electronHandler.wrapTree(tin);
    jetHandler.wrapTree(tin);
    if (useIsotrackVeto) isotrackHandler.wrapTree(tin);

    RunNumber_t* ptr_RunNumber = nullptr;
    LuminosityBlock_t* ptr_LuminosityBlock = nullptr;
    EventNumber_t* ptr_EventNumber = nullptr;
    float* ptr_genmet_pt = nullptr;
    float* ptr_genmet_phi = nullptr;
    tin->getValRef("event", ptr_EventNumber);
    if ((runSyncExercise || produce_trees) && !isData){
      tin->getValRef("GenMET_pt", ptr_genmet_pt);
      tin->getValRef("GenMET_phi", ptr_genmet_phi);
    }
    if (isData){
      tin->getValRef("run", ptr_RunNumber);
      tin->getValRef("luminosityBlock", ptr_LuminosityBlock);
    }

    unsigned int n_traversed = 0;
    unsigned int n_recorded = 0;
    int nEntries = tin->getNEvents();
    IVYout << "Looping over " << nEntries << " events from " << tin->sampleIdentifier << "..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      bool doAccumulate = true;
      if (isData){
        if (eventIndex_begin>0 || eventIndex_end>0) doAccumulate = (
          tin->updateBranch(ev, "run", false)
          &&
          (eventIndex_begin<0 || (*ptr_RunNumber)>=static_cast<RunNumber_t>(eventIndex_begin))
          &&
          (eventIndex_end<0 || (*ptr_RunNumber)<=static_cast<RunNumber_t>(eventIndex_end))
          );
      }
      else doAccumulate = (
        (eventIndex_begin<0 || eventIndex_tracker>=static_cast<int>(eventIndex_begin))
        &&
        (eventIndex_end<0 || eventIndex_tracker<static_cast<int>(eventIndex_end))
        );
      if (has_eventnumber_filter){
        if (eventnumber_filter.empty()) break;
        doAccumulate &= (
          tin->updateBranch(ev, "event", false)
          &&
          HelperFunctions::removeListVariable(eventnumber_filter, *ptr_EventNumber)
          );
      }

      eventIndex_tracker++;
      if (!doAccumulate) continue;

      tin->getEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);
      n_traversed++;

      genInfoHandler.constructGenInfo();
      auto const& genInfo = genInfoHandler.getGenInfo();
      auto const& genparticles = genInfoHandler.getGenParticles();

      std::vector<GenParticleObject*> genmatch_promptleptons;
      std::vector<GenParticleObject*> genmatch_photonproducts;
      for (auto const& part:genparticles){
        if (part->status()!=1) continue; // Only final states
        bool const isLepton = std::abs(part->pdgId())==11 || std::abs(part->pdgId())==13;
        bool const isPhoton = part->pdgId()==22;
        if (!(isLepton || isPhoton)) continue;

        GenParticleObject* photonMother = nullptr;
        if (isLepton) photonMother = dynamic_cast<GenParticleObject*>(ParticleObject::getFirstMotherInChain_matchPDGId(part, 22));

        if (isLepton && part->extras.isPromptFinalState) genmatch_promptleptons.push_back(part);
        else if ((isPhoton && part->extras.isPromptFinalState) || (photonMother && photonMother->extras.isPrompt)) genmatch_photonproducts.push_back(part);
      }
      std::vector<GenParticleObject*> genmatchparticles;
      HelperFunctions::appendVector(genmatchparticles, genmatch_promptleptons);
      HelperFunctions::appendVector(genmatchparticles, genmatch_photonproducts);

      simEventHandler.constructSimEvent();

      double genwgt = 1;
      double puwgt = 1;
      double l1prefiringwgt = 1;
      if (!isData){
        genwgt = genInfo->getGenWeight(theGlobalSyst);
        puwgt = simEventHandler.getPileUpWeight(theGlobalSyst);
        // When ignore_PU_rewgt=true, make sure puwgt is scaled correctly.
        if (ignore_PU_rewgt) puwgt *= norm_scale / norm_scale_noPU;
        // L1 prefiring weight for 2016 and 2017
        l1prefiringwgt = simEventHandler.getL1PrefiringWeight(theGlobalSyst);
      }
      double wgt = genwgt * (!ignore_L1Prefiring_rewgt ? l1prefiringwgt : 1.) * (!ignore_PU_rewgt ? puwgt*norm_scale : norm_scale_noPU);

      muonHandler.constructMuons();
      electronHandler.constructElectrons();
      jetHandler.constructJetMET(&genInfoHandler, &simEventHandler, theGlobalSyst);

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
        false;
        //HelperFunctions::checkListVariable(std::vector<int>{ 16603 }, ev);
        //HelperFunctions::checkListVariable(std::vector<int>{ 1902, 5855, 6073, 7046, 11794, 16603 }, ev);
        //HelperFunctions::checkListVariable(std::vector<int>{ 1233, 1475, 1546, 1696, 2011, 2103, 2801, 2922, 3378, 3407, 3575, 3645, 5021, 5127, 6994, 7000, 7046, 7341, 7351, 8050, 9931, 10063, 10390, 10423, 10623, 10691, 10791, 10796, 11127, 11141, 11279, 11794, 12231, 12996, 13115, 13294, 13550, 14002, 14319, 15062, 15754, 16153, 16166, 16316, 16896, 16911, 17164 }, ev);
        //HelperFunctions::checkListVariable(std::vector<int>{663, 1469, 3087, 3281}, ev);
        //HelperFunctions::checkListVariable(std::vector<int>{204, 353, 438, 1419}, ev);
        //HelperFunctions::checkListVariable(std::vector<int>{3, 15, 30, 31, 32, 41, 153, 154, 162, 197, 215, 284, 317, 360, 572, 615, 747, 1019, 1119, 1129}, ev);

    // Muon sync. write variables
#define SYNC_MUONS_BRANCH_VECTOR_COMMANDS \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_loose) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_fakeable) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_tight) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_genmatched_prompt) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_genmatched_photonProduct) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_genmatched_prompt_sameCharge) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_genmatched_prompt_chargeFlip) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, muons, is_genmatched_prompt_photonProduct) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(int, muons, genmatch_pdgId) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(int, muons, matchedJetIdentifier) \
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
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_genmatched_prompt) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_genmatched_photonProduct) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_genmatched_prompt_sameCharge) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_genmatched_prompt_chargeFlip) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(bool, electrons, is_genmatched_prompt_photonProduct) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(int, electrons, genmatch_pdgId) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(int, electrons, matchedJetIdentifier) \
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
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(unsigned int, ak4jets, uniqueIdenntifier) \
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

      // Get leptons and match them to gen. particles
      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();

      std::vector<ParticleObject*> leptons; leptons.reserve(muons.size()+electrons.size());
      for (auto const& part:muons) leptons.push_back(part);
      for (auto const& part:electrons) leptons.push_back(part);
      std::unordered_map<ParticleObject*, GenParticleObject*> lepton_genmatchpart_map;
      ParticleObjectHelpers::matchParticles(
        ParticleObjectHelpers::kMatchBy_DeltaR, 0.2,
        leptons.begin(), leptons.end(),
        genmatchparticles.begin(), genmatchparticles.end(),
        lepton_genmatchpart_map
      );

      std::unordered_map<ParticleObject*, bool> leptons_isGenMatched_map;
      std::unordered_map<ParticleObject*, bool> leptons_isGenMatched_photonProduct_map;
      std::unordered_map<ParticleObject*, int> leptons_genMatchPDGId_map;
      for (auto const& part:leptons){
        bool is_genmatched_prompt = false;
        bool is_genmatched_photonProduct = false;
        int genmatch_pdgId = -9000; // This PDG id means not a valid particle.
        GenParticleObject* genpart = nullptr;
        auto it_genpart = lepton_genmatchpart_map.find(part);
        if (it_genpart!=lepton_genmatchpart_map.end()) genpart = it_genpart->second;
        if (genpart){
          is_genmatched_prompt = genpart->extras.isPrompt;
          is_genmatched_photonProduct = HelperFunctions::checkListVariable(genmatch_photonproducts, genpart);
          genmatch_pdgId = genpart->pdgId();
        }
        leptons_isGenMatched_map[part] = is_genmatched_prompt;
        leptons_isGenMatched_photonProduct_map[part] = is_genmatched_photonProduct;
        leptons_genMatchPDGId_map[part] = genmatch_pdgId;
      }

      if (printObjInfo) IVYout << "Lepton info for event " << ev << ":" << endl;

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

        // The following two variables are guaranteed to exist from global matching to all leptons.
        bool is_genmatched_prompt = leptons_isGenMatched_map.find(dynamic_cast<ParticleObject*>(part))->second;
        bool is_genmatched_photonProduct = leptons_isGenMatched_photonProduct_map.find(dynamic_cast<ParticleObject*>(part))->second;
        int genmatch_pdgId = leptons_genMatchPDGId_map.find(dynamic_cast<ParticleObject*>(part))->second;

        bool is_genmatched_prompt_sameCharge = is_genmatched_prompt && genmatch_pdgId==part->pdgId();
        bool is_genmatched_prompt_chargeFlip = is_genmatched_prompt && std::abs(genmatch_pdgId)==std::abs(part->pdgId()) && genmatch_pdgId!=part->pdgId();
        bool is_genmatched_prompt_photonProduct = is_genmatched_prompt && is_genmatched_photonProduct;

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
        int matchedJetIdentifier = (mother ? mother->getUniqueIdentifier() : -1);

        if (printObjInfo){
          IVYout
            << "\t- PDG id = " << part->pdgId()
            << ", pt = " << pt << ", eta = " << eta << ", phi = " << phi
            << ", array index = " << part->getUniqueIdentifier()
            << ", matched jet index = " << (mother ? static_cast<int>(mother->getUniqueIdentifier()) : -1)
            << ", loose? " << is_loose
            << ", fakeable? " << is_fakeable
            << ", tight? " << is_tight
            << endl;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) << "\t\t- " << #NAME << ": " << static_cast<double>(part->extras.NAME) << "\n"
          IVYout
            MUON_EXTRA_VARIABLES
            << "\t\t- b score = " << bscore << "\n"
            << "\t\t- ptrel = " << ptrel_final << "\n"
            << "\t\t- ptratio = " << ptratio_final
            << endl;
#undef MUON_VARIABLE
          if (has_extMVAscore) IVYout << "\t\t- External MVA score = " << extMVAscore << endl;
        }

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

        // The following two variables are guaranteed to exist from global matching to all leptons.
        bool is_genmatched_prompt = leptons_isGenMatched_map.find(dynamic_cast<ParticleObject*>(part))->second;
        bool is_genmatched_photonProduct = leptons_isGenMatched_photonProduct_map.find(dynamic_cast<ParticleObject*>(part))->second;
        int genmatch_pdgId = leptons_genMatchPDGId_map.find(dynamic_cast<ParticleObject*>(part))->second;

        bool is_genmatched_prompt_sameCharge = is_genmatched_prompt && genmatch_pdgId==part->pdgId();
        bool is_genmatched_prompt_chargeFlip = is_genmatched_prompt && std::abs(genmatch_pdgId)==std::abs(part->pdgId()) && genmatch_pdgId!=part->pdgId();
        bool is_genmatched_prompt_photonProduct = is_genmatched_prompt && is_genmatched_photonProduct;

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
        int matchedJetIdentifier = (mother ? mother->getUniqueIdentifier() : -1);

        if (printObjInfo){
          IVYout
            << "\t- PDG id = " << part->pdgId()
            << ", pt = " << pt << ", eta = " << eta << ", phi = " << phi
            << ", array index = " << part->getUniqueIdentifier()
            << ", matched jet index = " << (mother ? static_cast<int>(mother->getUniqueIdentifier()) : -1)
            << ", loose? " << is_loose
            << ", fakeable? " << is_fakeable
            << ", tight? " << is_tight
            << endl;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) << "\t\t- " << #NAME << ": " << static_cast<double>(part->extras.NAME) << "\n"
          IVYout
            ELECTRON_EXTRA_VARIABLES
            << "\t\t- Fall 17 v2 no-iso. MVA transformed = " << mvaFall17V2noIso_raw << "\n"
            << "\t\t- b score = " << bscore << "\n"
            << "\t\t- ptrel = " << ptrel_final << "\n"
            << "\t\t- ptratio = " << ptratio_final
            << endl;
#undef ELECTRON_VARIABLE
          if (has_extMVAscore) IVYout << "\t\t- External MVA score = " << extMVAscore << endl;
        }

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

      bool hasGenPromptLepton_ChargeFlip = false;
      bool hasGenMatchedPromptPhotonProduct = false;
      bool hasNonPromptLepton = false;
      for (auto const& part:leptons_tight){
        bool is_genmatched_prompt = leptons_isGenMatched_map.find(part)->second;
        bool is_genmatched_photonProduct = leptons_isGenMatched_photonProduct_map.find(part)->second;
        int genmatch_pdgId = leptons_genMatchPDGId_map.find(part)->second;

        if (is_genmatched_prompt && std::abs(genmatch_pdgId)==std::abs(part->pdgId())) hasGenPromptLepton_ChargeFlip |= (genmatch_pdgId!=part->pdgId());
        else if (is_genmatched_prompt && is_genmatched_photonProduct) hasGenMatchedPromptPhotonProduct |= true;
        else hasNonPromptLepton |= true;
      }

      bool hasVetoIsotrack = false;
      if (useIsotrackVeto){
        isotrackHandler.constructIsotracks(&muons, &electrons);
        for (auto const& isotrack:isotrackHandler.getProducts()){
          if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
            hasVetoIsotrack = true;
            break;
          }
        }
      }

      if (printObjInfo) IVYout << "Jet info for event " << ev << ":" << endl;
      double event_wgt_SFs_btagging = 1;
      float HT_ak4jets=0;
      auto const& ak4jets = jetHandler.getAK4Jets();
      std::vector<AK4JetObject*> ak4jets_tight_selected;
      std::vector<AK4JetObject*> ak4jets_tight_selected_btagged;
      for (auto const& jet:ak4jets){
        float pt = jet->pt();
        float eta = jet->eta();
        float phi = jet->phi();
        float mass = jet->mass();

        bool is_tight = ParticleSelectionHelpers::isTightJet(jet);
        bool is_btagged = jet->testSelectionBit(bit_preselection_btag);
        constexpr bool is_clean = true;
        unsigned int uniqueIdenntifier = jet->getUniqueIdentifier();

        if (!isData){
          float theSF_btag = 1;
          float theEff_btag = 1;
          btagSFHandler.getSFAndEff(theGlobalSyst, jet, theSF_btag, &theEff_btag); theSF_btag = std::max(theSF_btag, 1e-5f); event_wgt_SFs_btagging *= theSF_btag;
          if (theSF_btag<=1e-5f){
            IVYout
              << "Jet has b-tagging SF<=1e-5:"
              << "\n\t- pt = " << jet->pt()
              << "\n\t- eta = " << jet->eta()
              << "\n\t- b kin = " << jet->testSelectionBit(AK4JetSelectionHelpers::kKinOnly_BTag)
              << "\n\t- Loose = " << jet->testSelectionBit(AK4JetSelectionHelpers::kBTagged_Loose)
              << "\n\t- Medium = " << jet->testSelectionBit(AK4JetSelectionHelpers::kBTagged_Medium)
              << "\n\t- Tight = " << jet->testSelectionBit(AK4JetSelectionHelpers::kBTagged_Tight)
              << "\n\t- Eff =  " << theEff_btag
              << "\n\t- SF = " << theSF_btag
              << endl;
          }
        }

        if (printObjInfo) IVYout
          << "\t- pt = " << pt << ", eta = " << eta << ", phi = " << phi
          << ", array index = " << jet->getUniqueIdentifier()
          << " | btagged? " << is_btagged
          << " | tight? " << is_tight
          << " | clean? " << is_clean
          << " ("
          << "jetId = " << jet->extras.jetId
          << ", electron indices = " << std::vector<int>{ jet->extras.electronIdx1, jet->extras.electronIdx2 }
        << ", muon indices = " << std::vector<int>{ jet->extras.muonIdx1, jet->extras.muonIdx2 }
        << ")"
          << endl;

        if (writeSyncObjects){
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.push_back(NAME);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) ak4jets_##NAME.push_back(jet->extras.NAME);
          SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND
        }

        if (is_tight && pt>=minpt_anyjet && std::abs(eta)<absEtaThr_ak4jets){
          if (is_btagged && pt>=minpt_bjets) ak4jets_tight_selected_btagged.push_back(jet);
          if (pt>=minpt_jets){
            ak4jets_tight_selected.push_back(jet);
            HT_ak4jets += pt;
          }
        }
      }
      if (writeSyncObjects || printObjInfo){
        for (auto const& jet:jetHandler.getMaskedAK4Jets()){
          float pt = jet->pt();
          float eta = jet->eta();
          float phi = jet->phi();
          float mass = jet->mass();

          bool is_tight = ParticleSelectionHelpers::isTightJet(jet);
          bool is_btagged = jet->testSelectionBit(bit_preselection_btag);
          constexpr bool is_clean = false;
          unsigned int uniqueIdenntifier = jet->getUniqueIdentifier();

          if (printObjInfo) IVYout
            << "\t- pt = " << pt << ", eta = " << eta << ", phi = " << phi
            << ", array index = " << jet->getUniqueIdentifier()
            << " | btagged? " << is_btagged
            << " | tight? " << is_tight
            << " | clean? " << is_clean
            << " ("
            << "jetId = " << jet->extras.jetId
            << ", electron indices = " << std::vector<int>{ jet->extras.electronIdx1, jet->extras.electronIdx2 }
          << ", muon indices = " << std::vector<int>{ jet->extras.muonIdx1, jet->extras.muonIdx2 }
          << ")"
            << endl;

          if (writeSyncObjects){
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.push_back(NAME);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) ak4jets_##NAME.push_back(jet->extras.NAME);
            SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND
          }
        }
      }
      unsigned int const nak4jets_tight_selected = ak4jets_tight_selected.size();
      unsigned int const nak4jets_tight_selected_btagged = ak4jets_tight_selected_btagged.size();

      auto const& eventmet = jetHandler.getPFMET();
      float const pTmiss = eventmet->pt();
      float const phimiss = eventmet->phi();

      // Write object sync. info.
      if (writeSyncObjects){
        rcd_sync_objects.setNamedVal("EventNumber", *ptr_EventNumber);
        if (!isData){
          rcd_sync_objects.setNamedVal("GenMET_pt", *ptr_genmet_pt);
          rcd_sync_objects.setNamedVal("GenMET_phi", *ptr_genmet_phi);
        }
        else{
          rcd_sync_objects.setNamedVal("RunNumber", *ptr_RunNumber);
          rcd_sync_objects.setNamedVal("LuminosityBlock", *ptr_LuminosityBlock);
        }
        rcd_sync_objects.setNamedVal("PFMET_pt_final", eventmet->pt());
        rcd_sync_objects.setNamedVal("PFMET_phi_final", eventmet->phi());

        rcd_sync_objects.setNamedVal("hasGenPromptLepton_ChargeFlip", hasGenPromptLepton_ChargeFlip);
        rcd_sync_objects.setNamedVal("hasGenMatchedPromptPhotonProduct", hasGenMatchedPromptPhotonProduct);
        rcd_sync_objects.setNamedVal("hasNonPromptLepton", hasNonPromptLepton);

#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) rcd_sync_objects.setNamedVal(Form("%s_%s", #COLLNAME, #NAME), COLLNAME##_##NAME);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, muons, NAME)
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, electrons, NAME)
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, ak4jets, NAME)
        SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef ELECTRON_VARIABLE
#undef MUON_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND

        if (firstSyncObjectsOutputEvent){
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
        firstSyncObjectsOutputEvent = false;
      }

      // BEGIN PRESELECTION
      seltracker.accumulate("Full sample", wgt, printObjInfo);

      bool const pass_Nj_geq_2 = nak4jets_tight_selected>=2; rcd_output.setNamedVal("pass_Nj_geq_2", pass_Nj_geq_2);
      bool const pass_Nb_geq_2 = nak4jets_tight_selected_btagged>=2; rcd_output.setNamedVal("pass_Nb_geq_2", pass_Nb_geq_2);
      if (!produce_trees && applyPreselection && (!pass_Nj_geq_2 || !pass_Nb_geq_2)) continue;
      seltracker.accumulate("Pass Nj>=2 and Nb>=2", wgt, printObjInfo);

      bool const pass_pTmiss = pTmiss>=minpt_miss; rcd_output.setNamedVal("pass_pTmiss", pass_pTmiss);
      if (!produce_trees && applyPreselection && !pass_pTmiss) continue;
      seltracker.accumulate("Pass pTmiss", wgt, printObjInfo);

      bool const pass_HTjets = HT_ak4jets>=minHT_jets; rcd_output.setNamedVal("pass_HTjets", pass_HTjets);
      if (!produce_trees && applyPreselection && !pass_HTjets) continue;
      seltracker.accumulate("Pass HT", wgt, printObjInfo);

      bool const pass_Nleptons = (nleptons_tight>=2 && nleptons_tight<5); rcd_output.setNamedVal("pass_Nleptons", pass_Nleptons);
      if (!produce_trees && nleptons_tight<2) continue; // Skims are required to apply this selection, so no additional test on applyPreselection.
      if (!produce_trees && applyPreselection && !pass_Nleptons) continue;
      seltracker.accumulate("Has >=2 and <=4 leptons, >=2 of which are tight", wgt, printObjInfo);

      bool const pass_pTl1 = (nleptons_tight<1 || leptons_tight.front()->pt()>=25.); rcd_output.setNamedVal("pass_pTl1", pass_pTl1);
      bool const pass_pTl2 = (nleptons_tight<2 || leptons_tight.at(1)->pt()>=20.); rcd_output.setNamedVal("pass_pTl2", pass_pTl2);
      if (!produce_trees && applyPreselection && (!pass_pTl1 || !pass_pTl2)) continue;
      seltracker.accumulate("Pass pTl1 and pTl2", wgt, printObjInfo);

      bool const pass_pTl3 = (nleptons_tight<3 || leptons_tight.at(2)->pt()>=minpt_l3); rcd_output.setNamedVal("pass_pTl3", pass_pTl3);
      if (!produce_trees && applyPreselection && !pass_pTl3) continue;
      seltracker.accumulate("Pass pTl3 if >=3 tight leptons", wgt, printObjInfo);

      // Construct all possible dilepton pairs
      int nQ = 0;
      for (auto const& part:leptons_tight) nQ += (part->pdgId()>0 ? -1 : 1);
      bool const pass_trileptonSameCharge = (std::abs(nQ)<(6-static_cast<int>(nleptons_tight))); rcd_output.setNamedVal("pass_trileptonSameCharge", pass_trileptonSameCharge); // This req. necessarily vetoes Nleps>=5 because the actual number of same-sign leptons will always be >=3.
      if (!produce_trees && applyPreselection && !pass_trileptonSameCharge) continue;
      seltracker.accumulate("Pass 3-lepton same charge veto", wgt, printObjInfo);

      dileptonHandler.constructDileptons((!only_tight_dileptons ? &muons_selected : &muons_tight), (!only_tight_dileptons ? &electrons_selected : &electrons_tight));
      auto const& dileptons = dileptonHandler.getProducts();
      unsigned int ndileptons_SS = 0;
      unsigned int ndileptons_OS = 0;
      for (auto const& dilepton:dileptons){
        if (dilepton->isOS()) ndileptons_OS++;
        else ndileptons_SS++;
      }
      // If uncommented, these numbers could later tell you the average number of OS and SS dileptons before the rest of the selections
      // once you divide the numbers by the previous gen. weight sum.
      //seltracker.accumulate("nOS", wgt*static_cast<double>(ndileptons_OS), printObjInfo);
      //seltracker.accumulate("nSS", wgt*static_cast<double>(ndileptons_SS), printObjInfo);

      if (printObjInfo) IVYout << "Dilepton info for event " << ev << ":" << endl;
      bool fail_vetos = false;
      DileptonObject* dilepton_SS_tight = nullptr;
      DileptonObject* dilepton_OS_DYCand_tight = nullptr;
      for (auto const& dilepton:dileptons){
        bool isSS = !dilepton->isOS();
        bool isTight = dilepton->nTightDaughters()==2;
        bool isSF = dilepton->isSF();
        bool is_LowMass = dilepton->m()<12.;
        bool is_ZClose = std::abs(dilepton->m()-91.2)<15.;
        bool is_DYClose = std::abs(dilepton->m()-91.2)<15. || is_LowMass;
        if (printObjInfo){
          IVYout
            << "\t- Dilepton: pt = " << dilepton->pt() << ", mass = " << dilepton->mass()
            << ", daughter ids = { " << dilepton->getDaughter(0)->pdgId() << " [" << dilepton->getDaughter(0)->getUniqueIdentifier()
            << "], " << dilepton->getDaughter(1)->pdgId() << " [" << dilepton->getDaughter(1)->getUniqueIdentifier()  << "] }"
            << '\n';
          IVYout << "\t\t- isSS = " << isSS << '\n';
          IVYout << "\t\t- isTight = " << isTight << '\n';
          IVYout << "\t\t- isSF = " << isSF << '\n';
          IVYout << "\t\t- is_LowMass = " << is_LowMass << '\n';
          IVYout << "\t\t- is_ZClose = " << is_ZClose << '\n';
          IVYout << "\t\t- is_DYClose = " << is_DYClose << endl;
        }
        if (isSS && isSF && is_LowMass && std::abs(dilepton->getDaughter(0)->pdgId())==11){
          fail_vetos = true;
          if (applyPreselection) break;
        }
        if (isSS && isTight && !dilepton_SS_tight) dilepton_SS_tight = dilepton;
        if (!isSS && isSF && is_DYClose){
          if (isTight && is_ZClose && !dilepton_OS_DYCand_tight) dilepton_OS_DYCand_tight = dilepton;
          else{
            fail_vetos = true;
            if (applyPreselection) break;
          }
        }
      }

      bool const pass_dileptonVetos = !fail_vetos; rcd_output.setNamedVal("pass_dileptonVetos", pass_dileptonVetos);
      if (!produce_trees && applyPreselection && fail_vetos) continue;
      seltracker.accumulate("Pass dilepton vetos", wgt, printObjInfo);

      bool const has_dilepton_SS_tight = (dilepton_SS_tight!=nullptr); rcd_output.setNamedVal("has_dilepton_SS_tight", has_dilepton_SS_tight);
      bool const has_dilepton_OS_DYCand_tight = (dilepton_OS_DYCand_tight!=nullptr); rcd_output.setNamedVal("has_dilepton_OS_DYCand_tight", has_dilepton_OS_DYCand_tight);
      if (!produce_trees && applyPreselection && !has_dilepton_SS_tight) continue;
      seltracker.accumulate("Has at least one tight SS dilepton", wgt, printObjInfo);

      // Put event filters to the last because data has unique event tracking enabled.
      eventFilter.constructFilters(&simEventHandler);
      constexpr bool pass_HEMveto = true;
      //bool const pass_HEMveto = eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets);
      //bool const pass_HEMveto = eventFilter.test2018HEMFilter(&simEventHandler, &electrons, nullptr, &ak4jets);
      rcd_output.setNamedVal("pass_HEMveto", pass_HEMveto);
      if (!produce_trees && applyPreselection && !pass_HEMveto) continue; // Test for 2018 partial HEM failure
      seltracker.accumulate("Pass HEM veto", wgt, printObjInfo);
      bool const pass_METFilters = eventFilter.passMETFilters(); rcd_output.setNamedVal("pass_METFilters", pass_METFilters);
      if (!produce_trees && applyPreselection && !pass_METFilters) continue; // Test for MET filters
      seltracker.accumulate("Pass MET filters", wgt, printObjInfo);
      if (!eventFilter.isUniqueDataEvent()) continue; // Test if the data event is unique (i.e., dorky). Does not do anything in the MC.
      seltracker.accumulate("Pass unique event check", wgt, printObjInfo);

      // Triggers
      float event_wgt_triggers_dilepton = eventFilter.getTriggerWeight(hltnames_Dilepton);
      bool const pass_triggers_dilepton = (event_wgt_triggers_dilepton!=0.f); rcd_output.setNamedVal("pass_triggers_dilepton", pass_triggers_dilepton);
      if (!produce_trees && applyPreselection && !pass_triggers_dilepton) continue; // Test if any triggers passed at all
      seltracker.accumulate("Pass any trigger", wgt, printObjInfo);
      float event_wgt_triggers_dilepton_matched = eventFilter.getTriggerWeight(
        triggerPropsCheckList_Dilepton,
        &muons, &electrons, nullptr, &ak4jets, nullptr, nullptr
      );
      bool const pass_triggers_dilepton_matched = (event_wgt_triggers_dilepton_matched!=0.f); rcd_output.setNamedVal("pass_triggers_dilepton_matched", pass_triggers_dilepton_matched);
      seltracker.accumulate("Pass triggers after matching", wgt*static_cast<double>(pass_triggers_dilepton_matched), printObjInfo);

      // Optional check for an isotrack veto
      if (useIsotrackVeto){
        rcd_output.setNamedVal("pass_isotrack_vetos", !hasVetoIsotrack);
        if (!produce_trees && applyPreselection && hasVetoIsotrack) continue;
        seltracker.accumulate("Pass isolated track veto", wgt, printObjInfo);
      }

      /*************************************************/
      /* NO MORE CALLS TO SELECTION BEYOND THIS POINT! */
      /*************************************************/
      constexpr int idx_CRW = 15;
      int iCRZ = (dilepton_OS_DYCand_tight ? 1 : 0);
      int icat = 0;
      if (nak4jets_tight_selected_btagged>=2){
        if (nleptons_tight==2){
          switch (nak4jets_tight_selected_btagged){
          case 2:
            if (nak4jets_tight_selected<6) icat = idx_CRW;
            else icat = 1 + std::min(nak4jets_tight_selected, static_cast<unsigned int>(8))-6;
            break;
          case 3:
            if (nak4jets_tight_selected>=5) icat = 4 + std::min(nak4jets_tight_selected, static_cast<unsigned int>(8))-5;
            break;
          default: // Nb>=4
            if (nak4jets_tight_selected>=5) icat = 8;
            break;
          }
        }
        else{
          switch (nak4jets_tight_selected_btagged){
          case 2:
            if (nak4jets_tight_selected>=5) icat = 9 + std::min(nak4jets_tight_selected, static_cast<unsigned int>(7))-5;
            break;
          default: // Nb>=3
            if (nak4jets_tight_selected>=4) icat = 12 + std::min(nak4jets_tight_selected, static_cast<unsigned int>(6))-4;
            break;
          }
        }
      }

      if (tout){
        rcd_output.setNamedVal<float>("event_wgt", wgt);
        if (!isData){
          // If PU or L1 prefiring weights are ignored from the event reweighting scheme,
          // we should still be able to access them using event weight adjustments later.
          if (ignore_PU_rewgt) rcd_output.setNamedVal<float>("event_wgt_adjustment_PU", puwgt);
          if (ignore_L1Prefiring_rewgt) rcd_output.setNamedVal<float>("event_wgt_adjustment_L1Prefiring", l1prefiringwgt);
        }
        rcd_output.setNamedVal<float>("event_wgt_triggers_dilepton", event_wgt_triggers_dilepton);
        rcd_output.setNamedVal<float>("event_wgt_triggers_dilepton_matched", event_wgt_triggers_dilepton_matched);
        rcd_output.setNamedVal<float>("event_wgt_SFs_btagging", event_wgt_SFs_btagging);
        rcd_output.setNamedVal("EventNumber", *ptr_EventNumber);
        if (!isData){
          rcd_output.setNamedVal("GenMET_pt", *ptr_genmet_pt);
          rcd_output.setNamedVal("GenMET_phi", *ptr_genmet_phi);
        }
        else{
          rcd_output.setNamedVal("RunNumber", *ptr_RunNumber);
          rcd_output.setNamedVal("LuminosityBlock", *ptr_LuminosityBlock);
        }
        rcd_output.setNamedVal("nak4jets_tight_selected", nak4jets_tight_selected);
        rcd_output.setNamedVal("nak4jets_tight_selected_btagged", nak4jets_tight_selected_btagged);
        rcd_output.setNamedVal("nleptons_tight", nleptons_tight);
        rcd_output.setNamedVal("nleptons_fakeable", nleptons_fakeable);
        rcd_output.setNamedVal("nleptons_loose", nleptons_loose);
        rcd_output.setNamedVal("iCRZ_Run2Analysis", iCRZ);
        rcd_output.setNamedVal("icat_Run2Analysis", icat);
        rcd_output.setNamedVal("HT_ak4jets", HT_ak4jets);
        rcd_output.setNamedVal("pTmiss", pTmiss);
        rcd_output.setNamedVal("phimiss", phimiss);

        // Record leptons
        {
#define BRANCH_VECTOR_COMMANDS \
          BRANCH_VECTOR_COMMAND(float, pt) \
          BRANCH_VECTOR_COMMAND(float, eta) \
          BRANCH_VECTOR_COMMAND(float, phi) \
          BRANCH_VECTOR_COMMAND(float, mass) \
          BRANCH_VECTOR_COMMAND(int, pdgId)

#define BRANCH_VECTOR_COMMAND(TYPE, NAME) std::vector<TYPE> leptons_##NAME;
          BRANCH_VECTOR_COMMANDS;
#undef BRANCH_VECTOR_COMMAND
          std::vector<bool> leptons_is_genmatched_prompt;
          std::vector<bool> leptons_is_genmatched_photonProduct;
          std::vector<int> leptons_genmatch_pdgId;

          for (auto const& part:leptons_tight){
            bool is_genmatched_prompt = leptons_isGenMatched_map.find(part)->second;
            bool is_genmatched_photonProduct = leptons_isGenMatched_photonProduct_map.find(part)->second;
            int genmatch_pdgId = leptons_genMatchPDGId_map.find(part)->second;

            leptons_is_genmatched_prompt.push_back(is_genmatched_prompt);
            leptons_is_genmatched_photonProduct.push_back(is_genmatched_photonProduct);
            leptons_genmatch_pdgId.push_back(genmatch_pdgId);
#define BRANCH_VECTOR_COMMAND(TYPE, NAME) leptons_##NAME.push_back(part->NAME());
            BRANCH_VECTOR_COMMANDS;
#undef BRANCH_VECTOR_COMMAND
          }

#define BRANCH_VECTOR_COMMAND(TYPE, NAME) rcd_output.setNamedVal(Form("leptons_%s", #NAME), leptons_##NAME);
          BRANCH_VECTOR_COMMANDS;
#undef BRANCH_VECTOR_COMMAND
          rcd_output.setNamedVal("leptons_is_genmatched_prompt", leptons_is_genmatched_prompt);
          rcd_output.setNamedVal("leptons_is_genmatched_photonProduct", leptons_is_genmatched_photonProduct);
          rcd_output.setNamedVal("leptons_genmatch_pdgId", leptons_genmatch_pdgId);

#undef BRANCH_VECTOR_COMMANDS
        }

        // Record ak4jets
        {
#define BRANCH_VECTOR_COMMANDS \
          BRANCH_VECTOR_COMMAND(float, pt) \
          BRANCH_VECTOR_COMMAND(float, eta) \
          BRANCH_VECTOR_COMMAND(float, phi) \
          BRANCH_VECTOR_COMMAND(float, mass)

#define BRANCH_VECTOR_COMMAND(TYPE, NAME) std::vector<TYPE> ak4jets_##NAME;
          BRANCH_VECTOR_COMMANDS;
#undef BRANCH_VECTOR_COMMAND
          std::vector<bool> ak4jets_is_tight;
          std::vector<bool> ak4jets_is_btagged;
          std::vector<int> ak4jets_hadronFlavour;
          std::vector<int> ak4jets_partonFlavour;

          for (auto const& jet:ak4jets){
            bool is_tight = ParticleSelectionHelpers::isTightJet(jet);
            bool is_btagged = jet->testSelectionBit(bit_preselection_btag);
            float pt = jet->pt();
            float eta = jet->eta();
            float phi = jet->phi();
            float mass = jet->mass();

#define BRANCH_VECTOR_COMMAND(TYPE, NAME) ak4jets_##NAME.push_back(jet->NAME());
            BRANCH_VECTOR_COMMANDS;
#undef BRANCH_VECTOR_COMMAND
            ak4jets_is_tight.push_back(is_tight);
            ak4jets_is_btagged.push_back(is_btagged);
            ak4jets_hadronFlavour.push_back(jet->extras.hadronFlavour);
            ak4jets_partonFlavour.push_back(jet->extras.partonFlavour);
          }

#define BRANCH_VECTOR_COMMAND(TYPE, NAME) rcd_output.setNamedVal(Form("ak4jets_%s", #NAME), ak4jets_##NAME);
          BRANCH_VECTOR_COMMANDS;
#undef BRANCH_VECTOR_COMMAND
          rcd_output.setNamedVal("ak4jets_is_tight", ak4jets_is_tight);
          rcd_output.setNamedVal("ak4jets_is_btagged", ak4jets_is_btagged);
          rcd_output.setNamedVal("ak4jets_hadronFlavour", ak4jets_hadronFlavour);
          rcd_output.setNamedVal("ak4jets_partonFlavour", ak4jets_partonFlavour);

#undef BRANCH_VECTOR_COMMANDS
        }

        if (firstOutputEvent){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_output.named##name_t##s.begin(); itb!=rcd_output.named##name_t##s.end(); itb++) tout->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_output.namedV##name_t##s.begin(); itb!=rcd_output.namedV##name_t##s.end(); itb++) tout->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_output.namedVV##name_t##s.begin(); itb!=rcd_output.namedVV##name_t##s.end(); itb++) tout->putBranch(itb->first, &(itb->second));
          SIMPLE_DATA_OUTPUT_DIRECTIVES;
          VECTOR_DATA_OUTPUT_DIRECTIVES;
          DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES;
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
        }
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_output.named##name_t##s.begin(); itb!=rcd_output.named##name_t##s.end(); itb++) tout->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_output.namedV##name_t##s.begin(); itb!=rcd_output.namedV##name_t##s.end(); itb++) tout->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=rcd_output.namedVV##name_t##s.begin(); itb!=rcd_output.namedVV##name_t##s.end(); itb++) tout->setVal(itb->first, &(itb->second));
        SIMPLE_DATA_OUTPUT_DIRECTIVES;
        VECTOR_DATA_OUTPUT_DIRECTIVES;
        DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES;
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

        tout->fill();
      }

      if (icat>=0){
        unsigned short idx_gencat = 0;
        if (hasGenPromptLepton_ChargeFlip) idx_gencat = 1;
        else if (hasGenMatchedPromptPhotonProduct) idx_gencat = 2;
        else if (hasNonPromptLepton) idx_gencat = 3;

        if (!hCats.empty()) hCats.at(idx_gencat)->Fill(static_cast<double>(icat-1)+0.5, static_cast<double>(iCRZ)+0.5, wgt);
        if (printObjInfo){
          if (dilepton_OS_DYCand_tight) IVYout
            << "OS Z candidate found to have "
            << "pt = " << dilepton_OS_DYCand_tight->pt()
            << ", mass = " << dilepton_OS_DYCand_tight->mass()
            << ", daughter array indices = ["
            << (std::abs(dilepton_OS_DYCand_tight->getDaughter(0)->pdgId())==13 ? "muon" : "electron") << " " << dilepton_OS_DYCand_tight->getDaughter(0)->getUniqueIdentifier()
            << ", " << (std::abs(dilepton_OS_DYCand_tight->getDaughter(1)->pdgId())==13 ? "muon" : "electron") << " " << dilepton_OS_DYCand_tight->getDaughter(1)->getUniqueIdentifier()
            << "]"
            << endl;
          if (dilepton_SS_tight) IVYout
            << "SS(TT) candidate found to have "
            << "pt = " << dilepton_SS_tight->pt()
            << ", mass = " << dilepton_SS_tight->mass()
            << ", daughter array indices = ["
            << (std::abs(dilepton_SS_tight->getDaughter(0)->pdgId())==13 ? "muon" : "electron") << " " << dilepton_SS_tight->getDaughter(0)->getUniqueIdentifier()
            << ", " << (std::abs(dilepton_SS_tight->getDaughter(1)->pdgId())==13 ? "muon" : "electron") << " " << dilepton_SS_tight->getDaughter(1)->getUniqueIdentifier()
            << "]"
            << endl;
        }
      }
      if (printObjInfo) IVYout << "Nb = " << nak4jets_tight_selected_btagged << ", Nj = " << nak4jets_tight_selected << ", iCRZ = " << iCRZ << ", icat = " << icat << endl;

      if (runSyncExercise && (!applyPreselection || iCRZ==1 || icat>0)) foutput_sync
        << *ptr_EventNumber << ","
        << ev << ","
        << eventmet->pt() << ","
        << eventmet->phi() << ","
        << nleptons_selected << ","
        << nleptons_fakeable + nleptons_tight << ","
        << nleptons_tight << ","
        << (dilepton_SS_tight ? "SS" : "!SS") << ","
        << HT_ak4jets << ","
        << nak4jets_tight_selected << ","
        << nak4jets_tight_selected_btagged << ","
        << (iCRZ ? "Z" : (icat==idx_CRW ? "W" : (icat==0 ? "N/A" : std::to_string(icat).data())))
        << endl;

      n_recorded++;

      if (firstOutputEvent) firstOutputEvent = false;
    }

    IVYout << "Number of events recorded: " << n_recorded << " / " << n_traversed << " / " << nEntries << endl;
    nevents_total_traversed += n_traversed;
  }
  seltracker.print();

  if (is_sim_data_flag==0 && nevents_total_traversed>0){
    double const sum_ntotal_req = (nchunks>0 ? (eventIndex_end - eventIndex_begin) : nevents_total);
    double const sum_ntotal_traversed = nevents_total_traversed;
    double const sum_ntotal_scale = sum_ntotal_req/sum_ntotal_traversed;
    IVYout << "Total number of events traversed / requested: " << sum_ntotal_traversed << " / " << sum_ntotal_req << endl;
    IVYout << "\t- Scaling the yield by " << sum_ntotal_scale << "..." << endl;
    for (auto& hCat:hCats) hCat->Scale(sum_ntotal_scale);
  }
  for (unsigned short igencat=0; igencat<hCats.size(); igencat++){
    auto const& hCat = hCats.at(igencat);
    IVYout << "Event counts for gen. cat. '" << gencats.at(igencat) << "':" << endl;
    double integral_error=0;
    for (int ix=1; ix<=hCat->GetNbinsX(); ix++) IVYout << "\t- " << hCat->GetXaxis()->GetBinLabel(ix) << ": " << hCat->GetBinContent(ix, 1) << " +- " << hCat->GetBinError(ix, 1) << endl;
    IVYout << "\t- CRZ: " << HelperFunctions::getHistogramIntegralAndError(hCat, 0, hCat->GetNbinsX()+1, 2, 2, false, &integral_error) << " +- " << integral_error << endl;
    IVYout << "\t- Total SR: " << HelperFunctions::getHistogramIntegralAndError(hCat, 1, hCat->GetNbinsX()-1, 1, 1, false, &integral_error) << " +- " << integral_error << endl;
    IVYout << "\t- Failed: " << hCat->GetBinContent(0, 1) << " +- " << hCat->GetBinError(0, 1) << endl;
  }

  if (tout_sync_objects){
    tout_sync_objects->writeToFile(foutput_sync_objects);
    delete tout_sync_objects;
    foutput_sync_objects->Close();
  }

  if (runSyncExercise){
    foutput_sync.close();
    SampleHelpers::splitFileAndAddForTransfer(stroutput_sync);
  }

  for (auto& hCat:hCats){
    foutput->WriteTObject(hCat);
    delete hCat;
  }
  if (tout){
    tout->writeToFile(foutput);
    delete tout;
  }
  foutput->Close();

  curdir->cd();
  for (auto& tin:tinlist) delete tin;

  // Split large files, and add them to the transfer queue from Condor to the target site
  // Does nothing if you are running the program locally because your output is already in the desired location.
  SampleHelpers::splitFileAndAddForTransfer(stroutput);

  // Close the output and error log files
  IVYout.close();
  IVYerr.close();

  return 0;
}

int main(int argc, char** argv){
  // argv[0]==[Executable name]
  constexpr int iarg_offset=1;

  // Switches that do not need =true.
  std::vector<std::string> const extra_argument_flags{
    "ignore_PU_rewgt",
    "ignore_L1Prefiring_rewgt",
    "no_FOs_phys_checks",
    "only_tight_dileptons",
    "produce_trees",
    "run_sync",
    "write_sync_objects",
    "force_sync_preselection",
    "shorthand_Run2_UL_proposal_config"
  };

  bool print_help=false, has_help=false;
  int ichunk=0, nchunks=0;
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
      else if (HelperFunctions::checkListVariable(extra_argument_flags, value)) extra_arguments.setNamedVal<bool>(value, true);
      else{
        IVYerr << "ERROR: Unknown argument " << value << endl;
        print_help=true;
      }
    }
    else if (HelperFunctions::checkListVariable(extra_argument_flags, wish)){
      bool tmpval;
      HelperFunctions::castStringToValue(value, tmpval);
      extra_arguments.setNamedVal(wish, tmpval);
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
    else if (wish=="muon_id" || wish=="electron_id" || wish=="btag" || wish=="output_file" || wish=="sync_evtno_filter_input") extra_arguments.setNamedVal(wish, value);
    else if (wish=="xsec"){
      if (xsec<0.) xsec = 1;
      xsec *= std::stod(value);
    }
    else if (wish=="BR"){
      if (xsec<0.) xsec = 1;
      xsec *= std::stod(value);
    }
    else if (wish=="nchunks") nchunks = std::stoi(value);
    else if (wish=="ichunk") ichunk = std::stoi(value);
    else if (wish=="minpt_jets" || wish=="minpt_bjets" || wish=="minpt_l3" || wish=="minpt_miss" || wish=="minHT_jets"){
      double tmpval;
      HelperFunctions::castStringToValue(value, tmpval);
      extra_arguments.setNamedVal(wish, tmpval);
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
  if (!print_help && (nchunks>0 && (ichunk<0 || ichunk>=nchunks))){
    IVYerr << "ERROR: Invalid arguments." << endl;
    print_help=true;
  }

  // Ensure that nchunks and ichunk make sense.
  if (nchunks<0) nchunks = 0;
  if (nchunks<=0) ichunk = 0;

  if (print_help){
    IVYout << argv[0] << " options:\n\n";
    IVYout << "- help: Prints this help message.\n";
    IVYout << "- dataset: Data set name. Mandatory.\n";
    IVYout << "- short_name: Process short name. Mandatory.\n";
    IVYout << "- period: Data period. Mandatory.\n";
    IVYout << "- input_tag: Version of the input skims. Mandatory.\n";
    IVYout << "- output_tag: Version of the output. Mandatory.\n";
    IVYout << "- xsec: Cross section value. Mandatory in the MC.\n";
    IVYout << "- BR: BR value. Mandatory in the MC.\n";
    IVYout << "- nchunks: Number of splits of the input. Optional. Input splitting is activated only if nchunks>0.\n";
    IVYout << "- ichunk: Index of split input split to run. Optional. The index should be in the range [0, nchunks-1].\n";
    IVYout
      << "- output_file: Identifier of the output file if different from 'short_name'. Optional.\n"
      << "  The full output file name is '[identifier](_[ichunk]_[nchunks]).root.'\n";
    IVYout << "- input_files: Input files to run. Optional. Default is to run on all files.\n";
    IVYout
      << "- muon_id: Can be 'Cutbased_Run2', 'TopMVA_Run2', or 'TopMVAv2_Run2'.\n"
      << "  Default is whatever is in MuonSelectionHelpers (currently 'Cutbased_Run2') if no value is given.\n";
    IVYout
      << "- electron_id: Can be 'Cutbased_Run2', 'TopMVA_Run2', or 'TopMVAv2_Run2'.\n"
      << "  Default is whatever is in ElectronSelectionHelpers (currently 'Cutbased_Run2') if no value is given.\n";
    IVYout << "- btag: Name of the b-tagging WP. Can be 'medium' or 'loose' (case-insensitive). Default='medium'.\n";
    IVYout << "- minpt_jets: Minimum pT of jets in units of GeV. Default=40.\n";
    IVYout << "- minpt_bjets: Minimum pT of b-tagged jets in units of GeV. Default=25.\n";
    IVYout << "- minpt_l3: Minimum pT of third lepton in units of GeV. Default=20.\n";
    IVYout << "- minpt_miss: Minimum pTmiss in units of GeV. Default=50.\n";
    IVYout << "- minHT_jets: Minimum HT over ak4 jets in units of GeV. Default=300.\n";
    IVYout << "- ignore_PU_rewgt: Flag to turn off PU reweighting. Optional. Default is to include it.\n";
    IVYout << "- ignore_L1Prefiring_rewgt: Flag to turn off L1 prefiring reweighting in 2016 and 2017. Optional. Default is to include it.\n";
    IVYout << "- no_FOs_phys_checks: Flag to turn off the usage of fakeable leptons in physics object checks such as jet cleaning. Optional. Default is to take FOs into account.\n";
    IVYout << "- only_tight_dileptons: Flag to use only tight leptons in general dilepton pair construction. Optional. Default is to use all loose, fakeable, or tight leptons.\n";
    IVYout << "- produce_trees: Flag to switch to tree output instead of histograms. Optional. Default is to produce histograms.\n";
    IVYout << "  Note that tree production keeps all events, so one needs to use the selection bit map.\n";
    IVYout << "- run_sync: Turn on synchronization output. Optional. Default is to run without synchronization output.\n";
    IVYout << "- write_sync_objects: Create a file that contains the info. for all leptons and jets, and event identifiers. Ignored if run_sync=false. Optional. Default is to not produce such a file.\n";
    IVYout << "- force_sync_preselection: When sync. mode is on, also apply SR/CR preselection.  Ignored if run_sync=false. Optional. Default is to run without preselection.\n";
    IVYout << "- sync_evtno_filter_input: Input file to specify filter on event numbers. Optional. Processed only if run_sync=true.\n";
    IVYout
      << "- shorthand_Run2_UL_proposal_config: Shorthand flag for the switches for the Run 2 UL analysis proposal:\n"
      << "  * muon_id='TopMVA_Run2'\n"
      << "  * electron_id='TopMVA_Run2'\n"
      << "  * btag='loose'\n"
      << "  * minpt_jets=25\n"
      << "  * minpt_bjets=25\n"
      << "  * minpt_l3=10\n"
      << "  The use of this shorthand will ignore the user-defined setting of these options above.\n";

    IVYout << endl;
    return (has_help ? 0 : 1);
  }

  SampleHelpers::configure(str_period, Form("skims:%s", str_tag.data()), HostHelpers::kUCSDT2);

  return ScanChain(str_outtag, str_dset, str_proc, xsec, ichunk, nchunks, extra_arguments);
}
