#include <string>
#include <iostream>
#include <iomanip>
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
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
#include "BtagScaleFactorHandler.h"
#include "SamplesCore.h"
#include "FourTopTriggerHelpers.h"
#include "DileptonHandler.h"
#include "InputChunkSplitter.h"
#include "SplitFileAndAddForTransfer.h"


using namespace std;
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


int ScanChain(std::string const& strdate, std::string const& dset, std::string const& proc, double const& xsec, int const& ichunk, int const& nchunks, SimpleEntry const& extra_arguments){
  bool const isCondorRun = SampleHelpers::checkRunOnCondor();
  if (!isCondorRun) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  TDirectory* curdir = gDirectory;

  // Systematics hypothesis: For now, keep as "Nominal" for the purposes of this example. This could be made an argument, or the looper could go through each relevant systematic hypothesis in each event.
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  // Data period quantities
  auto const& theDataPeriod = SampleHelpers::getDataPeriod();
  auto const& theDataYear = SampleHelpers::getDataYear();

  // Keep this as true for checks as cleaning jets from fakeable leptons
  constexpr bool useFakeableIdForPhysicsChecks = true;
  ParticleSelectionHelpers::setUseFakeableIdForPhysicsChecks(useFakeableIdForPhysicsChecks);

  // b-tagging |eta| threshold should be 2.4 in year<=2016 (Phase-0 tracker acceptance) and 2.5 afterward (Phase-1 tracker)
  float const absEtaThr_ak4jets = (theDataYear<=2016 ? AK4JetSelectionHelpers::etaThr_btag_Phase0Tracker : AK4JetSelectionHelpers::etaThr_btag_Phase1Tracker);

  // Integrated luminosity for the data period
  double const lumi = SampleHelpers::getIntegratedLuminosity(theDataPeriod);
  // If theDataPeriod is not a specific era, print the data eras that are actually contained.
  // lumi is computed for the contained eras only.
  IVYout << "Valid data periods for " << theDataPeriod << ": " << SampleHelpers::getValidDataPeriods() << endl;
  IVYout << "Integrated luminosity: " << lumi << endl;

  // This is the output directory.
  // Output should always be recorded as if you are running the job locally.
  // We will inform the Condor job later on that some files would need transfer if we are running on Condor.
  TString coutput_main = TString("output/Analysis_TTJetRadiation/") + strdate.data() + "/" + SampleHelpers::getDataPeriod();
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

  // Turn on individual files
  std::string input_files;
  extra_arguments.getNamedVal("input_files", input_files);

  // Shorthand option for the Run 2 UL analysis proposal
  bool use_shorthand_Run2_UL_proposal_config;
  extra_arguments.getNamedVal("shorthand_Run2_UL_proposal_config", use_shorthand_Run2_UL_proposal_config);

  // Options to set alternative muon and electron IDs, or b-tagging WP
  std::string muon_id_name;
  std::string electron_id_name;
  if (use_shorthand_Run2_UL_proposal_config){
    muon_id_name = electron_id_name = "TopMVA_Run2";
  }
  else{
    extra_arguments.getNamedVal("muon_id", muon_id_name);
    extra_arguments.getNamedVal("electron_id", electron_id_name);
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

  // Trigger configuration
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton{
    TriggerHelpers::kMuEle
  };
  // These PFHT triggers were used in the 2016 analysis. They are a bit more efficient.
  if (SampleHelpers::getDataYear()==2016){
    requiredTriggers_Dilepton = std::vector<TriggerHelpers::TriggerType>{
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

  // SF handlers
  BtagScaleFactorHandler btagSFHandler;

  // These are called handlers, but they are more like helpers.
  DileptonHandler dileptonHandler;
  ParticleDisambiguator particleDisambiguator;

  // Some advanced event filters
  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(has_multiple_dsets);
  eventFilter.setCheckHLTPathRunRanges(true);

  curdir->cd();

  // We conclude the setup of event processing specifications and move on to I/O configuration next.

  // Open the output ROOT file
  SimpleEntry rcd_output;
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("Events");
  tout->setAutoSave(0);
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
    IVYout << "Acquired a sum of weights of " << sum_wgts << ". Overall normalization will be " << norm_scale << "." << endl;

    nevents_total += tin->getNEvents();

    curdir->cd();

    // Book the necessary branches
    genInfoHandler.bookBranches(tin);
    simEventHandler.bookBranches(tin);
    eventFilter.bookBranches(tin);
    muonHandler.bookBranches(tin);
    electronHandler.bookBranches(tin);
    jetHandler.bookBranches(tin);

    // Book a few additional branches
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, NANONAME) tin->bookBranch<TYPE>(#NANONAME, 0);
    EVENT_VARIABLE;
    if (!isData){
      tin->bookBranch<float>("GenMET_pt", 0);
      tin->bookBranch<float>("GenMET_phi", 0);
    }
    else{
      RUNLUMI_VARIABLES;
    }
#undef RUNLUMIEVENT_VARIABLE
  }

  curdir->cd();

  // Prepare to loop!
  // Keep track of sums of weights
  SelectionTracker seltracker;

  // Keep track of the traversed events
  bool firstOutputEvent = true;
  int eventIndex_begin = -1;
  int eventIndex_end = -1;
  int eventIndex_tracker = 0;
  splitInputEventsIntoChunks((is_sim_data_flag==1), nevents_total, ichunk, nchunks, eventIndex_begin, eventIndex_end);

  for (auto const& tin:tinlist){
    if (SampleHelpers::doSignalInterrupt==1) break;

    auto const& norm_scale = tin_normScale_map.find(tin)->second;
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

#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, NANONAME) TYPE* ptr_##NAME = nullptr;
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    float* ptr_genmet_pt = nullptr;
    float* ptr_genmet_phi = nullptr;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, NANONAME) tin->getValRef(#NANONAME, ptr_##NAME);
    EVENT_VARIABLE;
    if (!isData){
      tin->getValRef("GenMET_pt", ptr_genmet_pt);
      tin->getValRef("GenMET_phi", ptr_genmet_phi);
    }
    else{
      RUNLUMI_VARIABLES;
    }
#undef RUNLUMIEVENT_VARIABLE

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

      eventIndex_tracker++;
      if (!doAccumulate) continue;

      tin->getEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);
      n_traversed++;

      genInfoHandler.constructGenInfo();
      auto const& genInfo = genInfoHandler.getGenInfo();
      auto const& genparticles = genInfoHandler.getGenParticles();

      std::vector<GenParticleObject*> genmatch_promptleptons;
      for (auto const& part:genparticles){
        if (part->status()!=1) continue; // Only final states
        if (!part->extras.isPromptFinalState) continue; // Only prompt final states
        bool const isLepton = std::abs(part->pdgId())==11 || std::abs(part->pdgId())==13;
        if (!isLepton) continue; // Only leptons
        genmatch_promptleptons.push_back(part);
      }

      simEventHandler.constructSimEvent();

      double wgt = 1;
      if (!isData){
        // Regular gen. weight
        double genwgt = 1;
        genwgt = genInfo->getGenWeight(theGlobalSyst);

        // PU reweighting (time-dependent)
        double puwgt = 1;
        puwgt = simEventHandler.getPileUpWeight(theGlobalSyst);

        wgt = genwgt * puwgt;

        // Add L1 prefiring weight for 2016 and 2017
        wgt *= simEventHandler.getL1PrefiringWeight(theGlobalSyst);
      }
      // Overall sample normalization
      wgt *= norm_scale;

      muonHandler.constructMuons();
      electronHandler.constructElectrons();
      jetHandler.constructJetMET(&genInfoHandler, &simEventHandler, theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, nullptr, &jetHandler);

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
        genmatch_promptleptons.begin(), genmatch_promptleptons.end(),
        lepton_genmatchpart_map
      );

      // Keep track of leptons
      std::vector<ParticleObject*> leptons_tight;

      std::vector<MuonObject*> muons_selected;
      std::vector<MuonObject*> muons_tight;
      std::vector<MuonObject*> muons_fakeable;
      std::vector<MuonObject*> muons_loose;
      for (auto const& part:muons){
        if (ParticleSelectionHelpers::isTightParticle(part)){
          muons_tight.push_back(part);
          leptons_tight.push_back(dynamic_cast<ParticleObject*>(part));
        }
        else if (ParticleSelectionHelpers::isFakeableParticle(part)) muons_fakeable.push_back(part);
        else if (ParticleSelectionHelpers::isLooseParticle(part)) muons_loose.push_back(part);
      }
      HelperFunctions::appendVector(muons_selected, muons_tight);
      HelperFunctions::appendVector(muons_selected, muons_fakeable);
      HelperFunctions::appendVector(muons_selected, muons_loose);

      std::vector<ElectronObject*> electrons_selected;
      std::vector<ElectronObject*> electrons_tight;
      std::vector<ElectronObject*> electrons_fakeable;
      std::vector<ElectronObject*> electrons_loose;
      for (auto const& part:electrons){
        if (ParticleSelectionHelpers::isTightParticle(part)){
          electrons_tight.push_back(part);
          leptons_tight.push_back(dynamic_cast<ParticleObject*>(part));
        }
        else if (ParticleSelectionHelpers::isFakeableParticle(part)) electrons_fakeable.push_back(part);
        else if (ParticleSelectionHelpers::isLooseParticle(part)) electrons_loose.push_back(part);
      }
      HelperFunctions::appendVector(electrons_selected, electrons_tight);
      HelperFunctions::appendVector(electrons_selected, electrons_fakeable);
      HelperFunctions::appendVector(electrons_selected, electrons_loose);

      unsigned int const nleptons_tight = muons_tight.size() + electrons_tight.size();
      unsigned int const nleptons_fakeable = muons_fakeable.size() + electrons_fakeable.size();
      unsigned int const nleptons_loose = muons_loose.size() + electrons_loose.size();
      unsigned int const nleptons_selected = nleptons_tight + nleptons_fakeable + nleptons_loose;

      ParticleObjectHelpers::sortByGreaterPt(leptons_tight);

      // Keep track of jets
      double event_wgt_SFs_btagging = 1;
      auto const& ak4jets = jetHandler.getAK4Jets();
      std::vector<AK4JetObject*> ak4jets_tight_selected;
      for (auto const& jet:ak4jets){
        bool is_btaggable = jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTaggable);
        constexpr bool is_clean = true;
        if (is_btaggable){
          ak4jets_tight_selected.push_back(jet);
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
        }
      }
      unsigned int const njets_tight = ak4jets_tight_selected.size();

      // MET info
      auto const& eventmet = jetHandler.getPFMET();

      // BEGIN PRESELECTION
      seltracker.accumulate("Full sample", wgt);

      if (nleptons_fakeable+nleptons_loose>0) continue;
      seltracker.accumulate("Has exactly 0 loose and 0 fakeable leptons", wgt);

      if (electrons_tight.size()!=1) continue;
      seltracker.accumulate("Has exactly one tight electron", wgt);

      if (muons_tight.size()!=1) continue;
      seltracker.accumulate("Has exactly one tight muon", wgt);

      // check oppososite charge
      if (electrons_tight.front()->pdgId() * muons_tight.front()->pdgId() > 0) continue;
      seltracker.accumulate("Has opposite charge", wgt);

      // check lep1 pt>25 lep2 pt>20
      if (leptons_tight.front()->pt() < 25. || leptons_tight.back()->pt() < 20.) continue;
      seltracker.accumulate("Has pTl1>=25 GeV and pTl2>=20 GeV", wgt);

      // check MET
      float const pTmiss = eventmet->pt();
      float const phimiss = eventmet->phi();
      if (pTmiss<25.f) continue;
      seltracker.accumulate("pTmiss>=25 GeV", wgt);

      // Put event filters to the last because data has unique event tracking enabled.
      eventFilter.constructFilters(&simEventHandler);
      constexpr bool pass_HEMveto = true;
      //bool const pass_HEMveto = eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets);
      //bool const pass_HEMveto = eventFilter.test2018HEMFilter(&simEventHandler, &electrons, nullptr, &ak4jets);
      if (!pass_HEMveto) continue;
      seltracker.accumulate("Pass HEM veto", wgt);

      if (!eventFilter.passMETFilters()) continue;
      seltracker.accumulate("Pass MET filters", wgt);

      if (!eventFilter.isUniqueDataEvent()) continue;
      seltracker.accumulate("Pass unique event check", wgt);

      // Triggers
      float event_wgt_triggers_dilepton = eventFilter.getTriggerWeight(hltnames_Dilepton);
      if (event_wgt_triggers_dilepton==0.f) continue;
      seltracker.accumulate("Pass any trigger", wgt);

      float event_wgt_triggers_dilepton_matched = eventFilter.getTriggerWeight(
        triggerPropsCheckList_Dilepton,
        &muons, &electrons, nullptr, &ak4jets, nullptr, nullptr
      );
      bool const pass_triggers_dilepton_matched = (event_wgt_triggers_dilepton_matched!=0.f);
      // Do not skip the event. Instead, record a flag for HLT object matching.
      rcd_output.setNamedVal("pass_triggers_dilepton_matched", pass_triggers_dilepton_matched);
      seltracker.accumulate("Pass triggers after matching", wgt*static_cast<double>(pass_triggers_dilepton_matched));


      /*************************************************/
      /* NO MORE CALLS TO SELECTION BEYOND THIS POINT! */
      /*************************************************/
      // Write output
      rcd_output.setNamedVal<float>("event_wgt", wgt);
      rcd_output.setNamedVal<float>("event_wgt_triggers_dilepton", event_wgt_triggers_dilepton);
      rcd_output.setNamedVal<float>("event_wgt_triggers_dilepton_matched", event_wgt_triggers_dilepton_matched);
      rcd_output.setNamedVal<float>("event_wgt_SFs_btagging", event_wgt_SFs_btagging);
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, NANONAME) rcd_output.setNamedVal(#NAME, *ptr_##NAME);
      EVENT_VARIABLE;
      if (!isData){
        rcd_output.setNamedVal("GenMET_pt", *ptr_genmet_pt);
        rcd_output.setNamedVal("GenMET_phi", *ptr_genmet_phi);
        rcd_output.setNamedVal("genTtbarId", genInfo->extras.genTtbarId);
      }
      else{
        RUNLUMI_VARIABLES;
      }
#undef RUNLUMIEVENT_VARIABLE
      rcd_output.setNamedVal("nak4jets_tight_pt25", njets_tight);
      rcd_output.setNamedVal("pTmiss", pTmiss);
      rcd_output.setNamedVal("phimiss", phimiss);

      {
        // make vectors of pt, eta, phi, mass, pdgId for leptons
        std::vector<bool> leptons_is_genmatched_prompt;
        std::vector<int> leptons_pdgId;
        std::vector<float> leptons_pt;
        std::vector<float> leptons_eta;
        std::vector<float> leptons_phi;
        std::vector<float> leptons_mass;
        for (auto const& part:leptons_tight){
          leptons_pdgId.push_back(part->pdgId());
          leptons_pt.push_back(part->pt());
          leptons_eta.push_back(part->eta());
          leptons_phi.push_back(part->phi());
          leptons_mass.push_back(part->mass());

          auto it_genmatch = lepton_genmatchpart_map.find(part);
          leptons_is_genmatched_prompt.push_back((it_genmatch!=lepton_genmatchpart_map.end() && it_genmatch->second!=nullptr));
        }

        // setNamedVal for pt, eta, phi, mass, pdgId for leptons
        rcd_output.setNamedVal("leptons_is_genmatched_prompt", leptons_is_genmatched_prompt);
        rcd_output.setNamedVal("leptons_pdgId", leptons_pdgId);
        rcd_output.setNamedVal("leptons_pt", leptons_pt);
        rcd_output.setNamedVal("leptons_eta", leptons_eta);
        rcd_output.setNamedVal("leptons_phi", leptons_phi);
        rcd_output.setNamedVal("leptons_mass", leptons_mass);
      }

      {
        // make vectors of pt, eta, phi, mass, is_btagged, and hadronFlavour for jets
        std::vector<unsigned char> ak4jets_pass_btagging;
        std::vector<int> ak4jets_hadronFlavour;
        std::vector<float> ak4jets_pt;
        std::vector<float> ak4jets_eta;
        std::vector<float> ak4jets_phi;
        std::vector<float> ak4jets_mass;
        for (auto const& jet:ak4jets_tight_selected){
          bool is_btagged_loose = jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Loose);
          bool is_btagged_medium = jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Medium);
          bool is_btagged_tight = jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Tight);
          unsigned char pass_btagging = int(is_btagged_loose) + int(is_btagged_medium) + int(is_btagged_tight);
          ak4jets_pass_btagging.push_back(pass_btagging);
          ak4jets_hadronFlavour.push_back(jet->extras.hadronFlavour);
          ak4jets_pt.push_back(jet->pt());
          ak4jets_eta.push_back(jet->eta());
          ak4jets_phi.push_back(jet->phi());
          ak4jets_mass.push_back(jet->mass());
        }

        // setNamedVal for pt, eta, phi, mass, is_btagged, and hadronFlavour for jets
        rcd_output.setNamedVal("ak4jets_pass_btagging", ak4jets_pass_btagging);
        rcd_output.setNamedVal("ak4jets_hadronFlavour", ak4jets_hadronFlavour);
        rcd_output.setNamedVal("ak4jets_pt", ak4jets_pt);
        rcd_output.setNamedVal("ak4jets_eta", ak4jets_eta);
        rcd_output.setNamedVal("ak4jets_phi", ak4jets_phi);
        rcd_output.setNamedVal("ak4jets_mass", ak4jets_mass);
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

      // Record whatever is in rcd_output into the tree.
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
      n_recorded++;

      if (firstOutputEvent) firstOutputEvent = false;
    }

    IVYout << "Number of events recorded: " << n_recorded << " / " << n_traversed << " / " << nEntries << endl;
  }
  seltracker.print();

  tout->writeToFile(foutput);
  delete tout;
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
  constexpr int iarg_offset=1; // argv[0]==[Executable name]

  // Switches that do not need =true.
  std::vector<std::string> const extra_argument_flags{
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
    HelperFunctions::splitOption(strarg, wish, value, '=');

    if (wish.empty()){
      if (value=="help"){ print_help=has_help=true; }
      else if (HelperFunctions::checkListVariable(extra_argument_flags, value)) extra_arguments.setNamedVal<bool>(value, true);
      else{
        IVYerr << "ERROR: Unknown argument " << value << endl;
        print_help=true;
      }
    }
    else if (HelperFunctions::checkListVariable(extra_argument_flags, wish)){
      // Case where the user runs '[executable] [flag]=true/false' as a variation of the above.
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
    else if (wish=="muon_id" || wish=="electron_id" || wish=="output_file") extra_arguments.setNamedVal(wish, value);
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
      << "- shorthand_Run2_UL_proposal_config: Shorthand flag for the switches for the Run 2 UL analysis proposal:\n"
      << "  * muon_id='TopMVA_Run2'\n"
      << "  * electron_id='TopMVA_Run2'\n"
      << "  The use of this shorthand will ignore the user-defined setting of these options above.\n";
    IVYout
      << "- muon_id: Can be 'Cutbased_Run2', 'TopMVA_Run2', or 'TopMVAv2_Run2'.\n"
      << "  Default is whatever is in MuonSelectionHelpers (currently 'Cutbased_Run2') if no value is given.\n";
    IVYout
      << "- electron_id: Can be 'Cutbased_Run2', 'TopMVA_Run2', or 'TopMVAv2_Run2'.\n"
      << "  Default is whatever is in ElectronSelectionHelpers (currently 'Cutbased_Run2') if no value is given.\n";

    IVYout << endl;
    return (has_help ? 0 : 1);
  }

  SampleHelpers::configure(str_period, Form("skims:%s", str_tag.data()), HostHelpers::kUCSDT2);

  return ScanChain(str_outtag, str_dset, str_proc, xsec, ichunk, nchunks, extra_arguments);
}
