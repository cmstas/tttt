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
#include "InputChunkSplitter.h"
#include "SplitFileAndAddForTransfer.h"


using namespace std;
using namespace HelperFunctions;
using namespace IvyStreamHelpers;
using namespace SystematicsHelpers;


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

int ScanChain(std::string const& strdate, std::string const& dset, std::string const& proc, int const& ichunk, int const& nchunks, SimpleEntry const& extra_arguments){
  bool const isCondorRun = SampleHelpers::checkRunOnCondor();
  if (!isCondorRun) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  TDirectory* curdir = gDirectory;

  // Data period quantities
  auto const& theDataPeriod = SampleHelpers::getDataPeriod();
  auto const& theDataYear = SampleHelpers::getDataYear();

  // Configure analysis-specific stuff
  constexpr bool useFakeableIdForPhysicsChecks = true;
  ParticleSelectionHelpers::setUseFakeableIdForPhysicsChecks(useFakeableIdForPhysicsChecks);

  float const absEtaThr_ak4jets = (theDataYear<=2016 ? AK4JetSelectionHelpers::etaThr_btag_Phase0Tracker : AK4JetSelectionHelpers::etaThr_btag_Phase1Tracker);

  // This is the output directory.
  // Output should always be recorded as if you are running the job locally.
  // We will inform the Condor job later on that some files would need transfer if we are running on Condor.
  TString coutput_main = TString("output/SimJetEffs/") + strdate.data() + "/" + theDataPeriod;
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
  if (has_multiple_dsets){
    IVYerr << "Only a single data set should be passed." << endl;
    assert(0);
  }

  // Configure output file name
  std::string output_file;
  extra_arguments.getNamedVal("output_file", output_file);
  if (output_file=="") output_file = proc;
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
  bool use_shorthand_Run2_UL_proposal_config = false;
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

  // No triggers are requested.
  std::vector<TriggerHelpers::TriggerType> const requestedTriggers;

  // Declare handlers
  GenInfoHandler genInfoHandler;
  SimEventHandler simEventHandler;
  EventFilterHandler eventFilter(requestedTriggers);
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;

  // These are called handlers, but they are more like helpers.
  DileptonHandler dileptonHandler;
  ParticleDisambiguator particleDisambiguator;

  // Some advanced event filters
  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(false);
  eventFilter.setCheckHLTPathRunRanges(false);

  curdir->cd();

  // We conclude the setup of event processing specifications and move on to I/O configuration next.

  // Open the output ROOT file
  SimpleEntry rcd_output;
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
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

  std::vector<SystematicsHelpers::SystematicVariationTypes> const allowedSysts_GenSim{
    sNominal,
    ePUDn, ePUUp
  };
  std::vector<SystematicsHelpers::SystematicVariationTypes> const allowedSysts_MomScale{
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp
  };
  /*
  std::vector<SystematicsHelpers::SystematicVariationTypes> const allowedSysts_JetEff{
    sNominal,
    ePUJetIdEffDn, ePUJetIdEffUp
  };
  */


  int nevents_total = 0;
  std::unordered_map< BaseTree*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> > tin_normScale_map;
  std::vector<BaseTree*> tinlist; tinlist.reserve(dset_proc_pairs.size());
  for (auto const& dset_proc_pair:dset_proc_pairs){
    TString strinput = SampleHelpers::getInputDirectory() + "/" + strinputdpdir + "/" + dset_proc_pair.second.data();
    TString cinput = (input_files=="" ? strinput + "/*.root" : strinput + "/" + input_files.data());
    IVYout << "Accessing input files " << cinput << "..." << endl;
    TString const sid = SampleHelpers::getSampleIdentifier(dset_proc_pair.first);
    bool const isData = SampleHelpers::checkSampleIsData(sid);
    if (isData){
      IVYerr << "This executable should only run on the MC." << endl;
      assert(0);
    }
    BaseTree* tin = new BaseTree(cinput, "Events", "", (isData ? "" : "Counters"));
    if (!tin->isValid()){
      IVYout << "An error occured while acquiring the input from " << cinput << ". Aborting..." << endl;
      assert(0);
    }
    tin->sampleIdentifier = sid;

    tin_normScale_map[tin] = std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double>();
    auto& normScales = tin_normScale_map.find(tin)->second;
    {
      TH2D* hCounters = (TH2D*) tin->getCountersHistogram();
      double const sumN = hCounters->GetBinContent(0, 2);
      for (auto const& ss:allowedSysts_GenSim){
        int ix = 1;
        switch (ss){
        case ePUDn:
          ix = 2;
          break;
        case ePUUp:
          ix = 3;
          break;
        default:
          break;
        }
        double sum_wgts = hCounters->GetBinContent(ix, 1);;
        normScales[ss] = (sum_wgts!=0. ? sumN / sum_wgts : 1.);
        IVYout << "Acquired a sum of weights of " << sum_wgts << " for systematic " << SystematicsHelpers::getSystName(ss) << ". Overall normalization will be " << normScales[ss] << "." << endl;
      }
    }

    // Add the tree to the list of trees to process
    tinlist.push_back(tin);

    nevents_total += tin->getNEvents();

    curdir->cd();

    // Book the necessary branches
    genInfoHandler.bookBranches(tin);
    simEventHandler.bookBranches(tin);
    eventFilter.bookBranches(tin);
    muonHandler.bookBranches(tin);
    electronHandler.bookBranches(tin);
    jetHandler.bookBranches(tin);
    isotrackHandler.bookBranches(tin);
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
  splitInputEventsIntoChunks(false, nevents_total, ichunk, nchunks, eventIndex_begin, eventIndex_end);

  for (auto const& tin:tinlist){
    if (SampleHelpers::doSignalInterrupt==1) break;

    auto const& normScales = tin_normScale_map.find(tin)->second;
    constexpr bool isData = false;

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
    isotrackHandler.wrapTree(tin);

    unsigned int n_traversed = 0;
    unsigned int n_recorded = 0;
    int nEntries = tin->getNEvents();
    IVYout << "Looping over " << nEntries << " events from " << tin->sampleIdentifier << "..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      bool doAccumulate = (
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

      simEventHandler.constructSimEvent();

      double wgt_nominal = 1;
      for(auto const& syst:allowedSysts_GenSim){
        double wgt = normScales.find(syst)->second;
        wgt *= genInfo->getGenWeight(syst);
        wgt *= simEventHandler.getPileUpWeight(syst);
        wgt *= simEventHandler.getL1PrefiringWeight(syst);
        if (syst == sNominal) wgt_nominal = wgt;
        rcd_output.setNamedVal(Form("event_wgt_%s", SystematicsHelpers::getSystName(syst).data()), wgt);
      }

      // ak4jet output variables
#define SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(unsigned short, ak4jets, btagcat) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, pt) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, eta) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, phi) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, mass) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, genjet_pt) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, genjet_eta) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, genjet_phi) \
      SYNC_OBJ_BRANCH_VECTOR_COMMAND(float, ak4jets, genjet_mass) \
      AK4JET_EXTRA_INPUT_VARIABLES
      // All outputs
#define SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS \
      SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) std::vector<TYPE> COLLNAME##_##NAME;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, ak4jets, NAME)
      SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND

      // Since we are not processing data, we can put event filters before anything else.
      eventFilter.constructFilters(&simEventHandler);

      seltracker.accumulate("Full sample", wgt_nominal);
      if (!eventFilter.passMETFilters()) continue; // Test for MET filters
      seltracker.accumulate("Pass MET filters", wgt_nominal);
      if (!eventFilter.isUniqueDataEvent()) continue; // Test if the data event is unique (i.e., dorky). Does not do anything in the MC.
      seltracker.accumulate("Pass unique event check", wgt_nominal);

      // Loop over jet momentum systematics
      bool passAnyMomSyst = false;
      for (auto const& syst:allowedSysts_MomScale){
        auto systname = SystematicsHelpers::getSystName(syst);

        muonHandler.constructMuons();
        electronHandler.constructElectrons();
        jetHandler.constructJetMET(&genInfoHandler, &simEventHandler, syst);

        particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, nullptr, &jetHandler);

        auto const& muons = muonHandler.getProducts();
        std::vector<MuonObject*> muons_selected;
        std::vector<MuonObject*> muons_tight;
        std::vector<MuonObject*> muons_fakeable;
        std::vector<MuonObject*> muons_loose;
        for (auto const& part:muons){
          if (ParticleSelectionHelpers::isTightParticle(part)) muons_tight.push_back(part);
          else if (ParticleSelectionHelpers::isFakeableParticle(part)) muons_fakeable.push_back(part);
          else if (ParticleSelectionHelpers::isLooseParticle(part)) muons_loose.push_back(part);
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
          if (ParticleSelectionHelpers::isTightParticle(part)) electrons_tight.push_back(part);
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

        auto const& ak4jets = jetHandler.getAK4Jets();
        unsigned int nak4jets_selected = 0;
        for (auto const& jet:ak4jets){
          if (!jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTaggable)) continue;

          float pt = jet->pt();
          float eta = jet->eta();
          float phi = jet->phi();
          float mass = jet->mass();

          unsigned short btagcat = 0;
          if (jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Loose)) btagcat++;
          if (jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Medium)) btagcat++;
          if (jet->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Tight)) btagcat++;

          float genjet_pt = -1;
          float genjet_eta = 0;
          float genjet_phi = 0;
          float genjet_mass = 0;
          {
            GenJetObject* genjet = nullptr;
            for (auto const& part:jet->getMothers()){
              genjet = dynamic_cast<GenJetObject*>(part);
              if (genjet) break;
            }
            if (genjet){
              genjet_pt = genjet->pt();
              genjet_eta = genjet->eta();
              genjet_phi = genjet->phi();
              genjet_mass = genjet->mass();
            }
          }

#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.push_back(NAME);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) ak4jets_##NAME.push_back(jet->extras.NAME);
          SYNC_AK4JETS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND

          nak4jets_selected++;
        }

        // PRESELECTION
        bool hasFakeableLepton = ((nleptons_tight + nleptons_fakeable)>0);
        if (syst==sNominal) seltracker.accumulate("Has at least one fakeable lepton", wgt_nominal*static_cast<float>(hasFakeableLepton));

        bool hasTightBJet = (nak4jets_selected>0);
        if (syst==sNominal) seltracker.accumulate("Has at least one tight, b-taggable jet", wgt_nominal*static_cast<float>((hasFakeableLepton && hasTightBJet)));

        passAnyMomSyst |= hasFakeableLepton && hasTightBJet;

        // No need to apply a HEM veto in the MC

        /****************************************************/
        /* NO MORE CALLS TO PRESELECTION BEYOND THIS POINT! */
        /****************************************************/
        // Write output
        rcd_output.setNamedVal(Form("nleptons_tight_%s", systname.data()), nleptons_tight);
        rcd_output.setNamedVal(Form("nleptons_fakeable_%s", systname.data()), nleptons_fakeable);
        rcd_output.setNamedVal(Form("nleptons_loose_%s", systname.data()), nleptons_loose);

#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) rcd_output.setNamedVal(Form("%s_%s_%s", #COLLNAME, #NAME, systname.data()), COLLNAME##_##NAME);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, ak4jets, NAME)
        SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND

        // Reset the caches because the next systematic should re-reconstruct the particles from scratch
        muonHandler.resetCache();
        electronHandler.resetCache();
        jetHandler.resetCache();

        // Clear collections
#define SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, COLLNAME, NAME) COLLNAME##_##NAME.clear();
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) SYNC_OBJ_BRANCH_VECTOR_COMMAND(TYPE, ak4jets, NAME)
        SYNC_ALLOBJS_BRANCH_VECTOR_COMMANDS;
#undef AK4JET_VARIABLE
#undef SYNC_OBJ_BRANCH_VECTOR_COMMAND
      }
      if (!passAnyMomSyst) continue;

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
  // argv[0]==[Executable name]
  constexpr int iarg_offset=1;

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
    else if (wish=="muon_id" || wish=="electron_id" || wish=="output_file") extra_arguments.setNamedVal(wish, value);
    else if (wish=="nchunks") nchunks = std::stoi(value);
    else if (wish=="ichunk") ichunk = std::stoi(value);
    else if (wish=="xsec" || wish=="BR"){ /* Do nothing. */}
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

    IVYout << endl;
    return (has_help ? 0 : 1);
  }

  SampleHelpers::configure(str_period, Form("skims:%s", str_tag.data()), HostHelpers::kUCSDT2);

  return ScanChain(str_outtag, str_dset, str_proc, ichunk, nchunks, extra_arguments);
}
