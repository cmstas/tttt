#include <cassert>
#include "HostHelpersCore.h"
#include "BtagHelpers.h"
#include "SamplesCore.h"
#include "IvyFramework/IvyDataTools/interface/HelperFunctions.h"
#include "IvyFramework/IvyDataTools/interface/IvyStreamHelpers.hh"


using namespace std;
using namespace IvyStreamHelpers;


namespace BtagHelpers{
  std::unordered_map<BtagWPType, float> btagwptype_btagwp_map;
}

// All of the WPs comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation.
void BtagHelpers::configureBtagWPs(){
  if (!SampleHelpers::runConfigure){
    IVYerr << "BtagHelpers::configureBtagWPs: SampleHelpers::configure needs to be called first." << endl;
    assert(0);
  }

  auto const& dp = SampleHelpers::getDataPeriod();
  auto const& dy = SampleHelpers::getDataYear();
  bool const isAPV2016Affected = SampleHelpers::isAPV2016Affected(dp);

  switch (dy){
  case 2016:
  {
    if (isAPV2016Affected){
      btagwptype_btagwp_map = std::unordered_map<BtagWPType, float>{
        { kDeepFlav_Loose, 0.0508 },
        { kDeepFlav_Medium, 0.2598 },
        { kDeepFlav_Tight, 0.6502 },
        { kDeepCSV_Loose, 0.2027 },
        { kDeepCSV_Medium, 0.6001 },
        { kDeepCSV_Tight, 0.8819 }
      };
    }
    else{
      btagwptype_btagwp_map = std::unordered_map<BtagWPType, float>{
        { kDeepFlav_Loose, 0.0480 },
        { kDeepFlav_Medium, 0.2489 },
        { kDeepFlav_Tight, 0.6377 },
        { kDeepCSV_Loose, 0.1918 },
        { kDeepCSV_Medium, 0.5847 },
        { kDeepCSV_Tight, 0.8767 }
      };
    }
    break;
  }
  case 2017:
  {
    btagwptype_btagwp_map = std::unordered_map<BtagWPType, float>{
      { kDeepFlav_Loose, 0.0532 },
      { kDeepFlav_Medium, 0.3040 },
      { kDeepFlav_Tight, 0.7476 },
      { kDeepCSV_Loose, 0.1355 },
      { kDeepCSV_Medium, 0.4506 },
      { kDeepCSV_Tight, 0.7738 }
    };
    break;
  }
  case 2018:
  case 2022:
  {
    btagwptype_btagwp_map = std::unordered_map<BtagWPType, float>{
      { kDeepFlav_Loose, 0.0490 },
      { kDeepFlav_Medium, 0.2783 },
      { kDeepFlav_Tight, 0.7100 },
      { kDeepCSV_Loose, 0.1208 },
      { kDeepCSV_Medium, 0.4168 },
      { kDeepCSV_Tight, 0.7665 }
    };
    if (dy==2022) IVYout << "BtagHelpers::configureBtagWPs: WARNING! Using the WPs for year 2018 in place of " << dy << "." << endl;
    break;
  }
  default:
    IVYerr << "BtagHelpers::configureBtagWPs: Data period " << dp << " is not defined." << endl;
    assert(0);
  }
}
float const& BtagHelpers::getBtagWP(BtagHelpers::BtagWPType type){
  auto it = btagwptype_btagwp_map.find(type);
  if (it == btagwptype_btagwp_map.cend()){
    IVYerr << "BtagHelpers::getBtagWP: WP for WP type  " << type << " is not found." << endl;
    assert(0);
  }
  return it->second;
}
std::vector<float> BtagHelpers::getBtagWPs(BtagHelpers::BtagWPType type){
  std::vector<float> res; res.reserve(3);
  switch (type){
  case kDeepFlav_Loose:
  case kDeepFlav_Medium:
  case kDeepFlav_Tight:
    res.push_back(getBtagWP(kDeepFlav_Loose));
    res.push_back(getBtagWP(kDeepFlav_Medium));
    res.push_back(getBtagWP(kDeepFlav_Tight));
    break;
  case kDeepCSV_Loose:
  case kDeepCSV_Medium:
  case kDeepCSV_Tight:
    res.push_back(getBtagWP(kDeepCSV_Loose));
    res.push_back(getBtagWP(kDeepCSV_Medium));
    res.push_back(getBtagWP(kDeepCSV_Tight));
    break;
  default:
    IVYerr << "BtagHelpers::getBtagWPs: No implementation for the b-tagging WP type " << type << ". Aborting..." << endl;
    assert(0);
  }
  return res;
}

TString BtagHelpers::getBtagSFFileName(BtagWPType type){
  auto const& dp = SampleHelpers::getDataPeriod();
  auto const& dy = SampleHelpers::getDataYear();
  bool const isAPV2016Affected = SampleHelpers::isAPV2016Affected(dp);
  TString const strdp = TString(Form("%i", dy)) + (dy!=2016 ? "" : (isAPV2016Affected ? "_APV" : "_NonAPV"));

  bool valid_empty = false;

  TString res;
  switch (dy){
  case 2016:
  {
    switch (type){
    case kDeepFlav_Loose:
    case kDeepFlav_Medium:
    case kDeepFlav_Tight:
      // FIXME: Per threads like https://cms-talk.web.cern.ch/t/a-possible-issue-in-the-b-tagging-correction-for-2016postvfp-with-wp-loose/18837/2,
      // one should apply the APV SFs on both data sets for now.
      res = (isAPV2016Affected ? "wp_deepJet_106XUL16preVFP_v2.csv" : "wp_deepJet_106XUL16preVFP_v2.csv");
      break;
    default:
      break;
    }
    break;
  }
  case 2017:
  {
    switch (type){
    case kDeepFlav_Loose:
    case kDeepFlav_Medium:
    case kDeepFlav_Tight:
      res = "wp_deepJet_106XUL17_v3.csv";
      break;
    default:
      break;
    }
    break;
  }
  case 2018:
  {
    switch (type){
    case kDeepFlav_Loose:
    case kDeepFlav_Medium:
    case kDeepFlav_Tight:
      res = "wp_deepJet_106XUL18_v2.csv";
      break;
    default:
      break;
    }
    break;
  }
  case 2022:
  {
    valid_empty = true;
    break;
  }
  default:
    break;
  }
  if (res == ""){
    (!valid_empty ? IVYerr : IVYout)
      << "BtagHelpers::getBtagSFFileName: "
      << (valid_empty ? "WARNING! " : "")
      << "WP " << type << " is not implemented for period " << strdp << "."
      << endl;
    assert(valid_empty);
  }
  else{
    res = ANALYSISPKGDATAPATH + Form("external/BTagging/ScaleFactors/") + res;
    HostHelpers::ExpandEnvironmentVariables(res);
    if (!HostHelpers::FileReadable(res.Data())){
      IVYerr << "BtagHelpers::getBtagSFFileName: File " << res << " is not readable." << endl;
      assert(0);
    }
  }

  return res;
}
TString BtagHelpers::getBtagEffFileName(){
  auto const& dp = SampleHelpers::getDataPeriod();
  auto const& dy = SampleHelpers::getDataYear();
  bool const isAPV2016Affected = SampleHelpers::isAPV2016Affected(dp);
  TString const strdp = TString(Form("%i", dy)) + (dy!=2016 ? "" : (isAPV2016Affected ? "_APV" : "_NonAPV"));

  if (dy>=2022){
    IVYout << "BtagHelpers::getBtagEffFileName: WARNING! b-tagging efficiency is not known yet for period " << strdp << "." << endl;
    return "";
  }

  TString res = ANALYSISPKGDATAPATH + Form("ScaleFactors/bTagging/%s/Final_bTag_Efficiencies_AllMC.root", strdp.Data());
  HostHelpers::ExpandEnvironmentVariables(res);
  if (!HostHelpers::FileReadable(res.Data())){
    IVYerr << "BtagHelpers::getBtagEffFileName: File " << res << " is not readable." << endl;
    assert(0);
  }
  return res;
}
TString BtagHelpers::getBtagEffHistName(BtagWPType type, const char* jet_type){
  switch (type){
  case kDeepFlav_Loose:
    return Form("DeepFlavor_LooseJets_%s", jet_type);
  case kDeepFlav_Medium:
    return Form("DeepFlavor_MediumJets_%s", jet_type);
  case kDeepFlav_Tight:
    return Form("DeepFlavor_TightJets_%s", jet_type);
  default:
    IVYerr << "BtagHelpers::getBtagEffHistName: WP " << type << " is not implemented." << endl;
    assert(0);
    break;
  }
  return "";
}
