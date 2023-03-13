/*
Originally from https://raw.githubusercontent.com/usarica/CMS3NtupleTools/combined/AnalysisTree/src/BtagScaleFactorHandler.cc (access: 23/01/01)
with modified calls for compatibility.
*/

#include <cassert>
#include "IvyFramework/IvyDataTools/interface/SampleHelpersCore.h"
#include "SamplesCore.h"
#include "BtagScaleFactorHandler.h"
#include "AK4JetSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "IvyFramework/IvyDataTools/interface/IvyStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace BtagHelpers;
using namespace IvyStreamHelpers;


BtagScaleFactorHandler::BtagScaleFactorHandler() : ScaleFactorHandlerBase()
{
  setup();
}

BtagScaleFactorHandler::~BtagScaleFactorHandler(){ this->reset(); }

void BtagScaleFactorHandler::evalEfficiencyFromHistogram(float& theSF, float const& pt, float const& eta, ExtendedHistogram_2D_f const& hist, bool etaOnY, bool useAbsEta) const{
  auto const* hh = hist.getHistogram();
  if (!hh){
    IVYerr << "BtagScaleFactorHandler::evalScaleFactorFromHistogram: Histogram is null." << endl;
    return;
  }

  float const abs_eta = std::abs(eta);
  float const& eta_to_use = (!useAbsEta ? eta : abs_eta);

  int ix, iy;
  int nbinsx = hh->GetNbinsX();
  int nbinsy = hh->GetNbinsY();
  if (!etaOnY){
    ix = hh->GetXaxis()->FindBin(eta_to_use);
    iy = hh->GetYaxis()->FindBin(pt);
  }
  else{
    ix = hh->GetXaxis()->FindBin(pt);
    iy = hh->GetYaxis()->FindBin(eta_to_use);
  }

  if (ix==0) ix=1;
  else if (ix>nbinsx) ix=nbinsx;
  if (iy==0) iy=1;
  else if (iy>nbinsy) iy=nbinsy;

  theSF = hh->GetBinContent(ix, iy);
}

bool BtagScaleFactorHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  this->reset();

  if (verbosity>=MiscUtils::INFO) IVYout << "BtagScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << SampleHelpers::getDataYear() << endl;

  TDirectory* curdir = gDirectory;
  TDirectory* uppermostdir = SampleHelpers::rootTDirectory;

  TString sf_fname_deepFlavor = BtagHelpers::getBtagSFFileName(kDeepFlav_Loose);
  WP_calib_map[kDeepFlav_Loose] = new BTagCalibration("DeepFlavor", sf_fname_deepFlavor.Data());

  for (int iwp=0; iwp<(int) nBtagWPTypes; iwp++){
    BtagWPType type = (BtagWPType) iwp;
    BTagEntry::OperatingPoint opPoint;
    BTagCalibration* calibration=nullptr;

    if (iwp>=3) break; // FIXME: Implement DeepCSV SFs once efficiencies become available fully.
    switch (type){
    case kDeepFlav_Loose:
      opPoint = BTagEntry::OP_LOOSE;
      calibration = WP_calib_map[kDeepFlav_Loose];
      break;
    case kDeepFlav_Medium:
      opPoint = BTagEntry::OP_MEDIUM;
      calibration = WP_calib_map[kDeepFlav_Loose];
      break;
    case kDeepFlav_Tight:
      opPoint = BTagEntry::OP_TIGHT;
      calibration = WP_calib_map[kDeepFlav_Loose];
      break;
    case kDeepCSV_Loose:
      opPoint = BTagEntry::OP_LOOSE;
      calibration = WP_calib_map[kDeepCSV_Loose];
      break;
    case kDeepCSV_Medium:
      opPoint = BTagEntry::OP_MEDIUM;
      calibration = WP_calib_map[kDeepCSV_Loose];
      break;
    case kDeepCSV_Tight:
      opPoint = BTagEntry::OP_TIGHT;
      calibration = WP_calib_map[kDeepCSV_Loose];
      break;
    default:
      IVYerr << "BtagScaleFactorHandler::setup: No operating point implementation for b tag WP " << type << ". Aborting..." << endl;
      assert(0);
      break;
    }

    WP_calibreader_map_nominal[type] = new BTagCalibrationReader(opPoint, "central");
    loadBTagCalibrations(WP_calibreader_map_nominal[type], calibration);

    WP_calibreader_map_dn[type] = new BTagCalibrationReader(opPoint, "down");
    loadBTagCalibrations(WP_calibreader_map_dn[type], calibration);

    WP_calibreader_map_up[type] = new BTagCalibrationReader(opPoint, "up");
    loadBTagCalibrations(WP_calibreader_map_up[type], calibration);
  }

  std::vector<SystematicsHelpers::SystematicVariationTypes> const allowedSysts{
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    ePUJetIdEffDn, ePUJetIdEffUp
  };
  std::vector< std::pair<BTagEntry::JetFlavor, TString> > const flavpairs{ { BTagEntry::FLAV_B, "b" },{ BTagEntry::FLAV_C, "c" },{ BTagEntry::FLAV_UDSG, "udsg" } };
  // Feature to be added if needed: The implementation could be generalized for PU jet cats. For now, keep it inclusive.
  std::vector<TString> const strpujetidcats{ "Inclusive" /*"T", "MnT", "LnM", "F"*/ };
  TFile* finput_eff = TFile::Open(BtagHelpers::getBtagEffFileName(), "read"); uppermostdir->cd();
  if (verbosity>=MiscUtils::INFO) IVYout << "BtagScaleFactorHandler::setup: Reading " << finput_eff->GetName() << " to acquire efficiency histograms..." << endl;
  {
    ExtendedHistogram_2D_f empty_hist; empty_hist.reset();
    TString hname;
    for (auto const& syst:allowedSysts){
      TString systname = SystematicsHelpers::getSystName(syst).data();
      syst_flav_pujetid_WP_mceffhist_map[syst] = std::unordered_map< BTagEntry::JetFlavor, std::vector<std::vector<ExtendedHistogram_2D_f>> >();
      for (auto const& flavpair:flavpairs){
        BTagEntry::JetFlavor jflav = flavpair.first;
        TString const& strflav = flavpair.second;
        syst_flav_pujetid_WP_mceffhist_map[syst][jflav] = std::vector<std::vector<ExtendedHistogram_2D_f>>(strpujetidcats.size(), std::vector<ExtendedHistogram_2D_f>(nBtagWPTypes, empty_hist));
        for (unsigned short ipujetidwp=0; ipujetidwp<strpujetidcats.size(); ipujetidwp++){
          TString const& strpujetidcat = strpujetidcats.at(ipujetidwp);
          if (strpujetidcat=="Inclusive" && (syst==ePUJetIdEffDn || syst==ePUJetIdEffUp)) continue;
          for (int iwp=0; iwp<(int) nBtagWPTypes; iwp++){
            if (iwp>=3) break; // FIXME: Implement DeepCSV SFs once efficiencies become available fully.
            BtagWPType wptype = (BtagWPType) iwp;
            hname = BtagHelpers::getBtagEffHistName(wptype, strflav.Data()); hname = hname + "_PUJetId_" + strpujetidcat + "_" + systname;
            if (verbosity>=MiscUtils::DEBUG) IVYout << "\t- Extracting MC efficiency histogram " << hname << "..." << endl;
            bool tmpres = getHistogram<TH2D, ExtendedHistogram_2D_f>(syst_flav_pujetid_WP_mceffhist_map[syst][jflav].at(ipujetidwp).at(iwp), finput_eff, hname);
            if (!tmpres && verbosity>=MiscUtils::DEBUG) IVYerr << "\t\t- FAILED!" << endl;
            res &= tmpres;
          }
        }
      }
    }
  }
  ScaleFactorHandlerBase::closeFile(finput_eff); curdir->cd();

  return res;
}
void BtagScaleFactorHandler::loadBTagCalibrations(BTagCalibrationReader* const& reader, BTagCalibration* const& calibration){
  reader->load(*calibration, BTagEntry::FLAV_B, "comb");
  reader->load(*calibration, BTagEntry::FLAV_C, "comb");
  reader->load(*calibration, BTagEntry::FLAV_UDSG, "incl");
}

void BtagScaleFactorHandler::reset(){
  for (auto it:WP_calibreader_map_up) delete it.second;
  for (auto it:WP_calibreader_map_dn) delete it.second;
  for (auto it:WP_calibreader_map_nominal) delete it.second;
  for (auto it:WP_calib_map) delete it.second;

  WP_calibreader_map_up.clear();
  WP_calibreader_map_dn.clear();
  WP_calibreader_map_nominal.clear();
  WP_calib_map.clear();
  syst_flav_pujetid_WP_mceffhist_map.clear();
}

float BtagScaleFactorHandler::getSFFromBTagCalibrationReader(
  BTagCalibrationReader const* calibReader, BTagCalibrationReader const* calibReader_Nominal,
  SystematicsHelpers::SystematicVariationTypes const& syst, BTagEntry::JetFlavor const& flav, float const& pt, float const& eta
) const{
  float myPt = pt;
  float const MaxJetEta = (SampleHelpers::getDataYear()<=2016 ? 2.4 : 2.5);
  if (std::abs(eta) > MaxJetEta) return 1; // Do not apply SF for jets with eta higher than the threshold

  std::pair<float, float> pt_min_max = calibReader->min_max_pt(flav, eta);
  if (pt_min_max.second<0.) return 1;
  bool DoubleUncertainty = false;
  if (myPt<pt_min_max.first){
    myPt = pt_min_max.first+1e-5;
    DoubleUncertainty = true;
  }
  else if (myPt>pt_min_max.second){
    myPt = pt_min_max.second-1e-5;
    DoubleUncertainty = true;
  }

  float SF = calibReader->eval(flav, eta, myPt);
  if (DoubleUncertainty && (syst==SystematicsHelpers::eBTagSFDn || syst==SystematicsHelpers::eBTagSFUp)){
    float SFcentral = calibReader_Nominal->eval(flav, eta, myPt);
    SF = 2.f*(SF - SFcentral) + SFcentral;
  }

  return SF;
}

void BtagScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, short const& pujetidcat, BTagEntry::JetFlavor const& flav, unsigned short const& btagcat, float& val, float* effval) const{
  using namespace SystematicsHelpers;
  auto const& btagWPType = AK4JetSelectionHelpers::btagger_type;

  val = 1;

  if (this->verbosity>=MiscUtils::DEBUG) IVYout
    << "BtagScaleFactorHandler::getSFAndEff: Calling for jet (pt, eta, PU jet id cat., flav, b-tag cat.) = ("
    << pt << ", " << eta << ", " << pujetidcat << ", " << flav << ", " << btagcat << "):"
    << endl;

  std::unordered_map<BtagHelpers::BtagWPType, BTagCalibrationReader*> const* WP_calibreader_map = nullptr;
  switch (syst){
  case SystematicsHelpers::eBTagSFDn:
    WP_calibreader_map = &WP_calibreader_map_dn;
    break;
  case SystematicsHelpers::eBTagSFUp:
    WP_calibreader_map = &WP_calibreader_map_up;
    break;
  default:
    WP_calibreader_map = &WP_calibreader_map_nominal;
    break;
  }

  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedJetSysts{
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp
  };
  if (pujetidcat>=0){
    allowedJetSysts.push_back(ePUJetIdEffDn);
    allowedJetSysts.push_back(ePUJetIdEffUp);
  }
  SystematicsHelpers::SystematicVariationTypes const jetsyst = (HelperFunctions::checkListVariable(allowedJetSysts, syst) ? syst : SystematicsHelpers::sNominal);

  std::vector<BTagCalibrationReader const*> calibReaders;
  std::vector<BTagCalibrationReader const*> calibReaders_Nominal;
  // Feature to be added if needed: The implementation could be generalized for PU jet cats. For now, keep it inclusive.
  // FIXME: If this feature is added, BE CAREFUL WITH THE OFFSET OF +1 ON THE LINE BELOW!
  std::vector<ExtendedHistogram_2D_f> const& effhists = syst_flav_pujetid_WP_mceffhist_map.find(jetsyst)->second.find(flav)->second.at(pujetidcat+1);
  unsigned short idx_offset_effmc = 0;
  switch (btagWPType){
  case kDeepFlav_Loose:
  case kDeepFlav_Medium:
  case kDeepFlav_Tight:
  {
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepFlav_Loose)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepFlav_Medium)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepFlav_Tight)->second);

    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepFlav_Loose)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepFlav_Medium)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepFlav_Tight)->second);

    idx_offset_effmc = (unsigned short) kDeepFlav_Loose;
    break;
  }
  case kDeepCSV_Loose:
  case kDeepCSV_Medium:
  case kDeepCSV_Tight:
  {
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepCSV_Loose)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepCSV_Medium)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepCSV_Tight)->second);

    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepCSV_Loose)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepCSV_Medium)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepCSV_Tight)->second);

    idx_offset_effmc = (unsigned short) kDeepCSV_Loose;
    break;
  }
  default:
    IVYerr << "BtagScaleFactorHandler::getSFAndEff: b tag calibration readers are not implemented." << endl;
    assert(0);
    break;
  }

  std::vector<float> effs_unscaled(calibReaders.size(), 1);
  std::vector<float> effs_scaled(calibReaders.size(), 1);
  std::vector<float> SFs(calibReaders.size(), 1);
  for (unsigned int i=0; i<calibReaders.size(); i++){
    SFs.at(i) = getSFFromBTagCalibrationReader(
      calibReaders.at(i), calibReaders_Nominal.at(i),
      syst, flav, pt, eta
    );
    if (i>0) SFs.at(i) /= SFs.at(i-1);
    evalEfficiencyFromHistogram(effs_unscaled.at(i), pt, eta, effhists.at(i+idx_offset_effmc), true, false);
    effs_scaled.at(i) = std::max(0.f, std::min(1.f, effs_unscaled.at(i) * SFs.at(i)));
  }

  float eff_unscaled = 1;
  float eff_scaled = 1;
  for (unsigned short iwp=0; iwp<calibReaders.size(); iwp++){
    float tmp_eff_unscaled = 1;
    float tmp_eff_scaled = 1;
    bool doContinue = true;
    if (btagcat>iwp){
      tmp_eff_unscaled *= effs_unscaled.at(iwp);
      tmp_eff_scaled *= effs_scaled.at(iwp);
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Passed";
    }
    else{
      tmp_eff_unscaled *= 1.f - effs_unscaled.at(iwp);
      tmp_eff_scaled *= 1.f - effs_scaled.at(iwp);
      doContinue = false;
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Failed";
    }
    if (this->verbosity>=MiscUtils::DEBUG) IVYout << " WP " << iwp << " with unscaled, scaled effs = " << tmp_eff_unscaled << ", " << tmp_eff_scaled << endl;
    eff_unscaled *= tmp_eff_unscaled;
    eff_scaled *= tmp_eff_scaled;
    if (!doContinue) break;
  }

  val = (eff_unscaled>0.f ? eff_scaled / eff_unscaled : 0.f);
  if (effval) *effval = (val>0.f ? eff_scaled : 0.f);

  if (this->verbosity>=MiscUtils::DEBUG){
    IVYout << "\t- Unscaled effs: " << effs_unscaled << endl;
    IVYout << "\t- Scaled effs: " << effs_scaled << endl;
    IVYout << "\t- SFs: " << SFs << endl;
    IVYout << "\t- Final SF: " << val << endl;
    IVYout << "\t- Final eff: " << eff_scaled << endl;
  }
}
void BtagScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, AK4JetObject const* obj, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 1;

  if (!obj) return;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTaggable)) return;

  // Feature to be added if needed: The implementation could be generalized for PU jet cats. For now, keep it inclusive.
  constexpr short pujetidcat = -1;
  /*
  short pujetidcat = 0;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kTightPUJetId)) pujetidcat++;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kMediumPUJetId)) pujetidcat++;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kLoosePUJetId)) pujetidcat++;
  */

  unsigned short btagcat = 0;
  if (obj->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Loose)) btagcat++;
  if (obj->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Medium)) btagcat++;
  if (obj->testSelectionBit(AK4JetSelectionHelpers::kPreselectionTight_BTagged_Tight)) btagcat++;

  getSFAndEff(syst, obj->pt(), obj->eta(), pujetidcat, obj->getBTagJetFlavor(), btagcat, val, effval);
}
