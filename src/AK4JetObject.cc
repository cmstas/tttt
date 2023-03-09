#include <cassert>
#include <algorithm>
#include <utility>
#include "AK4JetObject.h"
#include "MuonObject.h"
#include "BtagHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "IvyFramework/IvyDataTools/interface/IvyStreamHelpers.hh"


using namespace std;
using namespace IvyStreamHelpers;


AK4JetVariables::AK4JetVariables(){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  AK4JET_EXTRA_VARIABLES;
#undef AK4JET_VARIABLE
}
AK4JetVariables::AK4JetVariables(AK4JetVariables const& other){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  AK4JET_EXTRA_VARIABLES;
#undef AK4JET_VARIABLE
}
void AK4JetVariables::swap(AK4JetVariables& other){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  AK4JET_EXTRA_VARIABLES;
#undef AK4JET_VARIABLE
}
AK4JetVariables& AK4JetVariables::operator=(const AK4JetVariables& other){
  AK4JetVariables tmp(other);
  this->swap(tmp);
  return *this;
}

AK4JetObject::AK4JetObject() :
  ParticleObject(),
  mom_original(0, 0, 0, 0),
  mom_mucands(0, 0, 0, 0),
  currentSyst(SystematicsHelpers::sNominal),
  currentJEC(1),
  currentJER(1),
  currentSystScale(1),
  extras()
{}
AK4JetObject::AK4JetObject(LorentzVector_t const& momentum_) :
  ParticleObject(0, momentum_),
  mom_original(momentum_),
  mom_mucands(0, 0, 0, 0),
  currentSyst(SystematicsHelpers::sNominal),
  currentJEC(1),
  currentJER(1),
  currentSystScale(1),
  extras()
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  mom_original(other.mom_original),
  mom_mucands(other.mom_mucands),
  currentSyst(other.currentSyst),
  currentJEC(other.currentJEC),
  currentJER(other.currentJER),
  currentSystScale(other.currentSystScale),
  extras(other.extras)
{}
void AK4JetObject::swap(AK4JetObject& other){
  ParticleObject::swap(other);
  std::swap(mom_original, other.mom_original);
  std::swap(mom_mucands, other.mom_mucands);
  std::swap(currentSyst, other.currentSyst);
  std::swap(currentJEC, other.currentJEC);
  std::swap(currentJER, other.currentJER);
  std::swap(currentSystScale, other.currentSystScale);
  extras.swap(other.extras);
}
AK4JetObject& AK4JetObject::operator=(const AK4JetObject& other){
  AK4JetObject tmp(other);
  this->swap(tmp);
  return *this;
}
AK4JetObject::~AK4JetObject(){}

BTagEntry::JetFlavor AK4JetObject::getBTagJetFlavor() const{
  auto const& jetFlavor = extras.hadronFlavour;
  if (std::abs(jetFlavor)==5) return BTagEntry::FLAV_B;
  else if (std::abs(jetFlavor)==4) return BTagEntry::FLAV_C;
  else return BTagEntry::FLAV_UDSG;
}
float AK4JetObject::getBtagValue(BtagHelpers::BtagWPType type) const{
  if (type==BtagHelpers::nBtagWPTypes) type = AK4JetSelectionHelpers::btagger_type;
  switch (type){
  case BtagHelpers::kDeepFlav_Loose:
  case BtagHelpers::kDeepFlav_Medium:
  case BtagHelpers::kDeepFlav_Tight:
    return this->extras.btagDeepFlavB;
  case BtagHelpers::kDeepCSV_Loose:
  case BtagHelpers::kDeepCSV_Medium:
  case BtagHelpers::kDeepCSV_Tight:
    return this->extras.btagDeepB;
  default:
    IVYerr << "AK4JetObject::getBtagValue: b-tagger type " << type << " is not implemented. Aborting..." << endl;
    assert(0);
  }
  return -1;
}

void AK4JetObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  currentJEC = 1;
  currentJER = 1;
  momentum = mom_original;
  switch (syst){
  case eJECDn:
    currentJEC = extras.JECNominal*(1.f - extras.relJECUnc);
    currentJER = extras.JERNominal;
    break;
  case eJECUp:
    currentJEC = extras.JECNominal*(1.f + extras.relJECUnc);
    currentJER = extras.JERNominal;
    break;
  case eJERDn:
    currentJEC = extras.JECNominal;
    currentJER = extras.JERDn;
    break;
  case eJERUp:
    currentJEC = extras.JECNominal;
    currentJER = extras.JERUp;
    break;
  case sUncorrected:
    break;
  default:
    currentJEC = extras.JECNominal;
    currentJER = extras.JERNominal;
    break;
  }
  float scale = currentJEC * currentJER;
  // Test new pt
  float newpt = momentum.Pt() * scale;
  if (newpt<1e-5 && momentum.Pt()>0.f) scale = 1e-5 / momentum.Pt();
  momentum = momentum * scale;
  currentSystScale = scale;
  currentSyst = syst;
}
