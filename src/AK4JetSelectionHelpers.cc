#include <cassert>
#include <cmath>
#include "AK4JetSelectionHelpers.h"
#include "IvyFramework/IvyDataTools/interface/HelperFunctions.h"
#include "IvyFramework/IvyDataTools/interface/IvyStreamHelpers.hh"
#include "SamplesCore.h"


namespace AK4JetSelectionHelpers{
  BtagHelpers::BtagWPType btagger_type = BtagHelpers::kDeepFlav_Loose;

  bool testKin(AK4JetObject const& part);
  bool testKin_BTag(AK4JetObject const& part);

  bool testJetId(AK4JetObject const& part);
  bool testId(AK4JetObject const& part);

  bool testPreselectionTight(AK4JetObject const& part);
  bool testPreselectionTight_BTaggable(AK4JetObject const& part);

  void setBTagBits(AK4JetObject& part);
}


using namespace std;
using namespace IvyStreamHelpers;

void AK4JetSelectionHelpers::setBTaggerType(BtagHelpers::BtagWPType const& type){ AK4JetSelectionHelpers::btagger_type = type; }

bool AK4JetSelectionHelpers::testKin(AK4JetObject const& part){
  return (part.pt()>=ptThr && std::abs(part.eta())<etaThr_common);
}
bool AK4JetSelectionHelpers::testKin_BTag(AK4JetObject const& part){
  auto const& dy = SampleHelpers::getDataYear();

  float etaThr_btag = -1;
  if (dy<=2016) etaThr_btag = etaThr_btag_Phase0Tracker;
  else etaThr_btag = etaThr_btag_Phase1Tracker;

  return (part.pt()>=std::max(ptThr_btag_infimum, ptThr) && std::abs(part.eta())<etaThr_btag);
}

bool AK4JetSelectionHelpers::testJetId(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.jetId, jetIdBitPos);
}
bool AK4JetSelectionHelpers::testId(AK4JetObject const& part){
  return part.testSelectionBit(kJetIdOnly); // Could add PU jet ID here as well...
}

bool AK4JetSelectionHelpers::testPreselectionTight(AK4JetObject const& part){
  return (
    part.testSelectionBit(kKinOnly)
    &&
    testId(part)
    );
}
bool AK4JetSelectionHelpers::testPreselectionTight_BTaggable(AK4JetObject const& part){
  return (
    part.testSelectionBit(kKinOnly_BTag)
    &&
    part.testSelectionBit(kPreselectionTight)
    );
}

void AK4JetSelectionHelpers::setBTagBits(AK4JetObject& part){
  auto const btagWPs = BtagHelpers::getBtagWPs(btagger_type);
  assert(btagWPs.size()==3);

  auto const bscore = part.getBtagValue();
  part.setSelectionBit(kBTagged_Loose, bscore>btagWPs.front());
  part.setSelectionBit(kBTagged_Medium, bscore>btagWPs.at(1));
  part.setSelectionBit(kBTagged_Tight, bscore>btagWPs.back());

  bool const pass_BTag_preselection = part.testSelectionBit(kPreselectionTight_BTaggable);
  part.setSelectionBit(kPreselectionTight_BTagged_Loose, pass_BTag_preselection && part.testSelectionBit(kBTagged_Loose));
  part.setSelectionBit(kPreselectionTight_BTagged_Medium, pass_BTag_preselection && part.testSelectionBit(kBTagged_Medium));
  part.setSelectionBit(kPreselectionTight_BTagged_Tight, pass_BTag_preselection && part.testSelectionBit(kBTagged_Tight));
}

void AK4JetSelectionHelpers::setSelectionBits(AK4JetObject& part){
  static_assert(std::numeric_limits<ParticleObject::SelectionBitsType_t>::digits >= nSelectionBits);

  part.setSelectionBit(kKinOnly, testKin(part));
  part.setSelectionBit(kKinOnly_BTag, testKin_BTag(part));

  part.setSelectionBit(kJetIdOnly, testJetId(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
  part.setSelectionBit(kPreselectionTight_BTaggable, testPreselectionTight_BTaggable(part));

  setBTagBits(part);
}
