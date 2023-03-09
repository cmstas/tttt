#ifndef AK4JETSELECTIONHELPERS_H
#define AK4JETSELECTIONHELPERS_H

#include <vector>
#include "AK4JetObject.h"
#include "BtagHelpers.h"


namespace AK4JetSelectionHelpers{
  enum SelectionBits{
    kKinOnly,
    kKinOnly_BTag,

    kJetIdOnly,

    kBTagged_Loose,
    kBTagged_Medium,
    kBTagged_Tight,

    kPreselectionTight,

    kPreselectionTight_BTaggable,
    kPreselectionTight_BTagged_Loose,
    kPreselectionTight_BTagged_Medium,
    kPreselectionTight_BTagged_Tight,

    nSelectionBits
  };

  // Kinematic pT thresholds
  constexpr float ptThr_btag_infimum = 20.;
  constexpr float ptThr = 25.;
  constexpr float ptThr_HEMVeto = 30.;

  // Kinematic eta thresholds
  constexpr float etaThr_common = 4.7;
  constexpr float etaThr_btag_Phase0Tracker = 2.4;
  constexpr float etaThr_btag_Phase1Tracker = 2.5;

  constexpr float deltaRThr = 0.4;

  constexpr int jetIdBitPos = 1;

  void setSelectionBits(AK4JetObject& part);

  // Default is BtagHelpers::kDeepFlav_Loose.
  extern BtagHelpers::BtagWPType btagger_type;
  void setBTaggerType(BtagHelpers::BtagWPType const& type);

}


#endif
