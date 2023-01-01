#ifndef BTAGHELPERS_H
#define BTAGHELPERS_H

#include <vector>


namespace BtagHelpers{
  enum BtagWPType{
    /*
    kDeepCSV_Loose=0,
    kDeepCSV_Medium,
    kDeepCSV_Tight,

    kDeepFlav_Loose,
    */
    kDeepFlav_Loose=0,
    kDeepFlav_Medium,
    kDeepFlav_Tight,

    nBtagWPTypes
  };

  void configureBtagWPs();

  float const& getBtagWP(BtagWPType type);
  std::vector<float> getBtagWPs(BtagWPType type);

  TString getBtagSFFileName(BtagWPType type);
  TString getBtagEffFileName();
  TString getBtagEffHistName(BtagWPType type, const char* jet_type);

}


#endif
