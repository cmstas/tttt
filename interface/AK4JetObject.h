#ifndef AK4JETOBJECT_H
#define AK4JETOBJECT_H

#include "ParticleObject.h"
#include "BtagHelpers.h"
#include "BTagCalibrationStandalone.h"
#include "SystematicVariations.h"


#define AK4JET_GENINFO_VARIABLES \
AK4JET_VARIABLE(int, hadronFlavour, 0) \
AK4JET_VARIABLE(int, partonFlavour, 0)

#define AK4JET_COMMON_VARIABLES \
AK4JET_VARIABLE(int, jetId, 0) \
AK4JET_VARIABLE(int, electronIdx1, -1) \
AK4JET_VARIABLE(int, electronIdx2, -1) \
AK4JET_VARIABLE(int, muonIdx1, -1) \
AK4JET_VARIABLE(int, muonIdx2, -1) \
AK4JET_VARIABLE(float, area, 0) \
AK4JET_VARIABLE(float, muonSubtrFactor, 0) \
AK4JET_VARIABLE(float, muEF, 0) \
AK4JET_VARIABLE(float, rawFactor, 0) \
AK4JET_VARIABLE(float, btagDeepFlavB, 0) \
AK4JET_VARIABLE(float, btagDeepB, 0)

// These variables are a subset of AK4JET_COMMON_VARIABLES, so there is no need to define them twice in the extras.
// Nevertheless, define them through a different command set (i.e., AK4JET_LOWPT_VARIABLE instead of AK4JET_VARIABLE).
#define AK4JET_LOWPT_EXTRA_INPUT_VARIABLES \
AK4JET_LOWPT_VARIABLE(float, area, 0) \
AK4JET_LOWPT_VARIABLE(float, muonSubtrFactor, 0)

#define AK4JET_EXTRA_INPUT_VARIABLES \
AK4JET_GENINFO_VARIABLES \
AK4JET_COMMON_VARIABLES

#define AK4JET_EXTRA_VARIABLES \
AK4JET_EXTRA_INPUT_VARIABLES \
AK4JET_VARIABLE(bool, isLowPtJet, false) \
AK4JET_VARIABLE(float, JECNominal, 1) \
AK4JET_VARIABLE(float, JECL1Nominal, 1) \
AK4JET_VARIABLE(float, relJECUnc, 0) \
AK4JET_VARIABLE(float, JERNominal, 1) \
AK4JET_VARIABLE(float, JERDn, 1) \
AK4JET_VARIABLE(float, JERUp, 1)


class AK4JetVariables{
public:
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  AK4JET_EXTRA_VARIABLES;
#undef AK4JET_VARIABLE

  AK4JetVariables();
  AK4JetVariables(AK4JetVariables const& other);
  AK4JetVariables& operator=(const AK4JetVariables& other);

  void swap(AK4JetVariables& other);
};


class AK4JetObject : public ParticleObject{
protected:
  SystematicsHelpers::SystematicVariationTypes currentSyst;
  float currentJEC;
  float currentJER;
  float currentSystScale;
  LorentzVector_t mom_original;
  LorentzVector_t mom_mucands;

public:
  constexpr static float ConeRadiusConstant = 0.4;

  AK4JetVariables extras;

  AK4JetObject();
  AK4JetObject(LorentzVector_t const& mom_);
  AK4JetObject(const AK4JetObject& other);
  AK4JetObject& operator=(const AK4JetObject& other);
  ~AK4JetObject();

  void swap(AK4JetObject& other);

  float const& getJECValue() const{ return this->currentJEC; }
  float const& getJERValue() const{ return this->currentJER; }
  SystematicsHelpers::SystematicVariationTypes const& getCurrentSyst() const{ return this->currentSyst; }

  LorentzVector_t uncorrected_p4() const{ return this->mom_original; }
  void reset_uncorrected_p4(LorentzVector_t const& new_p4){ this->mom_original = new_p4; }

  LorentzVector_t const& p4_mucands() const{ return this->mom_mucands; }
  void reset_p4_mucands(LorentzVector_t const& new_p4){ this->mom_mucands = new_p4; }

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst);

  BTagEntry::JetFlavor getBTagJetFlavor() const;
  float getBtagValue(BtagHelpers::BtagWPType type = BtagHelpers::nBtagWPTypes) const;

};

#endif
