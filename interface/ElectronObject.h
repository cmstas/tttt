#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "ParticleObject.h"


#define ELECTRON_EXTRA_OPTIONAL_VARIABLES_COMMON \
ELECTRON_VARIABLE(bool, isPFcand, 0) \
ELECTRON_VARIABLE(unsigned char, seedGain, 0) \
ELECTRON_VARIABLE(int, cutBased, 0) \
ELECTRON_VARIABLE(float, dr03EcalRecHitSumEt, 0) \
ELECTRON_VARIABLE(float, dr03HcalDepth1TowerSumEt, 0) \
ELECTRON_VARIABLE(float, dr03TkSumPt, 0) \
ELECTRON_VARIABLE(float, pfRelIso03_chg, 0) \
ELECTRON_VARIABLE(float, ip3d, 0) \
ELECTRON_VARIABLE(float, r9, 0) \
ELECTRON_VARIABLE(float, scEtOverPt, 0)

#define ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN2

#define ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN3

#define ELECTRON_EXTRA_OPTIONAL_VARIABLES \
ELECTRON_EXTRA_OPTIONAL_VARIABLES_COMMON \
ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN2 \
ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN3


#define ELECTRON_EXTRA_MANDATORY_VARIABLES_COMMON \
ELECTRON_VARIABLE(bool, convVeto, 0) \
ELECTRON_VARIABLE(unsigned char, lostHits, 0) \
ELECTRON_VARIABLE(unsigned char, jetNDauCharged, 0) \
ELECTRON_VARIABLE(int, jetIdx, -1) \
ELECTRON_VARIABLE(int, tightCharge, 0) \
ELECTRON_VARIABLE(float, hoe, 0) \
ELECTRON_VARIABLE(float, eInvMinusPInv, 0) \
ELECTRON_VARIABLE(float, sieie, 0) \
ELECTRON_VARIABLE(float, dxy, 0) \
ELECTRON_VARIABLE(float, dz, 0) \
ELECTRON_VARIABLE(float, sip3d, 0) \
ELECTRON_VARIABLE(float, miniPFRelIso_all, 0) \
ELECTRON_VARIABLE(float, miniPFRelIso_chg, 0) \
ELECTRON_VARIABLE(float, pfRelIso03_all, 0) \
ELECTRON_VARIABLE(float, jetPtRelv2, 0) \
ELECTRON_VARIABLE(float, jetRelIso, 0) \
ELECTRON_VARIABLE(float, deltaEtaSC, 0)

#define ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN2 \
ELECTRON_VARIABLE(bool, mvaFall17V2noIso_WPL, 0) \
ELECTRON_VARIABLE(float, mvaFall17V2noIso, 0) \
ELECTRON_VARIABLE(float, mvaFall17V2Iso, 0)

#define ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN3 \
ELECTRON_VARIABLE(bool, mvaNoIso_WPL, 0) \
ELECTRON_VARIABLE(float, mvaNoIso, 0) \
ELECTRON_VARIABLE(float, mvaIso, 0)

#define ELECTRON_EXTRA_MANDATORY_VARIABLES \
ELECTRON_EXTRA_MANDATORY_VARIABLES_COMMON \
ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN2 \
ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN3


#define ELECTRON_EXTRA_VARIABLES \
ELECTRON_EXTRA_MANDATORY_VARIABLES \
ELECTRON_EXTRA_OPTIONAL_VARIABLES


class ElectronVariables{
public:
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  ELECTRON_EXTRA_VARIABLES;
#undef ELECTRON_VARIABLE

  ElectronVariables();
  ElectronVariables(ElectronVariables const& other);
  ElectronVariables& operator=(const ElectronVariables& other);

  void swap(ElectronVariables& other);
};

class ElectronObject : public ParticleObject{
public:
  ElectronVariables extras;

  ElectronObject();
  ElectronObject(int id_);
  ElectronObject(int id_, LorentzVector_t const& mom_);
  ElectronObject(const ElectronObject& other);
  ElectronObject& operator=(const ElectronObject& other);
  ~ElectronObject();

  void swap(ElectronObject& other);

  float etaSC() const{ return this->extras.deltaEtaSC + this->eta(); }

  LorentzVector_t::Scalar ptrel() const;
  LorentzVector_t::Scalar ptratio() const;

};

#endif
