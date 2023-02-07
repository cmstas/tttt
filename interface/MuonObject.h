#ifndef MUONOBJECT_H
#define MUONOBJECT_H

#include "ParticleObject.h"

#define MUON_EXTRA_OPTIONAL_VARIABLES_COMMON \
MUON_VARIABLE(bool, highPurity, 0) \
MUON_VARIABLE(bool, inTimeMuon, 0) \
MUON_VARIABLE(bool, isPFcand, 0) \
MUON_VARIABLE(bool, isTracker, 0) \
MUON_VARIABLE(bool, tightId, 0) \
MUON_VARIABLE(bool, triggerIdLoose, 0) \
MUON_VARIABLE(int, nStations, 0) \
MUON_VARIABLE(float, tkRelIso, 0) \
MUON_VARIABLE(float, pfRelIso03_chg, 0)

#define MUON_EXTRA_OPTIONAL_VARIABLES_RUN2

#define MUON_EXTRA_OPTIONAL_VARIABLES_RUN3

#define MUON_EXTRA_OPTIONAL_VARIABLES \
MUON_EXTRA_OPTIONAL_VARIABLES_COMMON \
MUON_EXTRA_OPTIONAL_VARIABLES_RUN2 \
MUON_EXTRA_OPTIONAL_VARIABLES_RUN3


#define MUON_EXTRA_MANDATORY_VARIABLES_COMMON \
MUON_VARIABLE(bool, isGlobal, 0) \
MUON_VARIABLE(bool, looseId, 0) \
MUON_VARIABLE(bool, mediumId, 0) \
MUON_VARIABLE(unsigned char, jetNDauCharged, 0) \
MUON_VARIABLE(int, jetIdx, -1) \
MUON_VARIABLE(float, segmentComp, 0) \
MUON_VARIABLE(float, miniPFRelIso_all, 0) \
MUON_VARIABLE(float, miniPFRelIso_chg, 0) \
MUON_VARIABLE(float, pfRelIso03_all, 0) \
MUON_VARIABLE(float, dxy, 0) \
MUON_VARIABLE(float, dz, 0) \
MUON_VARIABLE(float, sip3d, 0) \
MUON_VARIABLE(float, jetPtRelv2, 0) \
MUON_VARIABLE(float, jetRelIso, 0) \
MUON_VARIABLE(float, ptErr, 0)

#define MUON_EXTRA_MANDATORY_VARIABLES_RUN2

#define MUON_EXTRA_MANDATORY_VARIABLES_RUN3 \
MUON_VARIABLE(bool, isStandalone, 0)

#define MUON_EXTRA_MANDATORY_VARIABLES \
MUON_EXTRA_MANDATORY_VARIABLES_COMMON \
MUON_EXTRA_MANDATORY_VARIABLES_RUN2 \
MUON_EXTRA_MANDATORY_VARIABLES_RUN3


#define MUON_EXTRA_VARIABLES \
MUON_EXTRA_MANDATORY_VARIABLES \
MUON_EXTRA_OPTIONAL_VARIABLES


class MuonVariables{
public:
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  MUON_EXTRA_VARIABLES;
#undef MUON_VARIABLE

  MuonVariables();
  MuonVariables(MuonVariables const& other);
  MuonVariables& operator=(const MuonVariables& other);

  void swap(MuonVariables& other);
};

class MuonObject : public ParticleObject{
public:
  MuonVariables extras;

  MuonObject();
  MuonObject(int id_);
  MuonObject(int id_, LorentzVector_t const& mom_);
  MuonObject(const MuonObject& other);
  MuonObject& operator=(const MuonObject& other);
  ~MuonObject();

  void swap(MuonObject& other);

  LorentzVector_t::Scalar ptrel() const;
  LorentzVector_t::Scalar ptratio() const;

};

#endif
