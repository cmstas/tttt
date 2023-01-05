#ifndef JECSCALEFACTORHANDLER_H
#define JECSCALEFACTORHANDLER_H

#include <vector>
#include <utility>
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "ScaleFactorHandlerBase.h"
#include "JESRHelpers.h"
#include "ParticleObject.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectionUncertainty.h"


class JECScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  JESRHelpers::JetType const type;

protected:
  FactorizedJetCorrector* corrector_data;
  FactorizedJetCorrector* corrector_MC;
  JetCorrectionUncertainty* uncertaintyEstimator_MC;

  FactorizedJetCorrector* makeCorrector(std::vector<TString> const& fnames);
  JetCorrectionUncertainty* makeUncertaintyEstimator(TString const& fname);

public:
  JECScaleFactorHandler(JESRHelpers::JetType type_);
  ~JECScaleFactorHandler();

  bool setup();
  void reset();

  void applyJEC(ParticleObject* obj, float const& rho, bool isMC, bool doRecompute);

};


#endif
