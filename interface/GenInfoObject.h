#ifndef GENINFOOBJECT_H
#define GENINFOOBJECT_H

#include <unordered_map>
#include <string>
#include <vector>
#include "StdExtensions.h"
#include "SystematicVariations.h"


// NanoAOD variables
#define GENINFO_NANOAOD_SCALAR_VARIABLES \
GENINFO_NANOAOD_SCALAR_VARIABLE(int, genTtbarId, 0) \
GENINFO_NANOAOD_SCALAR_VARIABLE(float, genWeight, 1) \
GENINFO_NANOAOD_SCALAR_VARIABLE(float, LHEWeight_originalXWGTUP, 1)

#define GENINFO_NANOAOD_ARRAY_VARIABLES \
GENINFO_NANOAOD_ARRAY_VARIABLE(float, LHEScaleWeight, 0, 50) \
GENINFO_NANOAOD_ARRAY_VARIABLE(float, LHEPdfWeight, 0, 500) \
GENINFO_NANOAOD_ARRAY_VARIABLE(float, PSWeight, 0, 500)

#define GENINFO_NANOAOD_ALLVARIABLES \
GENINFO_NANOAOD_SCALAR_VARIABLES \
GENINFO_NANOAOD_ARRAY_VARIABLES

// Framework translation
#define GENINFO_EXTRA_VARIABLES \
GENINFO_VARIABLE(int, genTtbarId, 0) \
GENINFO_VARIABLE(float, genHEPMCweight, 1) \
GENINFO_VARIABLE(float, LHEweight_unscaledOriginalWeight, 1) /* = LHEEventProduct::originalXWGTUP() */ \
GENINFO_VARIABLE(float, LHEweight_defaultMemberZero, 1) /* = PDFweights[0] */ \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR1_muF1, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR1_muF2, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR1_muF0p5, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR2_muF1, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR2_muF2, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR2_muF0p5, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR0p5_muF1, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR0p5_muF2, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR0p5_muF0p5, 1) \
GENINFO_VARIABLE(float, LHEweight_PDFVariation_Up, 1) \
GENINFO_VARIABLE(float, LHEweight_PDFVariation_Dn, 1) \
GENINFO_VARIABLE(float, LHEweight_AsMZ_Up, 1) \
GENINFO_VARIABLE(float, LHEweight_AsMZ_Dn, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muRoneoversqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muRoneoversqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muRsqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muRsqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR0p5, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR0p5, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR0p25, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR0p25, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR4, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR4, 1)


class GenInfoVariables{
public:
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  GENINFO_EXTRA_VARIABLES;
#undef GENINFO_VARIABLE

  std::unordered_map<TString, float> LHE_ME_weights;
  std::unordered_map<TString, float> Kfactors;

  GenInfoVariables();
  GenInfoVariables(GenInfoVariables const& other);
  GenInfoVariables& operator=(const GenInfoVariables& other);

  void swap(GenInfoVariables& other);

};

class GenInfoObject{
public:
  GenInfoVariables extras;

public:
  GenInfoObject();
  GenInfoObject(const GenInfoObject& other);
  GenInfoObject& operator=(const GenInfoObject& other);
  ~GenInfoObject();

  void swap(GenInfoObject& other);

#define GENINFO_NANOAOD_SCALAR_VARIABLE(TYPE, NAME, DEFVAL) , TYPE const& NAME
#define GENINFO_NANOAOD_ARRAY_VARIABLE(TYPE, NAME, DEFVAL, MAXSIZE) , unsigned int const& n##NAME, TYPE const* arr_##NAME
  void acquireGenInfo(
    TString const& sampleIdentifier
    GENINFO_NANOAOD_ALLVARIABLES
  );
#undef GENINFO_NANOAOD_ARRAY_VARIABLE
#undef GENINFO_NANOAOD_SCALAR_VARIABLE

  float getGenWeight(SystematicsHelpers::SystematicVariationTypes const& syst) const;
};

#endif
