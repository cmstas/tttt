#include <cassert>
#include <algorithm>
#include <utility>
#include <cmath>
#include "GenInfoObject.h"
#include "HelperFunctions.h"
#include "IvyFramework/IvyDataTools/interface/IvyStreamHelpers.hh"


using namespace std;
using namespace IvyStreamHelpers;


GenInfoVariables::GenInfoVariables(){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  GENINFO_EXTRA_VARIABLES;
#undef GENINFO_VARIABLE
}
GenInfoVariables::GenInfoVariables(GenInfoVariables const& other){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  GENINFO_EXTRA_VARIABLES;
#undef GENINFO_VARIABLE

  this->LHE_ME_weights = other.LHE_ME_weights;
  this->Kfactors = other.Kfactors;
}
void GenInfoVariables::swap(GenInfoVariables& other){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  GENINFO_EXTRA_VARIABLES;
#undef GENINFO_VARIABLE

  std::swap(this->LHE_ME_weights, other.LHE_ME_weights);
  std::swap(this->Kfactors, other.Kfactors);
}
GenInfoVariables& GenInfoVariables::operator=(const GenInfoVariables& other){
  GenInfoVariables tmp(other);
  swap(tmp);
  return *this;
}

GenInfoObject::GenInfoObject() :
  extras()
{}
GenInfoObject::GenInfoObject(const GenInfoObject& other) :
  extras(other.extras)
{}
void GenInfoObject::swap(GenInfoObject& other){
  extras.swap(other.extras);
}
GenInfoObject& GenInfoObject::operator=(const GenInfoObject& other){
  GenInfoObject tmp(other);
  swap(tmp);
  return *this;
}
GenInfoObject::~GenInfoObject(){}

float GenInfoObject::getGenWeight(SystematicsHelpers::SystematicVariationTypes const& syst) const{
  using namespace SystematicsHelpers;
  float wgt = extras.genHEPMCweight;
  switch (syst){
  case tPDFScaleDn:
    wgt *= extras.LHEweight_QCDscale_muR1_muF0p5;
    break;
  case tPDFScaleUp:
    wgt *= extras.LHEweight_QCDscale_muR1_muF2;
    break;
  case tQCDScaleDn:
    wgt *= extras.LHEweight_QCDscale_muR0p5_muF1;
    break;
  case tQCDScaleUp:
    wgt *= extras.LHEweight_QCDscale_muR2_muF1;
    break;
  case tAsMZDn:
    wgt *= extras.LHEweight_AsMZ_Dn;
    break;
  case tAsMZUp:
    wgt *= extras.LHEweight_AsMZ_Up;
    break;
  case tPDFReplicaDn:
    wgt *= extras.LHEweight_PDFVariation_Dn;
    break;
  case tPDFReplicaUp:
    wgt *= extras.LHEweight_PDFVariation_Up;
    break;
  case tPythiaScaleDn:
    wgt *= extras.PythiaWeight_isr_muR0p5 * extras.PythiaWeight_fsr_muR0p5;
    break;
  case tPythiaScaleUp:
    wgt *= extras.PythiaWeight_isr_muR2 * extras.PythiaWeight_fsr_muR2;
    break;
  default:
    break;
  }

  if (!HelperFunctions::checkVarNanInf(wgt)){
    IVYout << "GenInfoObject::getGenWeight(" << syst << "): Weight is " << wgt << "." << endl;
    assert(0);
  }

  return wgt;
}

#define GENINFO_NANOAOD_SCALAR_VARIABLE(TYPE, NAME, DEFVAL) , TYPE const& NAME
#define GENINFO_NANOAOD_ARRAY_VARIABLE(TYPE, NAME, DEFVAL, MAXSIZE) , unsigned int const& n##NAME, TYPE const* arr_##NAME
void GenInfoObject::acquireGenInfo(
  TString const& sampleIdentifier
  GENINFO_NANOAOD_ALLVARIABLES
){
#undef GENINFO_NANOAOD_ARRAY_VARIABLE
#undef GENINFO_NANOAOD_SCALAR_VARIABLE

  extras.genTtbarId = genTtbarId;

  // FIXME: I don't know the conventions just yet, so for now, set all systematics to 1, and central weight to wgts[0]
  extras.genHEPMCweight = genWeight;
  extras.LHEweight_unscaledOriginalWeight = LHEWeight_originalXWGTUP;
  // Weight vairations are supposed to be ratios to genHEPMCweight!
  // BE CAREFUL ABOUT HOW TO TAKE RATIOS!!!
}
