#include <cassert>

#include "IvyFramework/IvyDataTools/interface/ParticleObjectHelpers.h"
#include "GlobalCollectionNames.h"
#include "SamplesCore.h"
#include "ElectronHandler.h"
#include "ElectronSelectionHelpers.h"


using namespace std;
using namespace IvyStreamHelpers;


#define ELECTRON_MOMENTUM_VARIABLES \
ELECTRON_VARIABLE(float, pt, 0) \
ELECTRON_VARIABLE(float, eta, 0) \
ELECTRON_VARIABLE(float, phi, 0) \
ELECTRON_VARIABLE(float, mass, 0) \
ELECTRON_VARIABLE(int, pdgId, -9000)


const std::string ElectronHandler::colName = GlobalCollectionNames::colName_electrons;

ElectronHandler::ElectronHandler() :
  IvyBase(),
  extractOptionalInfo(false)
{
  setConsumed(0);
}

void ElectronHandler::setConsumed(unsigned char i_used_unused){
  auto const& dy = SampleHelpers::getDataYear();
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE* const>(ElectronHandler::colName + "_" + #NAME);
  if (i_used_unused==0){
    this->addConsumed<GlobalCollectionNames::collsize_t>(Form("n%s", ElectronHandler::colName.data()));
    ELECTRON_MOMENTUM_VARIABLES;
    ELECTRON_EXTRA_MANDATORY_VARIABLES_COMMON;
    if (dy>=2015 && dy<=2018){
      ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN2;
    }
    else if (dy==2022){
      ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN3;
    }
    else{
      IVYerr << "ElectronHandler::setConsumed: Could not identify data year to determine year-dependent mandatory variable names." << endl;
      assert(0);
    }
  }
  else{
    ELECTRON_EXTRA_OPTIONAL_VARIABLES_COMMON;
    if (dy>=2015 && dy<=2018){
      ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN2;
    }
    else if (dy==2022){
      ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN3;
    }
    else{
      IVYerr << "ElectronHandler::setConsumed: Could not identify data year to determine year-dependent optional variable names." << endl;
      assert(0);
    }
  }
#undef ELECTRON_VARIABLE
}


bool ElectronHandler::constructElectrons(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = constructElectronObjects();

  if (res) this->cacheEvent();
  return res;
}

bool ElectronHandler::constructElectronObjects(){
  auto const& dy = SampleHelpers::getDataYear();

  GlobalCollectionNames::collsize_t nProducts;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) TYPE* const* arr_##NAME = nullptr;
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_EXTRA_VARIABLES;
#undef ELECTRON_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = this->getConsumedValue(Form("n%s", ElectronHandler::colName.data()), nProducts);
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed<TYPE* const>(ElectronHandler::colName + "_" + #NAME, arr_##NAME);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_EXTRA_MANDATORY_VARIABLES_COMMON;
  if (extractOptionalInfo){
    ELECTRON_EXTRA_OPTIONAL_VARIABLES_COMMON;
  }
  if (dy>=2015 && dy<=2018){
    ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN2;
    if (extractOptionalInfo){
      ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN2;
    }
  }
  else if (dy==2022){
    ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN3;
    if (extractOptionalInfo){
      ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN3;
    }
  }
  else{
    IVYerr << "ElectronHandler::constructElectronObjects: Could not identify data year to determine year-dependent variable names." << endl;
    assert(0);
  }
#undef ELECTRON_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "ElectronHandler::constructElectronObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "ElectronHandler::constructElectronObjects: All variables are set up!" << endl;

  if (nProducts==0) return true; // Construction is successful, it is just that no electrons exist.

  productList.reserve(nProducts);
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) TYPE* it_##NAME = (arr_##NAME ? &((*arr_##NAME)[0]) : nullptr);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_EXTRA_VARIABLES
#undef ELECTRON_VARIABLE
  {
    GlobalCollectionNames::collsize_t ip=0;
    while (ip != nProducts){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "ElectronHandler::constructElectronObjects: Attempting electron " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line...
      productList.push_back(new ElectronObject(*it_pdgId, momentum));
      ElectronObject*& obj = productList.back();

      // Set extras
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) if (it_##NAME) obj->extras.NAME = *it_##NAME;
      ELECTRON_EXTRA_VARIABLES;
#undef ELECTRON_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Set selection bits
      //ElectronSelectionHelpers::setSelectionBits(*obj); // Do not set them here. Migrated to ParticleDisambiguator.

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) if (it_##NAME) it_##NAME++;
      ELECTRON_MOMENTUM_VARIABLES;
      ELECTRON_EXTRA_VARIABLES;
#undef ELECTRON_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void ElectronHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  if (extractOptionalInfo) setConsumed(1);

  auto const& dy = SampleHelpers::getDataYear();

  tree->bookBranch<GlobalCollectionNames::collsize_t>(Form("n%s", ElectronHandler::colName.data()), 0);
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookArrayBranch<TYPE>(ElectronHandler::colName + "_" + #NAME, DEFVAL, GlobalCollectionNames::colMaxSize_electrons);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_EXTRA_MANDATORY_VARIABLES_COMMON;
  if (extractOptionalInfo){
    ELECTRON_EXTRA_OPTIONAL_VARIABLES_COMMON;
  }
  if (dy>=2015 && dy<=2018){
    ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN2;
    if (extractOptionalInfo){
      ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN2;
    }
  }
  else if (dy==2022){
    ELECTRON_EXTRA_MANDATORY_VARIABLES_RUN3;
    if (extractOptionalInfo){
      ELECTRON_EXTRA_OPTIONAL_VARIABLES_RUN3;
    }
  }
  else{
    IVYerr << "ElectronHandler::bookBranches: Could not identify data year to determine year-dependent variable names." << endl;
    assert(0);
  }
#undef ELECTRON_VARIABLE
}


#undef ELECTRON_MOMENTUM_VARIABLES
