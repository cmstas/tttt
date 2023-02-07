#include <cassert>

#include "IvyFramework/IvyDataTools/interface/ParticleObjectHelpers.h"
#include "GlobalCollectionNames.h"
#include "SamplesCore.h"
#include "MuonHandler.h"
#include "MuonSelectionHelpers.h"


using namespace std;
using namespace IvyStreamHelpers;


#define MUON_MOMENTUM_VARIABLES \
MUON_VARIABLE(float, pt, 0) \
MUON_VARIABLE(float, eta, 0) \
MUON_VARIABLE(float, phi, 0) \
MUON_VARIABLE(float, mass, 0) \
MUON_VARIABLE(int, pdgId, -9000)


const std::string MuonHandler::colName = GlobalCollectionNames::colName_muons;

MuonHandler::MuonHandler() :
  IvyBase(),
  extractOptionalInfo(false)
{
  setConsumed(0);
}

void MuonHandler::setConsumed(unsigned char i_used_unused){
  auto const& dy = SampleHelpers::getDataYear();
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE* const>(MuonHandler::colName + "_" + #NAME);
  if (i_used_unused==0){
    this->addConsumed<GlobalCollectionNames::collsize_t>(Form("n%s", MuonHandler::colName.data()));
    MUON_MOMENTUM_VARIABLES;
    MUON_EXTRA_MANDATORY_VARIABLES_COMMON;
    if (dy>=2015 && dy<=2018){
      MUON_EXTRA_MANDATORY_VARIABLES_RUN2;
    }
    else if (dy==2022){
      MUON_EXTRA_MANDATORY_VARIABLES_RUN3;
    }
    else{
      IVYerr << "MuonHandler::setConsumed: Could not identify data year to determine year-dependent mandatory variable names." << endl;
      assert(0);
    }
  }
  else{
    MUON_EXTRA_OPTIONAL_VARIABLES_COMMON;
    if (dy>=2015 && dy<=2018){
      MUON_EXTRA_OPTIONAL_VARIABLES_RUN2;
    }
    else if (dy==2022){
      MUON_EXTRA_OPTIONAL_VARIABLES_RUN3;
    }
    else{
      IVYerr << "MuonHandler::setConsumed: Could not identify data year to determine year-dependent optional variable names." << endl;
      assert(0);
    }
  }
#undef MUON_VARIABLE
}

bool MuonHandler::constructMuons(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = constructMuonObjects();

  if (res) this->cacheEvent();
  return res;
}

bool MuonHandler::constructMuonObjects(){
  auto const& dy = SampleHelpers::getDataYear();

  GlobalCollectionNames::collsize_t nProducts;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) TYPE* const* arr_##NAME = nullptr;
  MUON_MOMENTUM_VARIABLES;
  MUON_EXTRA_VARIABLES;
#undef MUON_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = this->getConsumedValue(Form("n%s", MuonHandler::colName.data()), nProducts);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed<TYPE* const>(MuonHandler::colName + "_" + #NAME, arr_##NAME);
  MUON_MOMENTUM_VARIABLES;
  MUON_EXTRA_MANDATORY_VARIABLES_COMMON;
  if (extractOptionalInfo){
    MUON_EXTRA_OPTIONAL_VARIABLES_COMMON;
  }
  if (dy>=2015 && dy<=2018){
    MUON_EXTRA_MANDATORY_VARIABLES_RUN2;
    if (extractOptionalInfo){
      MUON_EXTRA_OPTIONAL_VARIABLES_RUN2;
    }
  }
  else if (dy==2022){
    MUON_EXTRA_MANDATORY_VARIABLES_RUN3;
    if (extractOptionalInfo){
      MUON_EXTRA_OPTIONAL_VARIABLES_RUN3;
    }
  }
  else{
    IVYerr << "MuonHandler::constructMuonObjects: Could not identify data year to determine year-dependent variable names." << endl;
    assert(0);
  }
#undef MUON_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "MuonHandler::constructMuonObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "MuonHandler::constructMuonObjects: All variables are set up!" << endl;

  if (nProducts==0) return true; // Construction is successful, it is just that no muons exist.

  productList.reserve(nProducts);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) TYPE* it_##NAME = nullptr; if (arr_##NAME) it_##NAME = &((*arr_##NAME)[0]);
  MUON_MOMENTUM_VARIABLES;
  MUON_EXTRA_VARIABLES
#undef MUON_VARIABLE
  {
    GlobalCollectionNames::collsize_t ip=0;
    while (ip != nProducts){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "MuonHandler::constructMuonObjects: Attempting muon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line...
      productList.push_back(new MuonObject(*it_pdgId, momentum));
      MuonObject*& obj = productList.back();

      // Set extras
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) if (it_##NAME) obj->extras.NAME = *it_##NAME;
      MUON_EXTRA_VARIABLES;
#undef MUON_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Set selection bits
      //MuonSelectionHelpers::setSelectionBits(*obj); // Do not set them here. Migrated to ParticleDisambiguator.

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) if (it_##NAME) it_##NAME++;
      MUON_MOMENTUM_VARIABLES;
      MUON_EXTRA_VARIABLES;
#undef MUON_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void MuonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  if (extractOptionalInfo) setConsumed(1);

  auto const& dy = SampleHelpers::getDataYear();

  tree->bookBranch<GlobalCollectionNames::collsize_t>(Form("n%s", MuonHandler::colName.data()), 0);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookArrayBranch<TYPE>(MuonHandler::colName + "_" + #NAME, DEFVAL, GlobalCollectionNames::colMaxSize_muons);
  MUON_MOMENTUM_VARIABLES;
  MUON_EXTRA_MANDATORY_VARIABLES_COMMON;
  if (extractOptionalInfo){
    MUON_EXTRA_OPTIONAL_VARIABLES_COMMON;
  }
  if (dy>=2015 && dy<=2018){
    MUON_EXTRA_MANDATORY_VARIABLES_RUN2;
    if (extractOptionalInfo){
      MUON_EXTRA_OPTIONAL_VARIABLES_RUN2;
    }
  }
  else if (dy==2022){
    MUON_EXTRA_MANDATORY_VARIABLES_RUN3;
    if (extractOptionalInfo){
      MUON_EXTRA_OPTIONAL_VARIABLES_RUN3;
    }
  }
  else{
    IVYerr << "MuonHandler::bookBranches: Could not identify data year to determine year-dependent variable names." << endl;
    assert(0);
  }
#undef MUON_VARIABLE
}


#undef MUON_MOMENTUM_VARIABLES
