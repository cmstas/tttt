#ifndef MUONHANDLER_H
#define MUONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ParticleDisambiguator.h"


class MuonHandler : public IvyBase{
public:
  typedef MuonObject ProductType_t;

  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool extractOptionalInfo;
  std::vector<ProductType_t*> productList;

  void clear(){ this->resetCache(); for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  void setConsumed(unsigned char i_used_unused);

  bool constructMuonObjects();

public:
  // Constructors
  MuonHandler();

  // Destructors
  ~MuonHandler(){ clear(); }

  bool constructMuons();

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  void bookBranches(BaseTree* tree);

  void setExtractOptionalInfo(bool flag){ extractOptionalInfo = flag; }

};


#endif
