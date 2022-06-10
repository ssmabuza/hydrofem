// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IOFactory_HPP__
#define __Hydrofem_IOFactory_HPP__

#include "Hydrofem_IOBase.hpp"
#include "Hydrofem_FEBasis.hpp"

namespace hydrofem
{

class DofMapper;

class IOFactory
{
public:

  IOFactory() {}

  virtual ~IOFactory() {}
  
  std::shared_ptr<IOBase> buildIO(const std::shared_ptr<DofMapper>& dofmapper,
                                  const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                                  bool xml = false, bool steady = false);
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_IOFactory_HPP__ */
