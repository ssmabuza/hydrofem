// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_GlobalGather_HPP__
#define __Hydrofem_GlobalGather_HPP__

#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

class GlobalGather
{
public:

  GlobalGather(const std::shared_ptr<DofMapper>& dofmapper);

  void doGather(FEArray<double>::CellBasis& local,
                const std::shared_ptr<const FEVector>& global);

protected:

  std::shared_ptr<DofMapper> m_dofmapper;
  
};



}
// end namespace hydrofem


#endif /** __Hydrofem_GlobalGather_HPP__ */
