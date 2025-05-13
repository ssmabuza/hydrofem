// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_GlobalGather.hpp"

namespace hydrofem
{

GlobalGather::
GlobalGather(const std::shared_ptr<DofMapper>& dofmapper)
{
  m_dofmapper = dofmapper;
}

void GlobalGather::
doGather(FEArray<double>::CellBasis& local,
         const std::shared_ptr<const FEVector>& global)
{
  const auto& lids = m_dofmapper->getLocDofIndexes();
  for (int cell = 0; cell < m_dofmapper->nelements(); ++cell)
  {
    const auto& glids = m_dofmapper->getGlobDofIndexes(cell);
    for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
      local(cell,lids[i]) = (*global)[glids[i]];
  }
}

}
// end namespace hydrofem
