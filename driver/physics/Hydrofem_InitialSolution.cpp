// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_GlobalConstants.hpp"
#include "Hydrofem_InitialSolution.hpp"

namespace hydrofem
{

// TODO: replace type with Kokkos::View<MDArray*,PHX::Device>
std::shared_ptr<std::vector<InitialSolution::MDArray>>
InitialSolution::
get_evaluatedField() const
{
  assert(m_is_computed);
  return m_result_gathered;
}

std::shared_ptr<FEVector>
InitialSolution::
get_evaluatedFieldGlobal() const
{
  assert(m_is_computed);
  return m_result;
}

std::shared_ptr<FEMatrix>
InitialSolution::
getConsistentMassMatrix() const
{
  return m_mass;
}

std::shared_ptr<FEVector>
InitialSolution::
getLumpedMassMatrix() const
{
  return m_lumped_mass;
}

}
// end namespace valiant
