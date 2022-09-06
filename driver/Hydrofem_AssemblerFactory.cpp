// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_AssemblerFactory.hpp"

#include "Hydrofem_Assembler_Poisson.hpp"
#include "Hydrofem_Assembler_Bioseparation.hpp"

namespace hydrofem
{

std::shared_ptr<Assembler_Base>
AssemblerFactory::build() const
{
  std::shared_ptr<Assembler_Base> assembler;
  if (m_name == "poisson")
    assembler = std::make_shared<Assembler_Poisson>(m_problem,m_dofmapper,m_basis,m_quadrature);
  else if (m_name == "bioseparation")
    assembler = std::make_shared<Assembler_Bioseparation>(m_problem,m_dofmapper,m_basis,m_quadrature,m_afc_enabled);
  assert(assembler);
  return assembler;
}

}
// end namespace hydrofem
