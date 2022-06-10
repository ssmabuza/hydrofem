// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_AssemblerFactory.hpp"

#include "Hydrofem_Assembler_Poisson.hpp"
//#include "Hydrofem_Assembler_Hyperbolic.hpp"
#include "Hydrofem_Assembler_Bioseparation.hpp"

namespace hydrofem
{

std::shared_ptr<Assembler_Base>
AssemblerFactory::build() const
{
  std::shared_ptr<Assembler_Base> assembler;
  if (m_name == "poisson")
    assembler = std::make_shared<Assembler_Poisson>(m_option_handler);
  else if (m_name == "bioseparation")
    assembler = std::make_shared<Assembler_Bioseparation>(m_option_handler);
  assert(assembler);
  return assembler;
}
  
  
}
// end namespace hydrofem
