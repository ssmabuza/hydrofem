// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_BC.hpp"
#include "Hydrofem_BCFactory.hpp"

#include "Hydrofem_Mesh.hpp"

#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_BC_ReactiveFlow.hpp"
#include "Hydrofem_BC_EulerEquations.hpp"

namespace hydrofem 
{

std::shared_ptr<BC> BCFactory::build()
{
  std::shared_ptr<BC> bc;
  
  if (m_bc_name == "bc-scalar")
  {
    
    bc = std::make_shared<BC_Scalar>(m_mesh);
    
  } else if (m_bc_name == "bc-euler-equations") {
    
    bc = std::make_shared<BC_EulerEquations>(m_mesh);
    
  } else if (m_bc_name == "bc-reactive-flow") {
    
    bc = std::make_shared<BC_ReactiveFlow>(m_mesh);
    
  } else {
    
    throw std::logic_error("Error in BCFactory::buildBC(), no valid BC chosen.");
    
  }
    
  assert(bc);
  return bc;
  
}

  
  
}
// end namespace valiant
