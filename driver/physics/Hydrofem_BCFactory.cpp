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

#include "Hydrofem_BC_CDR.hpp"
#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_BC_Bioseparation.hpp"

namespace hydrofem 
{

std::shared_ptr<BC> BCFactory::build()
{
  std::shared_ptr<BC> bc;
  
  if (m_bc_name == "bc-scalar")
  {
    bc = std::make_shared<BC_Scalar>(m_mesh);
  } else if (m_bc_name == "bc-bioseparation") {
    bc = std::make_shared<BC_Bioseparation>(m_mesh);
  } else if (m_bc_name == "bc-cdr") {
    bc = std::make_shared<BC_CDR>(m_mesh);
  } else {
    throw std::logic_error("Error in BCFactory::buildBC(), no valid BC chosen.");
  }
    
  assert(bc);
  return bc;
  
}

}
// end namespace hydrofem
