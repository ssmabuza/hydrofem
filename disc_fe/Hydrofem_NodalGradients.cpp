// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_NodalGradients.hpp"

namespace hydrofem
{

void ComputeNodalGradients(std::vector<std::shared_ptr<FEVector>>& grad_u,
                           const std::vector<std::shared_ptr<const FEMatrix>>& conv_mats,
                           const std::shared_ptr<const FEVector>& lumped_mass_mat,
                           const std::shared_ptr<const FEVector>& u)
{
  
  assert(grad_u.size()==conv_mats.size());
  for (std::size_t dim = 0; dim < grad_u.size(); ++dim)
  {
    grad_u.at(dim)->setConstant(0.0);
    (*grad_u.at(dim)) = (*conv_mats.at(dim)) * (*u);
    (*grad_u.at(dim)) *= *lumped_mass_mat;
  }
  
}

}
// end namespace hydrofem
