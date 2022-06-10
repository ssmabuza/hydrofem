// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_NodalGradients_HPP__
#define __Hydrofem_NodalGradients_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

/**
 * \brief A routine that builds the lumped mass nodal gradients for u
 * 
 * \param [out] grad_u the lumped mass L2 projection of gradient of u
 *
 * \param [in] conv_mats the global convection matrices
 * \param [in] lumped_mass_mat the global lumped mass matrix
 * \param [in] u the solution
 */
void ComputeNodalGradients(std::vector<std::shared_ptr<FEVector>>& grad_u,
                           const std::vector<std::shared_ptr<const FEMatrix>>& conv_mats,
                           const std::shared_ptr<const FEVector>& lumped_mass_mat,
                           const std::shared_ptr<const FEVector>& u);


}
// end namespace hydrofem

#endif /** __Hydrofem_NodalGradients_HPP__ */
