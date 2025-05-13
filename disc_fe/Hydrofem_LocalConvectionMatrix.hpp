// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_LocalConvectionMatrix_HPP__
#define __Hydrofem_LocalConvectionMatrix_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

/**
 * \brief integral over K of \phi_i\grad\phi_j
 * 
 * \param [out] loc_K local convection matrix
 * \param [in]  mesh global 2D mesh
 * \param [in]  elem_ind element index
 * 
 */
void LocalConvectionMatrix(LMAT_<double>& loc_K,
                           const std::vector<int>& loc_ind,
                           const std::vector<SPoint>& element,
                           const std::shared_ptr<Quadrature>& quadrature,
                           const std::vector<std::shared_ptr<FEBasis>>& fe_basis,
                           const std::function<SPoint(SPoint)>& velocity);

}
// end namespace hydrofem

#endif /** __Hydrofem_LocalConvectionMatrix_HPP__ */
