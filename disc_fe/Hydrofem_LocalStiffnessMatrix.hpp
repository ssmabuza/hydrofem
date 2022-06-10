// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_LocalStiffnessMatrix_HPP__
#define __Hydrofem_LocalStiffnessMatrix_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

/**
 * \brief integral over K of \diff\grad\phi_i\grad\phi_j
 *
 * \param [out] loc_S local stiffness matrix
 * \param [in]  elem_ind element index
 *
 */
void LocalStiffnessMatrix(LMAT_<double>& loc_S,
                          const std::vector<int>& loc_ind,
                          const std::vector<SPoint>& element,
                          const std::shared_ptr<Quadrature>& quadrature,
                          const std::vector<std::shared_ptr<FEBasis>>& fe_basis,
                          const double diff);

/**
 * \brief integral over K of \grad\phi_i \D \grad\phi_j
 *
 * \param [out] loc_S local stiffness matrix
 * \param [in]  elem_ind element index
 *
 */
void LocalStiffnessMatrix(LMAT_<double>& loc_S,
                          const std::vector<int>& loc_ind,
                          const std::vector<SPoint>& element,
                          const std::shared_ptr<Quadrature>& quadrature,
                          const std::vector<std::shared_ptr<FEBasis>>& fe_basis,
                          const std::function<LMAT_<double>(SPoint)>& diff);

}
// end namespace hydrofem

#endif /** __Hydrofem_LocalStiffnessMatrix_HPP__ */
