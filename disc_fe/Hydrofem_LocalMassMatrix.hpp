// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_LocalMassMatrix_HPP__
#define __Hydrofem_LocalMassMatrix_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

void LocalMassMatrix(LMAT_<double>& loc_Mc,
                     const std::vector<int>& loc_ind,
                     const std::vector<SPoint>& element,
                     const std::shared_ptr<Quadrature>& quadrature,
                     const std::vector<std::shared_ptr<FEBasis>>& fe_basis);

void LocalLumpedMassMatrix(LMAT_<double>& loc_Ml,
                           const std::vector<int>& loc_ind,
                           const std::vector<SPoint>& element,
                           const std::shared_ptr<Quadrature>& quadrature,
                           const std::vector<std::shared_ptr<FEBasis>>& fe_basis);

}
// end namespace hydrofem

#endif /** __Hydrofem_LocalMassMatrix_HPP__ */
