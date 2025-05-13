// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_GlobalConvectionMatrix_HPP__
#define __Hydrofem_GlobalConvectionMatrix_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"

namespace hydrofem
{

void GlobalConvectionMatrices(const std::vector<std::shared_ptr<FEMatrix>>& conv_mats,
                              const std::shared_ptr<const DofMapper>& dofmapper,
                              const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                              const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature);

}
// end namespace hydrofem

#endif /** __Hydrofem_GlobalConvectionMatrix_HPP__ */
