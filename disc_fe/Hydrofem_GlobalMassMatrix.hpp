// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_GlobalMassMatrix_HPP__
#define __Hydrofem_GlobalMassMatrix_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

void GlobalMassMatrix(const std::shared_ptr<FEMatrix>& mass_mat, 
                      const std::shared_ptr<const DofMapper>& dofmapper,
                      const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                      const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature);

void GlobalLumpedMassMatrix(const std::shared_ptr<FEVector>& lumped_mass_mat,
                            const std::shared_ptr<const DofMapper>& dofmapper,
                            const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                            const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature);

}
// end namespace hydrofem

#endif /** __Hydrofem_GlobalMassMatrix_HPP__ */
