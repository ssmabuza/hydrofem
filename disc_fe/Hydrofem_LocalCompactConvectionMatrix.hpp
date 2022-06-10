// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_LocalCompactConvectionMatrix_HPP__
#define __Hydrofem_LocalCompactConvectionMatrix_HPP__

#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

void LocalCompactConvectionMatrix(LMAT_<double>& comp_loc_K,
                                  const LMAT_<double>& loc_K,
                                  const LMAT_<double>& loc_mass,
                                  const LMAT_<double>& loc_lumped_mass);

}
// end namespace hydrofem

#endif /** __Hydrofem_LocalCompactConvectionMatrix_HPP__ */
