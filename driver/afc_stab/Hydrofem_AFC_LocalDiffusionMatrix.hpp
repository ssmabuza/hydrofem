// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_AFC_LocalDiffusionMatrix_HPP__
#define __Hydrofem_AFC_LocalDiffusionMatrix_HPP__

#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

void AFC_LocalDiffusionMatrix(LMAT_<double>& loc_D,
                              const LMAT_<double>& loc_K);

}
// end namespace hydrofem

#endif /** __Hydrofem_AFC_LocalDiffusionMatrix_HPP__ */
