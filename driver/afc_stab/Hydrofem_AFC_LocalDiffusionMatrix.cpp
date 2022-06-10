// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_AFC_LocalDiffusionMatrix.hpp"

namespace hydrofem
{

void AFC_LocalDiffusionMatrix(LMAT_<double>& loc_D,
                              const LMAT_<double>& loc_K)
{
  for (Eigen::Index i = 0; i < loc_K.dimension(0); ++i)
  {
    for (Eigen::Index j = 0; j <= i; ++j)
    {
      if (i != j)
      {
        const double d_ij = std::max(-loc_K(i, j), std::max(0.0, -loc_K(j, i)));
        loc_D(i,i) -= d_ij;
        loc_D(j,j) -= d_ij;
        loc_D(i,j)  = d_ij;
        loc_D(j,i)  = d_ij;
      }
    }
  }
}

}
// end namespace hydrofem
