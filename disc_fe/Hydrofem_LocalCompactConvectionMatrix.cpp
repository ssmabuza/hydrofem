// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

// #include <KokkosBlas3_gemm.hpp>
// #include <Teuchos_SerialDenseSolver.hpp>
#include "Hydrofem_LocalCompactConvectionMatrix.hpp"

namespace hydrofem
{

void LocalCompactConvectionMatrix(LMAT_<double>& comp_loc_K,
                                  const LMAT_<double>& loc_K,
                                  const LMAT_<double>& loc_mass,
                                  const LMAT_<double>& loc_lumped_mass)
{
  // deep copy local mass matrix to teuchos serial dense matrix
//   Teuchos::SerialDenseMatrix<int,double> tsd_loc_mass(loc_mass.extent(0),loc_mass.extent(1));
//   for (int i = 0; i < tsd_loc_mass.numRows(); ++i)
//     for (int j = 0; j < tsd_loc_mass.numCols(); ++j)
//       tsd_loc_mass(i,j) = loc_mass(i,j);
//   auto inv_loc_mass = loc_mass.inverse();
  // invert the matrix
//   Teuchos::SerialDenseSolver<int,double> solver;
//   solver.setMatrix(Teuchos::rcpFromRef(tsd_loc_mass));
//   solver.invert();
  // copy into a kokkos view
//   LMAT_<double> inv_loc_mass("inv_loc_mass",tsd_loc_mass.numRows(),tsd_loc_mass.numCols());
//   for (int i = 0; i < tsd_loc_mass.numRows(); ++i)
//     for (int j = 0; j < tsd_loc_mass.numCols(); ++j)
//       inv_loc_mass(i,j) = tsd_loc_mass(i,j);
  // perform mat mat multiplication Kokkos-BLAS
  comp_loc_K = loc_lumped_mass * (loc_mass.inverse() * loc_K);
//   LMAT_<double> tmp("tmp",tsd_loc_mass.numRows(),tsd_loc_mass.numCols());
//   KokkosBlas::gemm("N","N",1.0, inv_loc_mass,loc_K,0.0,tmp);
//   KokkosBlas::gemm("N","N",1.0, loc_lumped_mass,tmp,0.0,comp_loc_K);
}

}
// end namespace hydrofem
