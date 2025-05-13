// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_Assembler_Bioseparation_MatrixForm.hpp"


namespace hydrofem
{

void Assembler_Bioseparation_MatrixForm::
buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                         const std::shared_ptr<const FEVector>& U_dot,
                         const std::shared_ptr<FEVector>& res_U,
                         const std::shared_ptr<FEMatrix>& jac_U,
                         const double time,
                         const double delta_t,
                         const double beta) const
{

  // compute/update the mass matrix based on the current solution U

  // compute the stiffness matrix once

  // compute the convection matrix once

  // build F(U) and J(U) without AFC using powerfull Eigen ops




}


}
// end namespace hydrofem
