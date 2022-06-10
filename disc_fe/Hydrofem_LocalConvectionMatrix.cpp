// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_LocalConvectionMatrix.hpp"

namespace hydrofem
{

void LocalConvectionMatrix(LMAT_<double>& loc_K,
                           const std::vector<int>& loc_ind,
                           const std::vector<SPoint>& element,
                           const std::shared_ptr<Quadrature>& quadrature,
                           const std::vector<std::shared_ptr<FEBasis>>& fe_basis,
                           const std::function<SPoint(SPoint)>& velocity)
{
  const auto& quad_e = *quadrature;
  const auto& q_points = quad_e.get_q_points();
  const auto& q_weights = quad_e.get_q_weights();
  const auto& fe_basis_e = fe_basis;
  
  for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
  {
    const int& loc_i = loc_ind.at(dof_i);
    for (std::size_t dof_j = 0; dof_j < loc_ind.size(); ++dof_j)
    {
      const int& loc_j = loc_ind.at(dof_j);
      double val = 0.0;
      for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
      {
        const double phi_j = (*fe_basis_e.at(loc_j))(q_points.at(quad_ind),element);
        const SPoint grad_phi_i = (*fe_basis_e.at(loc_i)).grad(q_points.at(quad_ind),element);
        val += - q_weights.at(quad_ind) * phi_j * (grad_phi_i * velocity(q_points.at(quad_ind)));
      }
      loc_K(loc_i,loc_j) = val;
    }
  }
}

}
// end namespace hydrofem
