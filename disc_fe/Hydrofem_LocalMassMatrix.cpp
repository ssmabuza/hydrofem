// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_LocalMassMatrix.hpp"

namespace hydrofem
{

void LocalMassMatrix(LMAT_<double>& loc_Mc,
                     const std::vector<int>& loc_ind,
                     const std::vector<SPoint>& element,
                     const std::shared_ptr<Quadrature>& quadrature,
                     const std::vector<std::shared_ptr<FEBasis>>& fe_basis)
{
  const auto& quad_e = *quadrature;
  const auto& q_points = quad_e.get_q_points();
  const auto& q_weights = quad_e.get_q_weights();
  const auto& fe_basis_e = fe_basis;
  
  for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
  {
    const int& loc_i = loc_ind.at(dof_i);
    for (std::size_t dof_j = 0; dof_j <= dof_i; ++dof_j)
    {
      const int& loc_j = loc_ind.at(dof_j);
      double val = 0.0;
      for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
        val += q_weights.at(quad_ind) * (*fe_basis_e.at(loc_i))(q_points.at(quad_ind),element) * (*fe_basis_e.at(loc_j))(q_points.at(quad_ind),element);
      loc_Mc(loc_i,loc_j) = val;
      if (loc_i!=loc_j)
        loc_Mc(loc_j,loc_i) = val;
    }
  }
}

void LocalLumpedMassMatrix(LMAT_<double>& loc_Ml,
                           const std::vector<int>& loc_ind,
                           const std::vector<SPoint>& element,
                           const std::shared_ptr<Quadrature>& quadrature,
                           const std::vector<std::shared_ptr<FEBasis>>& fe_basis)
{
  const auto& quad_e = *quadrature;
  const auto& q_points = quad_e.get_q_points();
  const auto& q_weights = quad_e.get_q_weights();
  const auto& fe_basis_e = fe_basis;
  
  for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
  {
    const int& loc_i = loc_ind.at(dof_i);
    double val = 0.0;
    for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
      val += q_weights.at(quad_ind) * (*fe_basis_e.at(loc_i))(q_points.at(quad_ind),element);
    loc_Ml(loc_i,loc_i) = val;
  }
}

}
// end namespace hydrofem
