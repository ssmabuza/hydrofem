// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_AFC_LPSLimiter.hpp"

namespace hydrofem
{

AFC_LPSLimiter::
AFC_LPSLimiter(const std::shared_ptr<Mesh>& mesh,
               const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
               const std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>>& basis)
  :
  m_eps(1.0e-16 * mesh->getMeshSize())
{
  m_mesh = mesh;
  m_quadrature = quadrature;
  m_basis = basis;
  m_P = std::make_shared<FEVector>();
  m_Q = std::make_shared<FEVector>();
}

void
AFC_LPSLimiter::
buildLimiter(const std::shared_ptr<FEVector>& res,
             const std::shared_ptr<const FEVector>& u) const
{
  auto& res_view = *res;
  const auto& u_view = *u;
  
  m_P->setZero();
  auto& m_P_view = *m_P;
  m_Q->setZero();
  auto& m_Q_view = *m_Q;

  for (int elem_ind = 0; elem_ind < m_mesh->numOfElements(); ++elem_ind)
  {
    const Quadrature& quadrature_e = * m_quadrature->at(elem_ind);
    const auto & qweights = quadrature_e.get_q_weights();
    const auto & qpoints = quadrature_e.get_q_points();
    
    const std::vector<std::shared_ptr<FEBasis>>& basis_e = m_basis->at(elem_ind);
    
    // build the average of u on K_e 
    double u_e = 0.0;
    const auto& point_ind = m_mesh->getElementGlobalIndexes(elem_ind);
    const auto element = m_mesh->getElementVertices(elem_ind);
    const auto& loc_ind = point_ind;
    
    for (std::size_t i = 0; i < point_ind.size(); ++i)
      u_e += u_view[loc_ind[i]];
    u_e /= static_cast<double>(point_ind.size());
    
    LVEC_<double> lps_value(point_ind.size());
    for (std::size_t i = 0; i < point_ind.size(); ++i)
      lps_value[i] = 0.0;
    
    // integrate to get the monotone local projection of the solution in K_e
    for (int qp = 0; qp < quadrature_e.size(); ++qp)
    {
      // build the solution value at qp
      double uh = 0.0;
      for (std::size_t i = 0; i < point_ind.size(); ++i)
        uh += u_view[loc_ind[i]] * basis_e[i]->operator()(qpoints[qp],element);
      
      for (std::size_t i = 0; i < point_ind.size(); ++i)
        lps_value[i] += qweights[qp] * (uh - u_e) * basis_e[i]->operator()(qpoints[qp],element);
    }
    
    for (std::size_t i = 0; i < point_ind.size(); ++i)
    {
      m_P_view[loc_ind[i]] += lps_value[i];
      m_Q_view[loc_ind[i]] += std::fabs(lps_value[i]);
    }
  }
  
  // compute the final shock detector 
  for (long int i = 0; i < res_view.size(); ++i)
    res_view[i] = (std::fabs(m_P_view[i])+m_eps)/(m_Q_view[i]+m_eps);
}

}
// end namespace valiant
