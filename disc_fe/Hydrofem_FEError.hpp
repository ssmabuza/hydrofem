// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEError_HPP__
#define __Hydrofem_FEError_HPP__

#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

class FEError
{
public:
  
  FEError(const std::shared_ptr<DofMapper>& dofmapper,
          const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
          const std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>>& basis)
  {
    m_dofmapper = dofmapper;
    m_quadrature = quadrature;
    m_basis = basis;
  }
  
  virtual ~FEError() {}
  
  /** @brief L^n norm  */
  template <int n>
  inline double compute(const FEArray<double>::CellBasis& u, const std::function<double(SPoint)>& ref_sol)
  {
    double res = 0.0;
    for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
    {
      const Quadrature& quad_e = * (*m_quadrature)[elem_ind];
      const auto& q_points = quad_e.get_q_points();
      const auto& q_weights = quad_e.get_q_weights();
      const std::vector<std::shared_ptr<FEBasis>>& fe_basis_e = (*m_basis)[elem_ind];
      const std::vector<int>& loc_ind = m_dofmapper->getLocDofIndexes();
      
      for (int qp = 0; qp < quad_e.size(); ++qp)
      {
        double uh = 0.0;
        for (std::size_t b = 0; b < loc_ind.size(); ++b)
          uh += u(elem_ind,b) * (*fe_basis_e[b])(q_points[qp]);
        res += q_weights[qp] * std::pow(std::fabs(uh - ref_sol(q_points[qp])),double(n));
      }
    }
    return std::pow(res,1.0/double(n));
  }
  
protected:
  
  std::shared_ptr<DofMapper> m_dofmapper;
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
  std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>> m_basis;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_FEError_HPP__ */
