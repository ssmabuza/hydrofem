// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_Problem_Poisson.hpp"
#include "Hydrofem_Assembler_Poisson.hpp"
#include "Hydrofem_LocalStiffnessMatrix.hpp"

namespace hydrofem
{

void   
Assembler_Poisson::
buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                         const std::shared_ptr<FEVector>& res_U,
                         const std::shared_ptr<FEMatrix>& jac_U,
                         const double beta) const
{
  if (!m_jac_applied)
  {
    buildStiffMatrix(jac_U,m_basis,m_quadrature,m_dofmapper);
    m_jac_applied = true;
  }
  auto prob = std::dynamic_pointer_cast<Problem_Poisson>(m_problem);
  buildRHSVector(m_rhs,prob->rhsFunction(),m_basis,m_quadrature,m_quadrature);
  *res_U = (*jac_U)*(*U) - (*m_rhs);
}

[[maybe_unused]] 
void Assembler_Poisson::
buildStiffMatrix(const std::shared_ptr<FEMatrix>& stiff,
                 const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                 const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                 const std::shared_ptr<const DofMapper>& dofmapper)
{
  stiff->setZero();
  const auto& loc_ind = dofmapper->getLocDofIndexes();
  auto local_stiff_matrix = createKArray<LMAT_<double>>(dofmapper->local_ndof(),dofmapper->local_ndof());
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(dofmapper->nelements() * dofmapper->local_ndof() * dofmapper->local_ndof());
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    LocalStiffnessMatrix(local_stiff_matrix,dofmapper->getLocDofIndexes(),dofmapper->mesh()->getElementVertices(elem_ind),quadrature->at(elem_ind),*basis,-1.0);
    const auto& glob_ind = dofmapper->getGlobDofIndexes(elem_ind);
    for (int i = 0; i < dofmapper->local_ndof(); ++i)
      for (int j = 0; j < dofmapper->local_ndof(); ++j)
      {
        tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[i], glob_ind[j], local_stiff_matrix(loc_ind[i], loc_ind[j])));
      }
  }
  stiff->setFromTriplets(tripletList.begin(), tripletList.end());  
}

void Assembler_Poisson::
buildRHSVector(const std::shared_ptr<FEVector>& rhs,
               const std::shared_ptr<ScalarAnalyticalExpression>& f,
               const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
               const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
               const std::shared_ptr<const DofMapper>& dofmapper)
{
  const auto mesh = dofmapper->mesh();
  const auto& loc_ind = dofmapper->getLocDofIndexes();
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    const auto& glob_ind = dofmapper->getGlobDofIndexes(elem_ind);
    
    const Quadrature& quadrature_e = * quadrature->at(elem_ind);
    const auto& qweights = quadrature_e.get_q_weights();
    const auto& qpoints = quadrature_e.get_q_points();
    const auto element = mesh->getElementVertices(elem_ind);
//     const std::vector<FEBasis::FEBasisPtr>& basis_e = *basis;
    
    for (int i = 0; i < dofmapper->local_ndof(); ++i)
    {
      double loc_val = 0.0;
      for (int qp = 0; qp < quadrature_e.size(); ++qp)
        loc_val += qweights.at(qp) * f(qpoints.at(qp)) * (*(basis->at(loc_ind.at(i))))(qpoints.at(qp),element);
      (*rhs)[glob_ind[i]] += loc_val;
    }
  }
}

void Assembler_Poisson::
applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                 const std::shared_ptr<FEMatrix>& jac_U) const
{
  auto bc_ = std::dynamic_pointer_cast<BC_Scalar>(m_problem->bc());
  if (bc_)
  {
    for (auto it = bc_->beginBCPoints(); it != bc_->endBCPoints(); ++it)
    {
      if (bc_->isDirichlet(it))
      { 
        // get point index 
        const auto i = it->first;
        // boundary value = 0 here
        (*res_U)[i] = 0.0;
        if (!m_jac_applied)
        {
          // zero out the boundary row
          jac_U->row(i) *= 0.0;
          // set diagonal entry in compressed matrix to 1
          jac_u->coeffRef(i,i) = 1.0;
        }
      }
    }
  }
}

}
// end namespace hydrofem
