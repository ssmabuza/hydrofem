// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_InitialSolution.hpp"

#include "Hydrofem_LinearSolvers.hpp"
#include "Hydrofem_GlobalMassMatrix.hpp"

namespace hydrofem
{

std::shared_ptr<FEVector>
InitialSolution::
get_evaluatedField() const
{
  assert(m_is_computed);
  return m_result;
}

void NodalProjection::
evaluate()
{
  const auto mesh = m_lob->dofMapper()->mesh();
  for (Eigen::Index i = 0; i < m_result->size(); ++i)
    (*m_result)[i] = m_ic->evaluate(mesh->getPoint(i));
  m_is_computed = true;
}

ConsistentMassProjection::
ConsistentMassProjection(const std::shared_ptr<LinearObjectBuilder>& lob,
                         const std::shared_ptr<ScalarInitialCondition>& ic_fnc,
                         const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                         const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature) : InitialSolution(lob,ic_fnc)
{
  m_result = std::make_shared<FEVector>(m_lob->dofMapper()->mesh()->numOfPoints());
  m_mass = m_lob->createSparseMatrix();
  m_basis = fe_basis;
  m_quadrature = quadrature;
  GlobalMassMatrix(m_mass,m_lob->dofMapper(),m_basis,m_quadrature);
  m_mass_matrix_built = true;
}

void ConsistentMassProjection::evaluate()
{
  const auto dofmapper = m_lob->dofMapper();
  auto rhs = m_lob->createVector(); //rhs->resize(dofmapper->global_ndof());
  
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    const auto& loc_inds = dofmapper->getLocDofIndexes();
    const auto& glob_inds = dofmapper->getGlobDofIndexes(elem_ind);
    const auto elem_verts = dofmapper->mesh()->getElementVertices(elem_ind);

    const auto& quad_e = *(m_quadrature->at(elem_ind));
    const auto& q_pts = quad_e.get_q_points();
    const auto& q_wts = quad_e.get_q_weights();

    for (std::size_t i = 0; i < loc_inds.size(); ++i)
    {
      const auto& phi_i = *(m_basis->at(loc_inds[i]));
      double val (0.);
      for (int qp = 0; qp < quad_e.size(); ++qp)
      {
        val += q_wts[qp] * phi_i(q_pts[qp],elem_verts) * m_ic->evaluate(q_pts[qp]);
      }
      (*rhs)[glob_inds[i]] += val;
    }
  }

  // build the linear solver 
  LinearSolverInterface lsi("cg");
  // do the projection
  lsi.solve(m_mass,rhs,m_result);
  m_is_computed = true;
}

LumpedMassProjection::
LumpedMassProjection(const std::shared_ptr<LinearObjectBuilder>& lob,
                         const std::shared_ptr<ScalarInitialCondition>& ic_fnc,
                         const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                         const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature) : InitialSolution(lob,ic_fnc)
{
  m_result = std::make_shared<FEVector>(m_lob->dofMapper()->mesh()->numOfPoints());
  m_lumped_mass = m_lob->createVector();
  m_basis = fe_basis;
  m_quadrature = quadrature;
  GlobalLumpedMassMatrix(m_lumped_mass,m_lob->dofMapper(),m_basis,m_quadrature);
  m_lumped_mass_matrix_built = true;
}

void LumpedMassProjection::evaluate()
{
  const auto dofmapper = m_lob->dofMapper();
  auto rhs = m_lob->createVector(); rhs->resize(dofmapper->global_ndof());
  
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    const auto& loc_inds = dofmapper->getLocDofIndexes();
    const auto& glob_inds = dofmapper->getGlobDofIndexes(elem_ind);
    const auto elem_verts = dofmapper->mesh()->getElementVertices(elem_ind);

    const auto& quad_e = *(m_quadrature->at(elem_ind));
    const auto& q_pts = quad_e.get_q_points();
    const auto& q_wts = quad_e.get_q_weights();

    for (std::size_t i = 0; i < loc_inds.size(); ++i)
    {
      const auto& phi_i = *(m_basis->at(loc_inds[i]));
      double val (0.);
      for (int qp = 0; qp < quad_e.size(); ++qp)
      {
        val += q_wts[qp] * phi_i(q_pts[qp],elem_verts) * m_ic->evaluate(q_pts[qp]);
      }
      (*rhs)[glob_inds[i]] += val/(*m_lumped_mass)[glob_inds[i]];
    }
  }
  m_is_computed = true;
}

}
// end namespace hydrofem
