// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "InitialSolution.hpp"
#include "GlobalMassMatrix.hpp"
#include "LinearSolvers_Tpetra.hpp"

#include "InitialSolution_HyperbolicSystem.hpp"

namespace valiant
{

void InitialSolution_HyperbolicSystem::evaluate()
{
  if (m_proj_type == "Lumped") {
    buildLumpedMassMatrix();
    evaluateLumpedMassProjection();
  } else if (m_proj_type == "Consistent") {
    buildConsistentMassMatrix();
    evaluateConsistentProjection();
  } else
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"The initial projection provided : \""
                               << m_proj_type << "\" is either not implemented or invalid.");
}

void
InitialSolution_HyperbolicSystem::
buildLumpedMassMatrix()
{
  m_lumped_mass = Teuchos::rcp(new TVector(m_lob->getOwnedMap()));
  GlobalLumpedMassMatrix(m_lumped_mass,m_dofmapper,m_fe_basis,m_quadrature);
}

void
InitialSolution_HyperbolicSystem::
buildConsistentMassMatrix()
{
  // create the mass matrix from an existing Crs graph
  m_mass = Teuchos::rcp(new TMatrix(m_lob->getOwnedGraph()));
  // assemble the matrix
  GlobalMassMatrix(m_mass,m_dofmapper,m_fe_basis,m_quadrature);
}

void
InitialSolution_HyperbolicSystem::
evaluateLumpedMassProjection()
{
  // typedefs
  using LVec = LVEC_<double>;
  // local view lumped mass matrix
  const auto lumped_mass_view = m_lumped_mass->getData();
  
  std::vector<TVector> initial_values_vecs(m_neq,TVector(m_lob->getOwnedMap()));
  
  const auto& loc_ind = m_dofmapper->getLocDofIndexes();
  
  // loop over the local elements
  for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
  {
    const auto& quadrature_e = * m_quadrature->at(elem_ind);
    const auto& fe_bases_e = *m_fe_basis;
    
    const auto& quad_p = quadrature_e.get_q_points();
    const auto& quad_w = quadrature_e.get_q_weights();
    
    const auto& glob_ind = m_dofmapper->getGlobDofIndexes(elem_ind);
    const auto element = m_dofmapper->mesh()->getElementVertices(elem_ind);
    
    for (int qp = 0; qp < quadrature_e.size(); ++qp)
    {
      const LVec u = m_ic->evaluate(quad_p.at(qp));
      for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
      {
        for (int eq = 0; eq < m_neq; ++eq)
        {
          const double value = quad_w.at(qp) * (*(fe_bases_e[loc_ind[j]]))(quad_p.at(qp),element) * u(eq) / lumped_mass_view[glob_ind[j]];
          initial_values_vecs[eq].sumIntoGlobalValue(glob_ind[j],value);
        }
      }
    }
  }
  
  for (int eq = 0; eq < m_neq; ++eq)
  {
    m_scalar_gather->doGather(m_result_gathered->at(eq),Teuchos::rcpFromRef(initial_values_vecs[eq]));
  }
  
  m_global_scatter->doScatter(*m_result_gathered,m_result);
  m_is_computed = true;
}

void
InitialSolution_HyperbolicSystem::
evaluateConsistentProjection()
{
  // typedefs
  using LVec = LVEC_<double>;
  // initial values to be computed
  std::vector<TVector> initial_values_vecs(m_neq,TVector(m_lob->getOwnedMap()));
  // loop over the local elements
  for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
  {
    const auto& quadrature_e = * m_quadrature->at(elem_ind);
    const auto& fe_bases_e = *m_fe_basis;
    
    const auto& quad_p = quadrature_e.get_q_points();
    const auto& quad_w = quadrature_e.get_q_weights();
    
    const auto& glob_ind = m_dofmapper->getGlobDofIndexes(elem_ind);
    const auto element = m_dofmapper->mesh()->getElementVertices(elem_ind);
    
    for (int qp = 0; qp < quadrature_e.size(); ++qp)
    {
      LVec u = m_ic->evaluate(quad_p.at(qp));
      for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
      {
        for (int eq = 0; eq < m_neq; ++eq)
        {
          initial_values_vecs[eq].sumIntoGlobalValue(glob_ind[j], quad_w.at(qp)*(*(fe_bases_e[j]))(quad_p.at(qp),element) * u(eq));
        }
      }
    }
  }
  
  // do the initial projection for each variable
  for (int eq = 0; eq < m_neq; ++eq)
  {
    auto initial_values_res = Teuchos::rcp(new TVector(initial_values_vecs[eq].getMap()));
    initial_values_res->putScalar(0.0);
    LinearSolvers_Tpetra::solveSystemCG(m_mass,Teuchos::rcpFromRef(initial_values_vecs[eq]),initial_values_res);
    m_scalar_gather->doGather(m_result_gathered->at(eq),initial_values_res);
  }
  // finally scatter the info to a system vector
  m_global_scatter->doScatter(*m_result_gathered,m_result);
  m_is_computed = true;
}

}
// end namespace valiant
