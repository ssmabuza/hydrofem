// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Assembler_Bioseparation.hpp"

#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_BC_Bioseparation.hpp"
#include "Hydrofem_Problem_Bioseparation.hpp"
#include "Hydrofem_LocalStiffnessMatrix.hpp"
#include "Hydrofem_AFC_LocalDiffusionMatrix.hpp"

namespace hydrofem
{

void   
Assembler_Bioseparation::
buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                         const std::shared_ptr<const FEVector>& U_dot,
                         const std::shared_ptr<FEVector>& res_U,
                         const std::shared_ptr<FEMatrix>& jac_U,
                         const double /*time*/,
                         const double delta_t,
                         const double /*beta*/) const
{
  Problem_Bioseparation::Ptr problem = std::dynamic_pointer_cast<Problem_Bioseparation>(m_problem); assert(problem);
  BC_Bioseparation::Ptr bc = std::dynamic_pointer_cast<BC_Bioseparation>(problem->getBoundaryCondition()); assert(bc);

  // physical problem data from input file or input script
  //@{
  const auto vel = bc->getFluidVelocity();

  const double pi     = M_PI;
  const double omega  = problem->omega();
  const double rho_s  = problem->rho_s();
  const double q_max  = problem->q_max();
  const double Keq    = problem->Keq();
  const double width  = problem->width();
  const double fr     = problem->flowrate();
  const double alphaL = problem->alphaL();
  const double alphaT = problem->alphaT();
  const double d0     = problem->d0();

  // diffusive flow speed
  double uavg = fr/((pi*width/2)*(pi*width/2));

  // diffusion tensor
  double d11, d12, d21, d22;
  d11=omega*d0+alphaT*uavg;
  d12=0.0;
  d21=0.0;
  d22=omega*d0+alphaL*uavg;
  //@}

  // compute the limiter first
  m_limiter->buildLimiter(m_nodal_limiter,U);

  // triplet list for Jacobian build
  // TODO: prebuild Jac sparse graph and do sums
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(m_dofmapper->nelements() * m_dofmapper->local_ndof() * m_dofmapper->local_ndof());

  for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
  {
    // local indexes
    const auto& loc_ind = m_dofmapper->getLocDofIndexes();
    // the global dof indexes in this element
    const auto& glob_ind = m_dofmapper->getGlobDofIndexes(elem_ind);

    // compute the element limiter
    double alpha_e;
    // compute alpha_f_e
    {
      alpha_e = std::numeric_limits<double>::max();
      for (int glob_i : glob_ind)
        alpha_e = std::min(alpha_e,(*m_nodal_limiter)[glob_i]);
    }

    // add boundary contribution on gamma_plus
//     const auto& element = m_mesh->getElement(elem_ind);
    // element vertices
    const auto elem_verts = m_mesh->getElementVertices(elem_ind);
    // local solution vector
    auto U_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    // local solution time derivative vector
    auto U_dot_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    // local residual vector
    auto res_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    // do a deep copy of values from global to local
    for (std::size_t i = 0; i < glob_ind.size(); ++i)
    {
      U_loc[i] = (*U)[glob_ind[i]];
      U_dot_loc[i] = (*U_dot)[glob_ind[i]];
    }
    // local mass matrix
    auto mat_M  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local lumped mass matrix
    auto mat_Ml = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local stiff matrix
    auto mat_S  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local conv matrix
    auto mat_K  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local artificial diffusion matrix
    auto mat_D  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local jacobian
    auto mat_J  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // get the element quadrature
    const auto& quad_e = *(m_quadrature->at(elem_ind));
    // get the quadrature points
    const auto& q_points = quad_e.get_q_points();
    // get the quadrature weights
    const auto& q_weights = quad_e.get_q_weights();
    // get the basis functions
    const auto& fe_basis_e = *m_basis;

    // build the system adsorption mass matrix
    for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
    {
      const int& loc_i = loc_ind.at(dof_i);
      for (std::size_t dof_j = 0; dof_j <= dof_i; ++dof_j)
      {
        const int& loc_j = loc_ind.at(dof_j);
        double val = 0.0;
        for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
        {
          double Uh = 0.0;
          for (std::size_t uh_i = 0; uh_i < loc_ind.size(); ++uh_i)
            Uh += U_loc[uh_i] * (*fe_basis_e.at(uh_i))(q_points.at(quad_ind),elem_verts);

          double phi_i = (*fe_basis_e.at(loc_i))(q_points.at(quad_ind),elem_verts);
          double phi_j = (*fe_basis_e.at(loc_j))(q_points.at(quad_ind),elem_verts);
          double factor = omega + (1.0 - omega) * rho_s * q_max * Keq /( (1.0 + Keq*Uh)*(1.0 + Keq*Uh) );

          val += q_weights.at(quad_ind) * phi_i * factor * phi_j;
        }
        mat_M(loc_i,loc_j) = val;
        if (loc_i!=loc_j) mat_M(loc_j,loc_i) = val;
      }
    }

    // lump the local mass matrix
    for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
    {
      double val = 0.0;
      for (std::size_t dof_j = 0; dof_j < loc_ind.size(); ++dof_j)
        val += mat_M(dof_i,dof_j);
      mat_Ml(dof_i,dof_i) = val;
    }

    // build the problem stiffness matrix
    for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
    {
      const int& loc_i = loc_ind.at(dof_i);
      for (std::size_t dof_j = 0; dof_j <= dof_i; ++dof_j)
      {
        const int& loc_j = loc_ind.at(dof_j);
        double val = 0.0;

        for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
        {
          auto dphi_i = (*fe_basis_e.at(loc_i)).grad(q_points.at(quad_ind),elem_verts);
          auto dphi_j = (*fe_basis_e.at(loc_j)).grad(q_points.at(quad_ind),elem_verts);
          val += q_weights.at(quad_ind) *(d11 * dphi_i.x() * dphi_j.x() + d12 * dphi_i.x() * dphi_j.y()
                                        + d21 * dphi_i.y() * dphi_j.x() + d22 * dphi_i.y() * dphi_j.y());
        }
        mat_S(loc_i,loc_j) = val;
        if (loc_i!=loc_j) mat_S(loc_j,loc_i) = val;
      }
    }

    // build the problem convection matrix
    for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
    {
      const int& loc_i = loc_ind.at(dof_i);
      for (std::size_t dof_j = 0; dof_j < loc_ind.size(); ++dof_j)
      {
        const int& loc_j = loc_ind.at(dof_j);
        double val = 0.0;
        for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
        {
          auto dphi_i = (*fe_basis_e.at(loc_i)).grad(q_points.at(quad_ind),elem_verts);
          double phi_j = (*fe_basis_e.at(loc_j))(q_points.at(quad_ind),elem_verts);
          auto velh = (*vel)(q_points.at(quad_ind));
          val += - q_weights.at(quad_ind) * dphi_i * velh * phi_j;
        }
        mat_K(loc_i,loc_j) = val;
      }
    }

    // compute the local artificial diffusion matrix
    AFC_LocalDiffusionMatrix(mat_D, mat_K);

    // compute the local residual contributions
    res_loc = mat_M * U_dot_loc + mat_S * U_loc + mat_K * U_loc;
    // compute the local jacobian contributions
    mat_J = (1.0/delta_t)*mat_M + mat_S + mat_K;

    // Jacobian triples set up
    for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
    {
      (*res_U)[glob_ind[i]] = res_loc[loc_ind[i]];
      for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
        tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[i], glob_ind[j], mat_J(loc_ind[i], loc_ind[j])));
    }

  }
  // assemble Jacobian
  jac_U->setFromTriplets(tripletList.begin(), tripletList.end());

}

void Assembler_Bioseparation::
applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                 const std::shared_ptr<FEMatrix>& jac_U) const
{
  const auto bc_ = std::dynamic_pointer_cast<BC_Scalar>(m_problem->getBoundaryCondition());
  if (bc_)
  {
    for (auto it = bc_->startBCPoints(); it != bc_->endBCPoints(); ++it)
    {
      if (bc_->isDirichlet(it))
      { 
        // get point index 
        const auto i = it->first;
        // boundary value = 0 here
        (*res_U)[i] = 0.0;
        // zero out the boundary row
        jac_U->row(i) *= 0.0;
        // set diagonal entry in compressed matrix to 1
        jac_U->coeffRef(i,i) = 1.0;
      }
    }
  }
}

}
// end namespace hydrofem
