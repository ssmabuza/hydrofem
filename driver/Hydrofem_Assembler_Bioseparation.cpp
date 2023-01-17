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

#include "Hydrofem_Quadrature_Edge2D.hpp"

namespace hydrofem
{

void   
Assembler_Bioseparation::
buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                         const std::shared_ptr<const FEVector>& U_dot,
                         const std::shared_ptr<FEVector>& res_U,
                         const std::shared_ptr<FEMatrix>& jac_U,
                         const double /*time*/,
                         const double /*delta_t*/,
                         const double beta) const
{
  res_U->setZero();
  jac_U->setZero();

  const Problem_Bioseparation::Ptr problem = std::dynamic_pointer_cast<Problem_Bioseparation>(m_problem); assert(problem);
  const BC_Bioseparation::Ptr bc = std::dynamic_pointer_cast<BC_Bioseparation>(problem->getBoundaryCondition()); assert(bc);

  // physical problem data from input file or input script
  //@{
  const auto bc_fnc   = problem->getBoundaryFunction();
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
  const double uavg   = fr/((pi*width/2)*(pi*width/2));

  // diffusion tensor
  const double d11 = omega*d0 + alphaT*uavg;
  const double d12 = 0.0;
  const double d21 = 0.0;
  const double d22 = omega*d0 + alphaL*uavg;
  //@}

  // fluid velocity profile
  const std::function<SPoint(SPoint)> vel = [=](SPoint x)->SPoint
  { return SPoint(0.0,-3*fr*(x.x() - width)*x.x()/(4*pi*std::pow(width/2.0,3))); };

  // compute the limiter first
  if (m_do_afc) m_limiter->buildLimiter(m_nodal_limiter,U);

  std::vector<Eigen::Triplet<double>> D_tripletList;
  std::vector<Eigen::Triplet<double>> M_tripletList;
  if (m_do_afc)
  {
    if (!m_built_d_mat_graph)
      D_tripletList.reserve(m_dofmapper->nelements() * m_dofmapper->local_ndof() * m_dofmapper->local_ndof());
    M_tripletList.reserve(m_dofmapper->nelements() * m_dofmapper->local_ndof() * m_dofmapper->local_ndof());
  }
  
  for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
  {
    // local indexes
    const auto& loc_ind = m_dofmapper->getLocDofIndexes();
    // the global dof indexes in this element
    const auto& glob_ind = m_dofmapper->getGlobDofIndexes(elem_ind);
    // add boundary contribution on gamma_plus
    const auto& element = m_mesh->getElement(elem_ind);
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
      U_loc[loc_ind[i]] = (*U)[glob_ind[i]];
      U_dot_loc[loc_ind[i]] = (*U_dot)[glob_ind[i]];
    }
    
    // local mass matrix
    auto mat_M  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local lumped mass matrix
    LMAT_<double> mat_Ml;// = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local artificial diffusion matrix
    LMAT_<double> mat_D;
    if (m_do_afc)
    {
      mat_Ml = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
      if (!m_built_d_mat_graph) mat_D = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    }
    // local stiff matrix
    auto mat_S  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    // local conv matrix
    auto mat_K  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
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
    if (m_do_afc)
    {
      for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
      {
        double val = 0.0;
        for (std::size_t dof_j = 0; dof_j < loc_ind.size(); ++dof_j)
          val += mat_M(dof_i,dof_j);
        mat_Ml(dof_i,dof_i) = val;
      }

      // form triplet
      for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
        for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
          M_tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[i], glob_ind[j], mat_M(loc_ind[i], loc_ind[j])));

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

    // compute wall contributions
    LVEC_<double> res_loc_bdry = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    for (std::size_t ledge_ind = 0; ledge_ind < element.getEdges().size(); ++ledge_ind)
    {
      const auto& edge_ind = element.m_edges.at(ledge_ind);
      if (m_mesh->getEdge(edge_ind).m_is_boundary)
      {
        if (bc->isDirichletEdge(edge_ind))
        {
          Quadrature_Edge2D quad_bdry(2,m_mesh->getEdgeVertices(edge_ind));
          const auto nout = m_mesh->evalEdgeNormal(elem_ind,ledge_ind);
          const auto& q_pts_b = quad_bdry.get_q_points();
          const auto& q_wts_b = quad_bdry.get_q_weights();
          for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
            for (int quad_ind = 0; quad_ind < quad_bdry.size(); ++quad_ind)
            {
              const int& loc_i = loc_ind.at(dof_i);
              double phi_i = (*fe_basis_e.at(loc_i))(q_pts_b.at(quad_ind),elem_verts);
              double v_dot_n = vel(q_pts_b.at(quad_ind)) * nout;
              double cin = (*bc_fnc)(q_pts_b.at(quad_ind));
              res_loc_bdry[dof_i] += q_wts_b.at(quad_ind) * phi_i * v_dot_n * cin;
            }

        } else if (bc->isGammaPlusEdge(edge_ind)) {

          Quadrature_Edge2D quad_bdry(2,m_mesh->getEdgeVertices(edge_ind));
          const auto nout = m_mesh->evalEdgeNormal(elem_ind,ledge_ind);
          const auto& q_pts_b = quad_bdry.get_q_points();
          const auto& q_wts_b = quad_bdry.get_q_weights();
          for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j <= dof_i; ++dof_j)
              for (int quad_ind = 0; quad_ind < quad_bdry.size(); ++quad_ind)
              {
                const int& loc_i = loc_ind.at(dof_i);
                const int& loc_j = loc_ind.at(dof_j);
                double phi_i = (*fe_basis_e.at(loc_i))(q_pts_b.at(quad_ind),elem_verts);
                double phi_j = (*fe_basis_e.at(loc_j))(q_pts_b.at(quad_ind),elem_verts);
                double v_dot_n = vel(q_pts_b.at(quad_ind)) * nout;
                res_loc_bdry[dof_i] += q_wts_b.at(quad_ind) * phi_i * v_dot_n * phi_j * U_loc[loc_j];
              }

        }
      }
    }

    // build the problem convection matrix and local conv operator
    std::vector<LMAT_<double>> c_loc(m_mesh->numOfDims());
    if (m_do_afc && !m_built_d_mat_graph)
      for (auto& c_loc_i : c_loc) c_loc_i = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
    for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
    {
      const int& loc_i = loc_ind.at(dof_i);
      for (std::size_t dof_j = 0; dof_j < loc_ind.size(); ++dof_j)
      {
        const int& loc_j = loc_ind.at(dof_j);
        double val = 0.0;
        for (int quad_ind = 0; quad_ind < quad_e.size(); ++quad_ind)
        {
          const auto dphi_i = (*fe_basis_e.at(loc_i)).grad(q_points.at(quad_ind),elem_verts);
          const double phi_j = (*fe_basis_e.at(loc_j))(q_points.at(quad_ind),elem_verts);
          const auto velh = vel(q_points.at(quad_ind));
          val += - q_weights.at(quad_ind) * dphi_i * velh * phi_j;

          if (m_do_afc && !m_built_d_mat_graph)
          {
            const auto dphi_j = (*fe_basis_e.at(loc_j)).grad(q_points.at(quad_ind),elem_verts).data();
            const double phi_i = (*fe_basis_e.at(loc_i))(q_points.at(quad_ind),elem_verts);
            for (std::size_t dim = 0; dim < c_loc.size(); ++dim)
              c_loc.at(dim)(loc_i,loc_j) += phi_i * dphi_j(dim);
          }
        }
        mat_K(loc_i,loc_j) = val;
      }
    }

    // compute the local artificial diffusion matrix
    if (m_do_afc && !m_built_d_mat_graph)
    {

      for (std::size_t dof_i = 0; dof_i < loc_ind.size(); ++dof_i)
      {
        const int& loc_i = loc_ind.at(dof_i);
        const int& glob_i = element.m_nodes[loc_i];
        const auto vel_i = vel(m_mesh->getPoint(glob_i)).data();
        for (std::size_t dof_j = 0; dof_j <= dof_i; ++dof_j)
        {
          const int& loc_j = loc_ind.at(dof_j);
          if (loc_i!=loc_j)
          {
            const int& glob_j = element.m_nodes[loc_j];
            const auto vel_j = vel(m_mesh->getPoint(glob_j)).data();

            std::vector<double> prod(4);
            for (std::size_t dim = 0; dim < c_loc.size(); ++dim)
            {
              prod[0] += c_loc.at(dim)(loc_i,loc_j) * vel_i(dim);
              prod[1] += c_loc.at(dim)(loc_i,loc_j) * vel_j(dim);
              prod[2] += c_loc.at(dim)(loc_j,loc_i) * vel_i(dim);
              prod[3] += c_loc.at(dim)(loc_j,loc_i) * vel_j(dim);
            }
            const double dij = *(std::max_element(prod.begin(),prod.end()));
            mat_D(loc_i,loc_j) = dij;
            mat_D(loc_j,loc_i) = dij;
            mat_D(loc_i,loc_i) -= dij;
            mat_D(loc_j,loc_j) -= dij;
          }
        }
      }

      // form triplet
      for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
        for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
          D_tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[i], glob_ind[j], mat_D(loc_ind[i], loc_ind[j])));

    }
    
    // compute the local residual contributions
    for (int i = 0; i < mat_M.dimension(0); ++i)
      for (int j = 0; j < mat_M.dimension(1); ++j)
      {
        if (m_do_afc)
          res_loc(i) += mat_Ml(i,j) * U_dot_loc(j) + (mat_S(i,j)+mat_K(i,j)+mat_D(i,j))*U_loc(j);
        else
          res_loc(i) += mat_M(i,j) * U_dot_loc(j) + (mat_S(i,j)+mat_K(i,j))*U_loc(j);
      }
    
    for (int i = 0; i < mat_M.dimension(0); ++i)
      res_loc(i) += res_loc_bdry[i];

    // compute the local jacobian contributions
    if (m_do_afc)
      mat_J = beta * mat_Ml + mat_S + mat_K;
    else
      mat_J = beta * mat_M + mat_S + mat_K;

    // Jacobian triples set up
    for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
    {
      (*res_U)[glob_ind[i]] += res_loc[loc_ind[i]];
      for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
        jac_U->coeffRef(glob_ind[i], glob_ind[j]) += mat_J(loc_ind[i], loc_ind[j]);
        //jac_tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[i], glob_ind[j], mat_J(loc_ind[i], loc_ind[j])));
    }

  }
  // assemble Jacobian
  //jac_U->setFromTriplets(jac_tripletList.begin(), jac_tripletList.end());
  if (m_do_afc && !m_built_d_mat_graph)
    m_D_mat->setFromTriplets(D_tripletList.begin(), D_tripletList.end());
  if (m_do_afc)
    m_M_mat->setFromTriplets(M_tripletList.begin(), M_tripletList.end());

  // Do limiting and add to residual & Jacobian
  if (m_do_afc)
  {
    // build limited internodal fluxes
    for (int edge_ind = 0; edge_ind < m_mesh->numOfEdges(); ++edge_ind)
    {
      const auto& edge = m_mesh->getEdge(edge_ind);
      const int& i = edge.m_nodes[0];
      const int& j = edge.m_nodes[1];

      const double mij = m_M_mat->coeff(i,j);
      const double dij = m_D_mat->coeff(i,j);
      
      const double& Ui = (*U)[i];
      const double& Uj = (*U)[j];
      const double& dot_Ui = (*U_dot)[i];
      const double& dot_Uj = (*U_dot)[j];

      // compute the edge limiter (average for now)
      double alpha_ij = 0.5*((*m_nodal_limiter)[i] + (*m_nodal_limiter)[j]);

      double bar_fij = - alpha_ij * mij * (dot_Ui - dot_Uj) - alpha_ij * dij * (Ui - Uj);
      (*m_afc_vec)[i] += bar_fij;
      (*m_afc_vec)[j] -= bar_fij;

      const double aij = beta * alpha_ij * mij + alpha_ij * dij;
      jac_U->coeffRef(i,j) += aij;
      jac_U->coeffRef(i,i) -= aij;
      jac_U->coeffRef(j,j) -= aij;
    }

    // close the sparse matrix
    //jac_U->makeCompressed();

    // check if the AFC contribution is conservative
    if (m_afc_vec->sum() > point_eps)
    {
      throw std::runtime_error("AFC contribution in \"Assembler_Bioseparation\" is not conservative.");
    }

    // sum into total residual
    (*res_U) += (*m_afc_vec);
  }

  // std::cout << "U = \n" << *U << std::endl;
  // std::cout << "U_dot = \n" << *U_dot << std::endl;
  // std::cout << "residual = \n" << *res_U << std::endl;
  // std::cin.get();
  
}

void Assembler_Bioseparation::
applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                 const std::shared_ptr<FEMatrix>& jac_U) const
{
  const auto bc_ = std::dynamic_pointer_cast<BC_Scalar>(m_problem->getBoundaryCondition());
  const auto bc_fnc = m_problem->getBoundaryFunction();
  assert(bc_fnc);
  if (bc_)
  {
    for (auto it = bc_->startBCPoints(); it != bc_->endBCPoints(); ++it)
    {
      if (bc_->isDirichlet(it))
      { 
        // get point index 
        const auto i = it->first;
        // Dirichlet boundary value 
        (*res_U)[i] = 0.0;//(*bc_fnc)(m_mesh->getPoint(i));
        // zero out the boundary row
        jac_U->row(i) *= 0.0;
        // set diagonal entry in compressed matrix to 1
        jac_U->coeffRef(i,i) = 1.0;
      }
    }
  }
  // compress the matrix if modified
  jac_U->makeCompressed();

  // print out the matrix
//  std::cout << *res_U << std::endl;

}

}
// end namespace hydrofem
