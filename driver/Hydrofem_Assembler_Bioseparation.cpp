// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_BC_Bioseparation.hpp"
#include "Hydrofem_Problem_Bioseparation.hpp"
#include "Hydrofem_Assembler_Bioseparation.hpp"
#include "Hydrofem_LocalStiffnessMatrix.hpp"

namespace hydrofem
{

void   
Assembler_Bioseparation::
buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                         const std::shared_ptr<const FEVector>& U_dot,
                         const std::shared_ptr<FEVector>& res_U,
                         const std::shared_ptr<FEMatrix>& jac_U,
                         const double time,
                         const double delta_t,
                         const double beta) const
{
  // available:
  // m_dofmapper (mesh)
  // basis_i
  // quadrature
  //

  // get the constants from the problem:
  // dynamic cast of the problem
  const auto problem = std::dynamic_pointer_cast<Problem_Bioseparation>(m_problem);
  assert(problem);
  const double omega = problem->omega();
  const double rho_s = problem->rho_s();
  const double q_max = problem->q_max();
  const double Keq = problem->Keq();

  const auto bc = std::dynamic_pointer_cast<BC_Bioseparation>(problem->getBC());
  assert(bc);
  const auto vel = bc->getFluidVelocity();

  double d11, d12, d21, d22;
  const double width = problem->width();
  const double fr = problem->flowrate();
  const double pi = M_PI;
  double uavg = fr/((pi*width/2)*(pi*width/2));
  // TODO : move to program options
  double alphaL = 1.37;
  double alphaT = 0.137;
  double d0 = 0.0000228;

  d11=omega*d0+alphaT*uavg;
  d12=0.0;
  d21=0.0;
  d22=omega*d0+alphaL*uavg;


  for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
  {

    auto U_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    auto U_dot_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    auto res_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());

    // the global dof indexes in this element
    const auto& glob_inds = m_dofmapper->mesh()->getElementGlobalIndexes(elem_ind);

    for (std::size_t i = 0; i < glob_inds.size(); ++i)
    {
      U_loc[i] = (*U)[glob_inds[i]];
      U_dot_loc[i] = (*U_dot)[glob_inds[i]];
    }

    auto mat_M  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());

    // local indexes
    const auto& loc_inds = m_dofmapper->getLocDofIndexes();

    const auto& quad_e = *(m_quadrature->at(elem_ind));
    const auto& q_points = quad_e.get_q_points();
    const auto& q_weights = quad_e.get_q_weights();
    const auto& fe_basis_e = *m_basis;

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
            Uh += U_loc[uh_i] * (*fe_basis_e.at(uh_i))(q_points.at(quad_ind),element);

          double phi_i = (*fe_basis_e.at(loc_i))(q_points.at(quad_ind),element);
          double phi_j = (*fe_basis_e.at(loc_j))(q_points.at(quad_ind),element);
          double factor = omega + (1.0 - omega) * rho_s * q_max * Keq /( (1.0 + Keq*Uh)*(1.0 + Keq*Uh) );

          val += q_weights.at(quad_ind) * phi_i * factor * phi_j;
        }
        mat_M(loc_i,loc_j) = val;
        if (loc_i!=loc_j)
          mat_M(loc_j,loc_i) = val;
      }
    }




  }






//   // compute the limiter from the solution values
//   m_limiter->buildLimiter(m_nodal_limiter,U);
//   // get the mesh
//   const auto mesh = m_dofmapper->mesh();
//   // set all the values in the residual to zero
//   res_U->setZero();
//   // set all the values in the Jacobian to zero
//   jac_U->setZero();
//
//   // pre-allocate some local matrices
//   auto mat_M  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
//   auto mat_Ml = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
//   auto mat_K  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
//   auto mat_D  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
//   auto mat_S  = createKArray<LMAT_<double>>(m_dofmapper->local_ndof(),m_dofmapper->local_ndof());
//
//   auto mat_sigma = createKArray<LMAT_<double>>(4,4);
//
//   auto U_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
//   auto U_dot_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
//   auto res_loc = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
//
//   // element limiters
//   double alpha_f (0.0);
//
//   // local indexes
//   const auto& loc_inds = m_dofmapper->getLocDofIndexes();
//
//   auto bc = std::dynamic_pointer_cast<BC_Scalar>(m_problem->bc());
//
//   // get BC info in the edges
//   const auto& bc_info_edges = bc->bcEdges();
//
//   // get BC info on the nodes
//   const auto& bc_info_pts = bc->bcPoints();
//
//   // assemble residual
//   for (int elem_ind = 0; elem_ind < mesh->numOfElements(); ++elem_ind)
//   {
//     // the global dof indexes in this element
//     const auto& glob_inds = mesh->getElementGlobalIndexes(elem_ind);
//     // add boundary contribution on gamma_plus
//     const auto& element = mesh->getElement(elem_ind);
//     // element vertices
//     const auto elemVertices = mesh->getElementVertices(elem_ind);
//
//     // zero out local vectors
//     zeroOutArray(U_loc);
//     zeroOutArray(U_dot_loc);
//     zeroOutArray(res_loc);
//
//     // zero out local matrices
//     zeroOutArray(mat_M);
//     zeroOutArray(mat_Ml);
//     zeroOutArray(mat_K);
//     zeroOutArray(mat_D);
//     zeroOutArray(mat_S);
//
//     for (std::size_t i = 0; i < glob_inds.size(); ++i)
//     {
//       U_loc[i] = (*U)[glob_inds[i]];
//       U_dot_loc[i] = (*U_dot)[glob_inds[i]];
//     }
//
//     // compute alpha_f_e
//     {
//       alpha_f = std::numeric_limits<double>::max();
//       for (int glob_ind : glob_inds)
//         alpha_f = std::min(alpha_f,(*m_nodal_limiter)[glob_ind]);
//     }
//
//     // compute the local conv matrix
//     LocalConvectionMatrix(mat_K,
//                           loc_inds,
//                           elemVertices,
//                           m_quadrature_f->at(elem_ind),
//                           *m_basis_f,
//                           *m_velocity_fnc);
//
//     // compute the local stiffness matrix
//     LocalStiffnessMatrix(mat_S,
//                          loc_inds,
//                          elemVertices,
//                          m_quadrature_f->at(elem_ind),
//                          *m_basis_f,
//                          m_diff);
//
//     LocalMassMatrix(mat_M,
//                     loc_inds,
//                     elemVertices,
//                     m_quadrature_f->at(elem_ind),
//                     *m_basis_f);
//
//     LocalLumpedMassMatrix(mat_Ml,
//                           loc_inds,
//                           elemVertices,
//                           m_quadrature_f->at(elem_ind),
//                           *m_basis_f);
//
//     // do calculations on the boundary plus limiting
//     for (std::size_t ledge_ind = 0; ledge_ind < element.getEdges().size(); ++ledge_ind)
//     {
//       const int& edge_ind = element.m_edges[ledge_ind];
//       if (bc_info_edges[edge_ind].m_boundary_condition_type == typeGammaPlus)
//       {
//         // get quadrature on this edge
//         const auto& quadrature = m_quadrature_bdry->at(edge_ind);
//         // integrate
//         const auto& q_points = quadrature->get_q_points();
//         const auto& q_weights = quadrature->get_q_weights();
//         const auto normal = mesh_f->evalEdgeNormal(elem_ind,static_cast<int>(ledge_ind));
//         for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
//         {
//           const auto& basis_i = *(m_basis_f->at(i));
//           for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
//           {
//             const auto& basis_j = *(m_basis_f->at(j));
//             double val = 0.0;
//             for (int qp = 0; qp < quadrature->size(); ++qp)
//             {
//               val += q_weights[qp] * ((*m_velocity_fnc)(q_points[qp]) * normal)
//                      * basis_i(q_points[qp], elemVertices)
//                      * basis_j(q_points[qp], elemVertices);
//             }
//             mat_K(i,j) += val;
//           }
//         }
//       }
//       // end calculation over gamma plus edge
//
//       if (bc_info_edges[edge_ind].m_boundary_condition_type == typeSigma)
//       {
//         zeroOutArray(mat_sigma);
//         const auto& point_inds = mesh_f->getEdgeGlobalIndexes(edge_ind);
//         const auto edge_verts = mesh_f->getEdgeVertices(edge_ind);
//
//         // compute the element limiter alpha_w here from the nodal limiters!
//         {
//           alpha_w = std::min(nodal_lim_w_view[bc_info_pts[point_inds[0]].m_index_in_system],
//                              nodal_lim_w_view[bc_info_pts[point_inds[1]].m_index_in_system]);
//         }
//         // get quadrature on this edge
//         const auto& quadrature = m_quadrature_bdry->at(edge_ind);
//         // integrate
//         const auto& q_points = quadrature->get_q_points();
//         const auto& q_weights = quadrature->get_q_weights();
//         for (int qp = 0; qp < quadrature->size(); ++qp)
//         {
//           // first compute conc_f (this is a hack on fixed mesh with a horizontal edge)
//           double uh_f = 0.0;
//           double uh_w = 0.0;
//           double uh_dot_w = 0.0;
//           for (std::size_t i = 0; i < point_inds.size(); ++i)
//           {
//             const int& ind_i_f = point_inds[i];
//             const int  ind_i_w = bc_info_pts[point_inds[i]].m_index_in_system + mesh_f->numOfPoints();
//             const auto& basis_i = *(m_basis_w->at(i));
//             uh_f += basis_i(q_points[qp],edge_verts) * U_view[ind_i_f];
//             uh_w += basis_i(q_points[qp],edge_verts) * U_view[ind_i_w];
//             uh_dot_w += basis_i(q_points[qp],edge_verts) * U_dot_view[ind_i_w];
//           }
//
//           for (std::size_t i = 0; i < point_inds.size(); ++i)
//           {
//             const int& ind_i_f = point_inds[i];
//             const int  ind_i_w = bc_info_pts[point_inds[i]].m_index_in_system + mesh_f->numOfPoints();
//             const auto& basis_i = *(m_basis_w->at(i));
//             res_view[ind_i_f] += k_d * q_weights[qp] * basis_i(q_points[qp],edge_verts) *
//               (alpha_f * m_lambda(uh_f) + (1.0 - alpha_f) * m_lambda_tilde(uh_f)*U_view[ind_i_f]
//               - alpha_w * uh_w - (1.0 - alpha_w) * U_view[ind_i_w]);
//             res_view[ind_i_w] -= k_d * q_weights[qp] * basis_i(q_points[qp],edge_verts) *
//                                  (alpha_f * m_lambda(uh_f) + (1.0 - alpha_f) * m_lambda_tilde(uh_f)*U_view[ind_i_f]
//                                   - alpha_w * uh_w - (1.0 - alpha_w) * U_view[ind_i_w])
//                                   // now add time derivative terms
//                                   - alpha_w * q_weights[qp] * basis_i(q_points[qp],edge_verts) * uh_dot_w
//                                   - (1.0 - alpha_w) * q_weights[qp] * basis_i(q_points[qp],edge_verts) * U_dot_view[ind_i_w];
//
//
//             for (std::size_t j = 0; j < point_inds.size(); ++j)
//             {
//               const int& ind_j_f = point_inds[j];
//               const int& ind_j_w = bc_info_pts[point_inds[j]].m_index_in_system + mesh_f->numOfPoints();
//               const auto& basis_j = *(m_basis_w->at(j));
//
//               mat_sigma(i,j)
//                 += alpha_f * k_d * q_weights[qp] * basis_i(q_points[qp],edge_verts) * m_lambda_tilde(uh_f) * basis_j(q_points[qp],edge_verts);
//
//               mat_sigma(i+point_inds.size(),j+point_inds.size())
//                 += alpha_w * k_d * q_weights[qp] * basis_i(q_points[qp],edge_verts) * basis_j(q_points[qp],edge_verts)
//                   // add the time derivative part
//                    - beta * alpha_w * q_weights[qp] * basis_i(q_points[qp],edge_verts) * basis_j(q_points[qp],edge_verts);
//
//               mat_sigma(i,j+point_inds.size()) -= alpha_w * k_d * q_weights[qp] * basis_i(q_points[qp],edge_verts) * basis_j(q_points[qp],edge_verts);
//
//               mat_sigma(i+point_inds.size(),j) -= alpha_f * k_d * q_weights[qp] * basis_i(q_points[qp],edge_verts) * m_lambda_tilde(uh_f) * basis_j(q_points[qp],edge_verts);
//
//               if (i == j)
//               {
//                 mat_sigma(i,j) +=
//                   (1.0 - alpha_f) * k_d * q_weights[qp] * basis_i(q_points[qp], edge_verts) * m_lambda_tilde(uh_f);
//
//                 mat_sigma(i+point_inds.size(),j+point_inds.size()) +=
//                   (1.0 - alpha_w) * k_d * q_weights[qp] * basis_i(q_points[qp], edge_verts)
//                   // add the time derivative
//                   - beta * (1.0 - alpha_w) * q_weights[qp] * basis_i(q_points[qp], edge_verts);
//
//                 mat_sigma(i,j+point_inds.size()) -= (1.0 - alpha_w) * k_d * q_weights[qp] * basis_i(q_points[qp], edge_verts);
//
//                 mat_sigma(i+point_inds.size(),j) -= (1.0 - alpha_f) * k_d * q_weights[qp] * basis_i(q_points[qp], edge_verts) * m_lambda_tilde(uh_f);
//               }
//             }
//           }
//
//         }
//         // end loop over quadrature
//
//         // add the entries into the Jacobian
//         for (std::size_t i = 0; i < point_inds.size(); ++i)
//         {
//           const int& ind_i_f = point_inds[i];
//           const int  ind_i_w = bc_info_pts[point_inds[i]].m_index_in_system + mesh_f->numOfPoints();
//           for (std::size_t j = 0; j < point_inds.size(); ++j)
//           {
//             const int& ind_j_f = point_inds[j];
//             const int& ind_j_w = bc_info_pts[point_inds[j]].m_index_in_system + mesh_f->numOfPoints();
//             jac->sumIntoGlobalValues(ind_i_f,1,&(mat_sigma(i,j)),&ind_j_f);
//             jac->sumIntoGlobalValues(ind_i_f,1,&(mat_sigma(i,j+point_inds.size())),&ind_j_w);
//             jac->sumIntoGlobalValues(ind_i_w,1,&(mat_sigma(i+point_inds.size(),j+point_inds.size())),&ind_j_w);
//             jac->sumIntoGlobalValues(ind_i_w,1,&(mat_sigma(i+point_inds.size(),j)),&ind_j_f);
//           }
//         }
//
//       }
//       // end calculation over sigma edge
//     }
//     // end loop over local edges
//     // use the conv matrix to compute artificial diffusion
//     AFC_LocalDiffusionMatrix(mat_D, mat_K);
//
//     KokkosBlas::gemv("N",1.0,mat_K,U_loc,0.0,res_loc_f);
//     KokkosBlas::gemv("N",1.0,mat_S,U_loc,1.0,res_loc_f);
//     KokkosBlas::gemv("N",1.0-alpha_f,mat_D,U_loc,1.0,res_loc_f);
//     KokkosBlas::gemv("N",-alpha_f,mat_M,U_dot_loc,1.0,res_loc_f);
//     KokkosBlas::gemv("N",-(1.0 - alpha_f),mat_Ml,U_dot_loc,1.0,res_loc_f);
//
//     // now assemble the residual & first diagonal block of the Jacobian
//     const std::vector<int>& cols = glob_inds;
//     std::vector<double> vals(m_dofmapper->local_ndof(),0.0);
//     for (int i = 0; i < m_dofmapper->local_ndof(); ++i)
//     {
//       res_view[glob_inds[i]] += res_loc_f[i];
//       for (int j = 0; j < m_dofmapper->local_ndof(); ++j)
//       {
//         vals[j] = - beta * (alpha_f * mat_M(i,j) + (1.0 - alpha_f) * mat_Ml(i,j))
//                   + (mat_K(i,j) + (1.0 - alpha_f) * mat_D(i,j) + mat_S(i,j));
//       }
//       jac->sumIntoGlobalValues(glob_inds[i],static_cast<int>(cols.size()),vals.data(),cols.data());
//     }
//   }
//   // end loop over elements
//   jac->fillComplete();
  
//   if (!m_jac_applied)
//   {
//     buildStiffMatrix(jac_U,m_basis,m_quadrature,m_dofmapper);
//     m_jac_applied = true;
//   }
//   auto prob = std::dynamic_pointer_cast<Problem_Bioseparation>(m_problem);
//   buildRHSVector(m_rhs,prob->rhsFunction(),m_basis,m_quadrature,m_quadrature);
//   *res_U = (*jac_U)*(*U) - (*m_rhs);



}

void Assembler_Bioseparation::
applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                 const std::shared_ptr<FEMatrix>& jac_U) const
{
  auto bc_ = std::dynamic_pointer_cast<BC_Bioseparation>(m_problem->bc());
  if (bc_)
  {
    for (auto it : bc_->bCPoints())
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
