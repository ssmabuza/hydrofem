// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Assembler_CDR_HPP__
#define __Hydrofem_Assembler_CDR_HPP__


#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_AFC_Limiter.hpp"
#include "Hydrofem_Assembler_Base.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem
{

class Assembler_CDR
  :
  public Assembler_Base
{
public:
  
  using LOB = Assembler_Base::LOB;

  /** \brief Ctor */
  Assembler_CDR(const std::shared_ptr<Problem>& problem,
                const std::shared_ptr<DofMapper>& dofmapper,
                const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature, bool /*afc_enabled*/)
  {
    m_problem = problem;
    m_dofmapper = dofmapper;
    m_basis = basis;
    m_quadrature = quadrature;
    m_mesh = dofmapper->mesh();
    m_problem = problem;
    m_lob = std::make_shared<LOB>(m_dofmapper);
    m_afc_vec = m_lob->createVector();
    //m_do_afc = afc_enabled;
  }
  
  /** \brief Dtor */
  ~Assembler_CDR() override = default;
  
  void buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                                const std::shared_ptr<const FEVector>& U_dot,
                                const std::shared_ptr<FEVector>& res_U,
                                const std::shared_ptr<FEMatrix>& jac_U,
                                const double time,
                                const double delta_t,
                                const double beta) const override;

  /** \brief  */
  void applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                        const std::shared_ptr<FEMatrix>& jac_U) const override;
                                
protected:

  // problem defining the Bioseparation equation and parameters
  std::shared_ptr<Problem> m_problem;
  // limiter 
  std::shared_ptr<AFC_Limiter> m_limiter;
  // storage for limiter
  std::shared_ptr<FEVector> m_nodal_limiter;
  // mesh from dof manager
  std::shared_ptr<Mesh> m_mesh;
  // the dof manager
  std::shared_ptr<DofMapper> m_dofmapper;
  // basis functions
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // global quadrature
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
  // AFC vector
  std::shared_ptr<FEVector> m_afc_vec;
  // do algebraic flux correction
  //bool m_do_afc;
  // AFC edge graphs
  // global mass matrix
  std::shared_ptr<FEMatrix> m_M_mat;
  // global artificial diffusion matrix
  std::shared_ptr<FEMatrix> m_D_mat;
  // flag for once off builds
  mutable bool m_built_d_mat_graph = false;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Assembler_CDR_HPP__ */
