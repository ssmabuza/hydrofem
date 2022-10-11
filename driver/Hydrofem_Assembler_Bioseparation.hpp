// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Assembler_Bioseparation_HPP__
#define __Hydrofem_Assembler_Bioseparation_HPP__


#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_AFC_Limiter.hpp"
#include "Hydrofem_Assembler_Base.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem
{

class Assembler_Bioseparation
  :
  public Assembler_Base
{
public:
  
  using LOB = Assembler_Base::LOB;

  /** \brief Ctor */
  Assembler_Bioseparation(const std::shared_ptr<Problem>& problem,
                          const std::shared_ptr<DofMapper>& dofmapper,
                          const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                          const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature, bool afc_enabled)
  {
    m_problem = problem;
    m_dofmapper = dofmapper;
    m_basis = basis;
    m_quadrature = quadrature;
    m_mesh = dofmapper->mesh();
    m_problem = problem;
    m_lob = std::make_shared<LOB>(m_dofmapper);
    m_afc_vec = m_lob->createVector();
    m_do_afc = afc_enabled;
  }
  
  /** \brief Dtor */
  ~Assembler_Bioseparation() override = default;
  
  void buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                                const std::shared_ptr<const FEVector>& U_dot,
                                const std::shared_ptr<FEVector>& res_U,
                                const std::shared_ptr<FEMatrix>& jac_U,
                                const double time,
                                const double delta_t,
                                const double beta) const override;

  void buildSteadyResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                                      const std::shared_ptr<const FEVector>& U_dot,
                                      const std::shared_ptr<FEVector>& res_U,
                                      const std::shared_ptr<FEMatrix>& jac_U,
                                      const double time,
                                      const double delta_t,
                                      const double beta) const override;

  /** \brief  */
  void applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                        const std::shared_ptr<FEMatrix>& jac_U) const override;
                                
private:

  // for alternative matrix based implementation
  // struct MatrixBasedImpl
  // {

  //   void init();

  //   void computeGlobalConvectionPlusOutletMatrix();

  //   void computeGlobalAdsorptionMatrix();

  //   void computeGlobalStiffnessMatrix();

  //   void computeGlobalLumpedMatrix();


  //   bool m_initialized = false;

  //   std::shared_ptr<FEMatrix> m_mat_M;
  //   std::shared_ptr<FEMatrix> m_mat_S;
  //   std::shared_ptr<FEMatrix> m_mat_K;

  // };

  // AFC calculations to be done
  // struct AFC_Group
  // {

  //   void computeFluxEdgeGraph();

  //   std::shared_ptr<FEMatrix> m_mat_D;
  //   std::shared_ptr<FEVector> m_mat_Ml;

  //   std::shared_ptr<FEVector> m_graph_M;
  //   std::shared_ptr<FEVector> m_graph_D;
  //   std::shared_ptr<FEVector> m_graph_F;

  // };

  // residual based implementation
  // void buildResidualAndJacobian_ResidualBased() const;

  // matrix based implementation
  // void buidResidualAndJacobian_MatrixAndGraphBased() const;

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
  bool m_do_afc;
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

#endif /** __Hydrofem_Assembler_Bioseparation_HPP__ */
