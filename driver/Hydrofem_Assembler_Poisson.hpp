// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Assembler_Poisson_HPP__
#define __Hydrofem_Assembler_Poisson_HPP__


#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_Assembler_Base.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{

class Assembler_Poisson
  :
  public Assembler_Base
{
public:
  
  using LOB = Assembler_Base::LOB;

  /** \brief Ctor */
  Assembler_Poisson(const std::shared_ptr<Problem>& problem,
                    const std::shared_ptr<DofMapper>& dofmapper,
                    const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                    const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature)
  {
    m_dofmapper = dofmapper;
    m_basis = basis;
    m_quadrature = quadrature;
    m_mesh = dofmapper->mesh();
    m_problem = problem;
    m_lob = std::make_shared<LOB>(m_dofmapper);
    m_rhs = m_lob->createVector();
    m_jac_applied = false;
    m_dirichlet_applied = false;
  }
  
  /** \brief Dtor */
  ~Assembler_Poisson() override = default;
  
  /**
   * \brief The non-transient assembly (steady solver)
   *
   * \param U     - current solution
   * \param res_U - residual
   * \param jac_U - Jacobian
   * \param beta  - system multiplier
   */
  void buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                                const std::shared_ptr<FEVector>& res_U,
                                const std::shared_ptr<FEMatrix>& jac_U,
                                const double beta) const override;
  
private:
  
  /** \brief  */
  [[maybe_unused]] void
  buildStiffMatrix(const std::shared_ptr<FEMatrix>& stiff,
                   const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                   const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                   const std::shared_ptr<const DofMapper>& dofmapper) const;
  
  /** \brief  */
  void buildRHSVector(const std::shared_ptr<FEVector>& rhs,
                      const std::shared_ptr<ScalarAnalyticalExpression>& f,
                      const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                      const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                      const std::shared_ptr<const DofMapper>& dofmapper) const;
  
  /** \brief  */
  void applyDirichletBC(const std::shared_ptr<FEVector>& res_U,
                        const std::shared_ptr<FEMatrix>& jac_U) const override;

  // problem defining the Poisson equation
  std::shared_ptr<Problem> m_problem;
  // mesh from dof manager
  std::shared_ptr<Mesh> m_mesh;
  // the dof manager
  std::shared_ptr<DofMapper> m_dofmapper;
  // basis functions
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // global quadrature
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
  // fixed RHS vector
  mutable std::shared_ptr<FEVector> m_rhs;
  // Fixed Jacobian applied
  mutable bool m_jac_applied;
  // Dirichlet BC applied to Jacobian
  mutable bool m_dirichlet_applied;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Assembler_Poisson_HPP__ */
