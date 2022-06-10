// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Valiant_InitialSolution_HyperbolicSystem_HPP__
#define __Valiant_InitialSolution_HyperbolicSystem_HPP__

#include "InitialSolution.hpp"

namespace valiant
{

class InitialSolution_HyperbolicSystem
  :
  public InitialSolution
{
public:
  
  using MDArray = FEArray<double>::CellBasis;
  
  /**
   * @brief Ctor
   *
   * @param n_eq - number of equations in the hyperbolic system
   * @param dofmapper
   * @param fe_basis
   * @param quadrature
   * @param ic
   * @param sys_lob
   * @param proj_type
   *
   */
  InitialSolution_HyperbolicSystem(const int n_eq,
                                   const Teuchos::RCP<DofMapper>& dofmapper,
                                   const Teuchos::RCP<std::vector<Teuchos::Ptr<FEBasis>>>& fe_basis,
                                   const Teuchos::RCP<std::vector<Teuchos::Ptr<Quadrature>>>& quadrature,
                                   const Teuchos::RCP<InitialCondition>& ic,
                                   const Teuchos::RCP<SystemLinearObjectBuilder<int>>& sys_lob,
                                   const std::string proj_type = "Lumped")
    :
    InitialSolution(sys_lob->getUnderLyingScalarLOB()),
    m_neq(n_eq),
    m_dofmapper(dofmapper),
    m_ic(ic),
    m_fe_basis(fe_basis),
    m_quadrature(quadrature),
    m_sys_lob(sys_lob),
    m_proj_type(proj_type)
  {
    m_is_computed = false;
  
    m_result_gathered
      = Teuchos::rcp(new std::vector<MDArray>(m_neq,
                                              createKArray<MDArray>(m_dofmapper->nelements(),
                                                                    m_dofmapper->local_ndof())));
  
    m_global_gather = Teuchos::rcp(new GlobalSystemGather(m_neq,dofmapper));
    m_global_scatter = Teuchos::rcp(new GlobalSystemScatter(m_neq,dofmapper));
    m_scalar_gather = Teuchos::rcp(new GlobalGather(dofmapper));
    m_scalar_scatter = Teuchos::rcp(new GlobalScatter(dofmapper));
    m_result = Teuchos::rcp(new TVector(m_sys_lob->getOwnedMap()));
    m_lob = m_sys_lob->getUnderLyingScalarLOB();
  }
  
  /**
   * @brief Dtor
   */
  virtual ~InitialSolution_HyperbolicSystem() {}
  
  //! \brief evaluate values
  virtual void evaluate();
  
  /** \brief gets the global gather for the system */
  Teuchos::RCP<GlobalSystemGather>
  getGlobalSystemGather() const
  { return m_global_gather; }
  
  /** \brief gets the global scatter for the system */
  Teuchos::RCP<GlobalSystemScatter>
  getGlobalSystemScatter() const
  { return m_global_scatter; }
  
  /** \brief gets the global gather for a scalar */
  Teuchos::RCP<GlobalGather>
  getGlobalScalarGather() const
  { return m_scalar_gather; }
  
  /** \brief gets the global scatter for a scalar */
  Teuchos::RCP<GlobalScatter>
  getGlobalScalarScatter() const
  { return m_scalar_scatter; }
  
  /** \brief gets the linear obj builder for the system */
  virtual Teuchos::RCP<LinearObjectBuilder<int>>
  getAppLOB() const
  { return m_sys_lob; }

private:
  
  // number of equations in the hyperbolic system
  int                                                 m_neq;
  // the dof manager for each var
  Teuchos::RCP<DofMapper>                             m_dofmapper;
  // the initial condition function
  Teuchos::RCP<InitialCondition>                      m_ic;
  // the FE basis
  Teuchos::RCP<std::vector<Teuchos::Ptr<FEBasis>>>    m_fe_basis;
  // the FE quadrature
  Teuchos::RCP<std::vector<Teuchos::Ptr<Quadrature>>> m_quadrature;
  // this is the system gather into an owned container
  Teuchos::RCP<GlobalSystemGather>                    m_global_gather;
  // this is the system scatter into an owned container
  Teuchos::RCP<GlobalSystemScatter>                   m_global_scatter;
  // this is the system gather into an owned container
  Teuchos::RCP<GlobalGather>                          m_scalar_gather;
  // this is the system scatter into an owned container
  Teuchos::RCP<GlobalScatter>                         m_scalar_scatter;
  // the system linear object builder
  Teuchos::RCP<SystemLinearObjectBuilder<int>>        m_sys_lob;
  // the scalar linear object builder
  Teuchos::RCP<LinearObjectBuilder<int>>              m_lob;
  // projection type for the initial condition: Lumped, Consistent, AFC
  std::string                                         m_proj_type;
  
  //! @brief the routine to build the lumped mass matrix
  void buildLumpedMassMatrix();

  //! @brief the routine to build the consistent mass matrix
  void buildConsistentMassMatrix();
  
  //! @brief
  void evaluateLumpedMassProjection();
  
  //! @brief
  void evaluateConsistentProjection();
  
};

}
// end namespace valiant

#endif /** __Valiant_InitialSolution_HyperbolicSystem_HPP__ */
