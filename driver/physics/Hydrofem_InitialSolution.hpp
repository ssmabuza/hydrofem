// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_InitialSolution_HPP__
#define __Hydrofem_InitialSolution_HPP__

#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_GlobalGather.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_InitialCondition.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem
{

//! TODO:   Specialize for scalars and vectors

/**
 * \brief Projection of the initial condition function into system vector (default: nodal)
 */
class InitialSolution 
{
public:

  using MDArray = FEArray<double>::CellBasis;
  
  /**
   * \brief default Ctor
   */
  InitialSolution()
  {
    m_is_computed = false;
  }
  
  /**
   * \brief Ctor
   */
  InitialSolution(const std::shared_ptr<LinearObjectBuilder>& lob,
                  const std::shared_ptr<ScalarInitialCondition>& ic_fnc) : m_lob(lob)
  {
    m_ic = ic_fnc;
    m_is_computed = false;
  }

  //! Dtor
  virtual ~InitialSolution() = default;
  
  //! \brief evaluate values 
  virtual void evaluate() = 0;
  
  virtual std::shared_ptr<FEVector>
  get_evaluatedField() const;

  /** \brief gets the app linear obj builder which is m_lob by default */
  virtual inline
  std::shared_ptr<LinearObjectBuilder>
  getAppLOB() const
  { return m_lob; }

protected:

  // initial cond function (scalar)
  std::shared_ptr<ScalarInitialCondition> m_ic;
  // the linear object builder
  std::shared_ptr<LinearObjectBuilder> m_lob;
  // true if initial solution is computed successfully
  mutable bool m_is_computed;
  // the initial solution as a distributed vector
  std::shared_ptr<FEVector> m_result;
  // flag for lumped mass matrix
  bool m_lumped_mass_matrix_built = false;
  // flag for consistent mass matrix
  bool m_mass_matrix_built = false;
  
};

/**
 * \brief nodal projection on nodal Dofs on a mesh i.e. mesh node is dof node
 */
class NodalProjection
  :
  public InitialSolution
{
public:

  NodalProjection(const std::shared_ptr<LinearObjectBuilder>& lob,
                  const std::shared_ptr<ScalarInitialCondition>& ic_fnc) : InitialSolution(lob,ic_fnc)
  {
    m_result = std::make_shared<FEVector>(m_lob->dofMapper()->mesh()->numOfPoints());
  }

  void evaluate() override;
  
};

class ConsistentMassProjection
  :
  public InitialSolution
{
public:

  ConsistentMassProjection(const std::shared_ptr<LinearObjectBuilder>& lob,
                           const std::shared_ptr<ScalarInitialCondition>& ic_fnc,
                           const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                           const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature);

  void evaluate() override;

private:

  // consistent mass matrix
  std::shared_ptr<FEMatrix> m_mass;
  // FE Basis
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // Quadrature Rule
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;

};

class LumpedMassProjection
  :
  public InitialSolution
{
public:

  LumpedMassProjection(const std::shared_ptr<LinearObjectBuilder>& lob,
                       const std::shared_ptr<ScalarInitialCondition>& ic_fnc,
                       const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                       const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature);

  void evaluate() override;

private:

  // lumped mass matrix
  std::shared_ptr<FEVector> m_lumped_mass;
  // FE Basis
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // Quadrature Rule
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_InitialSolution_HPP__ */
