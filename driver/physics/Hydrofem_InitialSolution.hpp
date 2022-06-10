// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_InitialSolution_HPP__
#define __Hydrofem_InitialSolution_HPP__

// #include "Term.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Equation.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_GlobalGather.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_GlobalScatter.hpp"
#include "Hydrofem_InitialCondition.hpp"
#include "Hydrofem_GlobalSystemGather.hpp"
#include "Hydrofem_GlobalSystemScatter.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem
{

/**
 * \brief Nodal projection of the initial data
 */
template <typename IC_Type>
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
  InitialSolution(const std::shared_ptr<LinearObjectBuilder>& lob) : m_lob(lob)
  {
    m_is_computed = false;
  }
  
  //! Dtor
  virtual ~InitialSolution() {}
  
  //! \brief evaluate values 
  virtual void evaluate() = 0;
  
  //! \brief get evaluated field
  virtual std::shared_ptr<std::vector<MDArray>>
  get_evaluatedField() const;

  virtual std::shared_ptr<FEVector>
  get_evaluatedFieldGlobal() const;

  /** \brief Fixed consistent mass matrix for explicit steppers */
  std::shared_ptr<FEMatrix>
  getConsistentMassMatrix() const;

  /** \brief Fixed lumped mass matrix for explicit steppers */
  std::shared_ptr<FEVector>
  getLumpedMassMatrix() const;
  
  /** \brief gets the app linear obj builder which is m_lob by default */
  virtual inline
  std::shared_ptr<LinearObjectBuilder>
  getAppLOB() const
  { return m_lob; }

  void buildConsistentMassMatrix()
  { throw std::runtime_error("Consistent mass matrix builder not implemented."); }

protected:
  
  // the linear object builder
  std::shared_ptr<LinearObjectBuilder>  m_lob;
  // true if initial solution is computed successfully
  mutable bool                          m_is_computed;
  // the initial solution as a gathered vector (comp x (cells * local_ndofs))
  std::shared_ptr<std::vector<MDArray>> m_result_gathered;
  // the initial solution as a distributed vector
  std::shared_ptr<FEVector>             m_result;
  // flag for consistent mass matrix
  bool                                  m_mass_matrix_built = false;
  // consistent mass matrix
  std::shared_ptr<FEMatrix>             m_mass;
  // flag for lumped mass matrix
  bool                                  m_lumped_mass_matrix = false;
  // lumped mass matrix
  mutable std::shared_ptr<FEVector>     m_lumped_mass;
  
};

}
// end namespace valiant

#endif /** __Valiant_InitialSolution_HPP__ */
