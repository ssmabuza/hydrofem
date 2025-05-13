// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Problem_Poisson_HPP__
#define __Hydrofem_Problem_Poisson_HPP__

#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_InitialCondition.hpp"

namespace hydrofem
{

class BC_Scalar;

class Problem_Poisson
  :
  public Problem
{
public:

  // creates a standard LPS stabilized problem
  explicit Problem_Poisson(const std::string& name = "poisson") : Problem()
  {
    setName(name);
    setDofNames({{"u"}});
  }

  ~Problem_Poisson() override = default;
  
  void init() override;
  
  void setExactSolution(const std::shared_ptr<ScalarAnalyticalExpression>& exact_)
  { m_exact = exact_; }
  
  void setBoundaryCondition(const std::shared_ptr<BC>& bc_)
  { m_bc = bc_; }
  
  [[nodiscard]] std::shared_ptr<BC> getBoundaryCondition() const override
  { return m_bc; }
  
  std::shared_ptr<ScalarAnalyticalExpression> getExactSolution() const override
  { return m_exact; }
  
  std::shared_ptr<ScalarAnalyticalExpression> getBoundaryFunction() const
  { return m_dirichlet_bc_fnc; }

  // RHS function in Poisson equation strong form
  std::shared_ptr<ScalarAnalyticalExpression> getRHSFunction() const
  { return m_rhs_fnc; }

  // steady problems do not need initial condition function
  
protected:

  // BC base
  std::shared_ptr<BC> m_bc;
  // exact soln
  std::shared_ptr<ScalarAnalyticalExpression> m_exact;
  // source function
  std::shared_ptr<ScalarAnalyticalExpression> m_rhs_fnc;
  // Dirichlet boundary values
  std::shared_ptr<ScalarAnalyticalExpression> m_dirichlet_bc_fnc;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_Poisson_HPP__ */
