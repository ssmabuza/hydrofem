// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Problem_Poisson_HPP__
#define __Hydrofem_Problem_Poisson_HPP__

#include "Hydrofem_BC.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_InitialCondition.hpp"

namespace hydrofem
{

class Problem_Poisson
  :
  public Problem
{
public:

  // creates a standard LPS stabilized problem
  explicit Problem_Poisson(const std::string& name = "poisson") : public Problem()
  {
    setName(name);
    setDofNames({{"u"}});
  }

  ~Problem_Poisson() override = default;
  
  void init() override;
  
  void set_exact(const std::shared_ptr<ScalarAnalyticalExpression>& exact_)
  { m_exact = exact_; }
  
  void set_bc(const std::shared_ptr<BC>& bc_)
  { m_bc = bc_; }
  
  std::shared_ptr<ScalarAnalyticalExpression> rhsFunction() const
  { return m_rhs_fnc; }
  
  std::shared_ptr<ScalarAnalyticalExpression> exact() const
  { return m_exact; }
  
  [[nodiscard]] std::shared_ptr<BC> bc() const override
  { return m_bc; }
  
  std::shared_ptr<DirichletBCFunction> dirichletBCFunction() const
  { return m_dirichlet_bc_fnc; }
  
protected:

  // BC base
  std::shared_ptr<BC> m_bc;
  // exact soln
  std::shared_ptr<ScalarAnalyticalExpression> m_exact;
  // source function
  std::shared_ptr<ScalarAnalyticalExpression> m_rhs_fnc;
  // Dirichlet boundary values
  std::shared_ptr<DirichletBCFunction> m_dirichlet_bc_fnc;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_Poisson_HPP__ */
