// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Problem_HyperbolicSystem_HPP__
#define __Hydrofem_Problem_HyperbolicSystem_HPP__

#include "Hydrofem_BC.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_Equation.hpp"
#include "Hydrofem_InitialCondition.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{

/**
 * \brief A class that describes a continuous problem : hyperbolic system
 */
class Problem_HyperbolicSystem
  :
  public Problem
{
public:

  Problem_HyperbolicSystem() = default;

  ~Problem_HyperbolicSystem() override = default;

  //! \brief get functions
  //@{
  //! get the inlet function
  virtual std::shared_ptr<AnalyticalExpression> Uin() const
  { return m_Uin; }

  //! \brief get the exact solution
  virtual std::shared_ptr<AnalyticalExpression> exact() const
  { return m_exact; }

  //! \brief get the equation
  virtual std::shared_ptr<Equation<RealType>> equation() const
  { return m_equation; }

  //! \brief get the initial solution
  virtual std::shared_ptr<InitialCondition> ic() const
  { return m_ic; }

  virtual std::shared_ptr<BC> bc()
  { return m_bc; }
  //@}

  //! \brief set functions
  //@{
  void set_bc(const std::shared_ptr<BC>& bc)
  { m_bc = bc; }
  
  void set_Uin(const std::shared_ptr<AnalyticalExpression>& uin)
  { m_Uin = uin; }

  void set_exact(const std::shared_ptr<AnalyticalExpression>& exact)
  { m_exact = exact; }

  void set_ic(const std::shared_ptr<InitialCondition>& ic)
  { m_ic = ic; }
  //@}

protected:

  // BC base
  mutable std::shared_ptr<BC>           m_bc;
  // inlet values
  std::shared_ptr<AnalyticalExpression> m_Uin;
  // exact soln
  std::shared_ptr<AnalyticalExpression> m_exact;
  // source function
  std::shared_ptr<AnalyticalExpression> m_source;
  // equation
  std::shared_ptr<Equation<RealType>>   m_equation;
  // IC
  std::shared_ptr<InitialCondition>     m_ic;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_HyperbolicSystem_HPP__ */
