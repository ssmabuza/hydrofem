// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Problem_Bioseparation_HPP__
#define __Hydrofem_Problem_Bioseparation_HPP__

#include "Hydrofem_Problem.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

class Problem_Bioseparation
  :
  public Problem, public Optionable
{
public:

  using Ptr = std::shared_ptr<Problem_Bioseparation>;

  // creates a standard LPS stabilized problem
  explicit Problem_Bioseparation(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Problem(),
    Optionable(option_handler)
  {
    setName("bioseparation");
    setDofNames({{"conc"}});
    option_handler->parse();
    m_is_initialized = false;
  }

  ~Problem_Bioseparation() override = default;
  
  void init() override;

  virtual std::shared_ptr<BC> getBoundaryCondition() const override
  { return m_bc; }
  
  virtual std::shared_ptr<ScalarInitialCondition> getInitialConditionFunction() const override
  { return m_ic; }

  virtual std::shared_ptr<ScalarAnalyticalExpression> getBoundaryFunction() const override
  { return m_u_in; }

  // exact solution is not known in this case

  // get model constants
  double omega() const { return m_omega; }
  double rho_s() const { return m_rho_s; }
  double q_max() const { return m_q_max; }
  double Keq() const { return m_K_eq; }
  double flowrate() const { return m_flowrate; }
  double width() const { return m_xf - m_x0; }
  double alphaT() const { return m_alphaT; }
  double alphaL() const { return m_alphaL; }
  double d0() const { return m_d0; }

private:

  virtual void addOptionsCallback(po::options_description &config);

  // model constants
  //@{
  double m_x0, m_xf;
  double m_omega;
  double m_rho_s;
  double m_q_max;
  double m_K_eq;
  double m_flowrate;
  double m_alphaL;
  double m_alphaT;
  double m_d0;
  //@}

  // boundary condition
  std::shared_ptr<BC> m_bc;
  // initial condition function
  std::shared_ptr<ScalarInitialCondition> m_ic;
  // boundary function
  std::shared_ptr<ScalarAnalyticalExpression> m_u_in;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_Bioseparation_HPP__ */
