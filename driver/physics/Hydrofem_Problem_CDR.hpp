// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Problem_CDR_HPP__
#define __Hydrofem_Problem_CDR_HPP__

#include "Hydrofem_Problem.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

class Problem_CDR
  :
  public Problem, public Optionable
{
public:

  using Ptr = std::shared_ptr<Problem_CDR>;

  // creates a standard LPS stabilized problem
  explicit Problem_CDR(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Problem(),
    Optionable(option_handler)
  {
    setName("cdr");
    setDofNames({{"conc"}});
    option_handler->parse();
    m_is_initialized = false;
  }

  ~Problem_CDR() override = default;
  
  void init() override;

  virtual std::shared_ptr<BC> getBoundaryCondition() const override
  { return m_bc; }
  
  virtual std::shared_ptr<ScalarInitialCondition> getInitialConditionFunction() const override
  { return m_ic; }

  virtual std::shared_ptr<ScalarAnalyticalExpression> getBoundaryFunction() const override
  { return m_u_in; }

  // exact solution is not known in this case

  // get model constants
  double L() const { return m_L; }
  double v_max() const { return m_v_max; }
  double diff() const { return m_diff; }
  double H() const { return m_H; }

private:

  virtual void addOptionsCallback(po::options_description &config);

  // model constants
  //@{
  double m_L;
  double m_v_max;
  double m_diff;
  double m_H;
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
