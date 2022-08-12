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

  // creates a standard LPS stabilized problem
  explicit Problem_Bioseparation(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Problem(),
    Optionable(option_handler)
  {
    setName("bioseparation");
    setDofNames({{"conc"}});
    parse();
  }

  ~Problem_Bioseparation() override = default;
  
  void init() override;

  double omega() const { return m_omega; }
  double rho_s() const { return m_rho_s; }
  double q_max() const { return m_q_max; }
  double Keq() const { return m_K_eq; }
  double flowrate() const { return m_flowrate; }
  double width() const { return m_width; }

  virtual std::shared_ptr<BC> getBC() const
  { return m_bc; }

  virtual std::shared_ptr<ScalarInitialCondition> ic_scalar() const
  { return m_ic; }

  virtual std::shared_ptr<ScalarAnalyticalExpression> uIn() const
  { return m_u_in; }

private:

  virtual void addOptionsCallback(po::options_description &config);

  // model constants
  //@{
  double m_omega;
  double m_rho_s;
  double m_q_max;
  double m_K_eq;
  double m_flowrate;
  double m_width;
  //@}

  std::shared_ptr<BC> m_bc;
  std::shared_ptr<ScalarInitialCondition> m_ic;
  std::shared_ptr<ScalarAnalyticalExpression> m_u_in;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_Bioseparation_HPP__ */
