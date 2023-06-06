// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Stepper_Classic_Theta_HPP__
#define __Hydrofem_Stepper_Classic_Theta_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Stepper.hpp"
#include "Hydrofem_NewtonSolver.hpp"
#include "Hydrofem_Assembler_Base.hpp"
#include "Hydrofem_InitialSolution.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem 
{

/**
 * \brief The theta scheme that builds a steady part of the residual and the
 *        manually forms the full residual and Jacobian
 */
class Stepper_Classic_Theta
  :
  public Stepper
{
public:
  
  //! Ctor
  Stepper_Classic_Theta(const std::shared_ptr<InitialSolution>& ic,
                        const std::shared_ptr<Assembler_Base>& assembler,
                        const std::shared_ptr<OptionHandler>& option_handler)
    :
    Stepper(option_handler,assembler)
  {
    // parse the options
    option_handler->parse();
    // get initial condition
    m_ic = ic;
    // get linear object builder 
    const auto lob = m_ic->getAppLOB();
    // build solution and residual vectors
    m_u_new = m_ic->get_evaluatedField();
    // previous solution
    m_u_old = lob->createVector(); m_u_old->setZero();
    // time derivative of solution
    m_u_dot = lob->createVector(); m_u_dot->setZero();
    // previous time derivative of solution
    m_u_dot_old = lob->createVector(); m_u_dot_old->setZero();
    // build Jacobian matrix
    m_jac = lob->createSparseMatrix();
    // initialize the matrix
    lob->buildSparseGraph(m_jac);
    // build the nonlinear solver
    m_nlsolver = std::make_shared<NewtonSolver>(option_handler,m_assembler);
    // initialize nonlinear solver
    m_nlsolver->initialize();
  }

  //! Dtor
  ~Stepper_Classic_Theta() override = default;
  
  /** \brief options to be parsed for stepper */
  virtual void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("stepper-theta",po::value<double>(&m_theta)->default_value(0.5),"Value for theta");
  }
  
  //! \brief solve one step and update solution and time
  void solveStep() override;
  
  //! \brief gets the current solution
  [[nodiscard]] inline std::shared_ptr<const FEVector> 
  getCurrentSolution() const override { return m_u_new; }

private:

  //! \brief routine for computing delta t from CFL
  [[maybe_unused]] void computeDeltaT() const override {}
  
  //! Nonlinear solver
  std::shared_ptr<NewtonSolver>    m_nlsolver;
  //! initial condition function
  std::shared_ptr<InitialSolution> m_ic;
  //! Jac op
  std::shared_ptr<FEMatrix>        m_jac;
  //! current solution
  std::shared_ptr<FEVector>        m_u_new;
  //! old solution t^n
  std::shared_ptr<FEVector>        m_u_old;
  //! time derivative
  std::shared_ptr<FEVector>        m_u_dot;
  //! old time derivative t^n
  std::shared_ptr<FEVector>        m_u_dot_old;
  //! theta scheme parameter
  double                           m_theta;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Stepper_Classic_Theta_HPP__ */

