// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Stepper_Theta_HPP__
#define __Hydrofem_Stepper_Theta_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Stepper.hpp"
#include "Hydrofem_NewtonSolver.hpp"
#include "Hydrofem_Assembler_Base.hpp"
#include "Hydrofem_InitialSolution.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem 
{

/**
 * \brief The theta scheme 
 */
class Stepper_Theta
  :
  public Stepper
{
public:
  
  //! Ctor
  Stepper_Theta(const std::shared_ptr<InitialSolution>& ic,
                const std::shared_ptr<Assembler_Base>& assembler,
                const std::shared_ptr<OptionHandler>& option_handler) : Stepper(option_handler)
  {
    // parse the options
    option_handler->parse();
    // get initial condition
    m_ic = ic;
    // get linear object builder 
    const auto lob = m_ic->getAppLOB();
    // build solution and residual vectors
    m_u_new = m_ic->get_evaluatedField(); //lob->createVector(); m_u_new->setZero();
    m_u_old = lob->createVector(); m_u_old->setZero();
    m_delta_u = lob->createVector(); m_delta_u->setZero();
    m_u_dot = lob->createVector(); m_u_dot->setZero();
    m_u_dot_old = lob->createVector(); m_u_dot_old->setZero();
    m_residual = lob->createVector(); m_residual->setZero();
    // build Jacobian matrix
    m_jac = lob->createSparseMatrix();
    // accept inexact solution
    m_inexact = true;
    // get the assembler
    m_assembler = assembler;
    // build the nonlinear solver
    m_nlsolver = std::make_shared<NewtonSolver>(option_handler,m_assembler);
    m_nlsolver->initialize();
  }

  //! Dtor
  ~Stepper_Theta() override = default;
  
  /** @brief options to be parsed for stepper */
  virtual void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("stepper-theta",po::value<double>(&m_theta)->default_value(0.5),"Value for theta")
      ("stepper-t0",po::value<double>(&m_t0)->default_value(0.0),"Initial time")
      ("stepper-tf",po::value<double>(&m_tf)->default_value(0.0),"Final time")
      ("stepper-dt",po::value<double>(&m_delta_t)->default_value(NAN),"Time step");
  }
  
  //! \brief solve one step and update solution and time
  void solveStep() override;
  
  //! \brief get the current time
  [[nodiscard]] inline double time() const override
  { return m_time; }
  
  //! \brief gets the current solution
  [[nodiscard]] inline std::shared_ptr<FEVector> getCurrentSolution() const override
  { return m_u_new; }

  //! \brief current time step size
  [[nodiscard]] inline double dt() const override { return m_delta_t; }

  //! \brief initial time
  [[nodiscard]] inline double t0() const override { return m_t0; }

  //! \brief final time
  [[nodiscard]] inline double tf() const override { return m_tf; }

private:

  //! \brief routine for computing delta t from CFL
  [[maybe_unused]] void computeDeltaT() const override {}
  
  //! tool for assembling residual and Jacobian
  std::shared_ptr<Assembler_Base>  m_assembler;
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
  //! total residual
  std::shared_ptr<FEVector>        m_residual;
  //! increment in NLS
  std::shared_ptr<FEVector>        m_delta_u;
  //! theta scheme parameter
  double                           m_theta;
  //! initial time
  double                           m_t0;
  //! final time
  double                           m_tf;
  //! delta t
  mutable double                   m_delta_t;
  //! current time
  double                           m_time;
  //! max number of fixed point its
  int                              m_its;
  //! tolerance for fixed point its 
  double                           m_tol;
  //! accept inexact solution
  bool                             m_inexact;
  //! time step number
  int                              m_stepnum = 0;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Stepper_Theta_HPP__ */

