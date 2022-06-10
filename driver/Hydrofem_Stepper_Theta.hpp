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
#include "Hydrofem_LinearSolvers.hpp"
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
    // get initial condition
    m_ic = ic;
    // get linear object builder 
    const auto lob = m_ic->getAppLOB();
    // build solution and residual vectors
    m_u_new = lob->createVector(); m_u_new->setZero();
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
  }

  //! Dtor
  ~Stepper_Theta() override = default;
  
  /** @brief options to be parsed for stepper */
  virtual void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("stepperTheta",po::value<double>(&m_theta)->default_value(0.5),"Value for theta")
      ("stepperInitialTime",po::value<double>(&m_t0)->default_value(0.0),"Initial time")
      ("stepperTimeStep",po::value<double>(&m_delta_t)->default_value(NAN),"Write the solution for output in ParaView")
      ("nonlinearSolverNumIterations",po::value<int>(&m_its)->default_value(5),"Write the solution for output in ParaView")
      ("nonlinearSolverTolerance",po::value<double>(&m_tol)->default_value(1.0e-7),"Write the solution for output in ParaView");
  }
  
  //! \brief solve one step and update solution and time
  void solveStep() override;
  
  //! \brief get the current time
  [[nodiscard]] inline double time() const override
  { return m_time; }
  
  //! \brief gets the current solution
  [[nodiscard]] inline std::shared_ptr<FEVector> getCurrentSolution() const override
  { return m_u_new; }

  //! \brief manual alteration of the fixed point iteration
  inline void setNonlinearSolverParams(const int num_its,
                                       const double tol,
                                       const bool inexact)
  {
    m_its = num_its;
    m_tol = tol;
    m_inexact = inexact;
  }

  //! \brief sets default nonlinear solver options
  [[maybe_unused]] inline void
  setSolverOptions(const int num_its = 5,
                   const double tol = 1.0e-7,
                   const bool inexact = true)
  {
    setNonlinearSolverParams(num_its,tol,inexact);
  }
  
  //! \brief routine for computing delta t from CFL
  [[maybe_unused]] void computeDeltaT(double& delta_t) {}
  
private:
  
  //! tool for assembling residual and Jacobian
  std::shared_ptr<Assembler_Base>  m_assembler;
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
  //! delta t
  double                           m_delta_t;
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

