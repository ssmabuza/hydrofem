// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Stepper_ExtrapolatedEuler_HPP__
#define __Hydrofem_Stepper_ExtrapolatedEuler_HPP__

#include "Hydrofem_Stepper.hpp"

#include "Hydrofem_NewtonSolver.hpp"
#include "Hydrofem_Assembler_Base.hpp"
#include "Hydrofem_InitialSolution.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem
{

class Stepper_ExtrapolatedEuler
  :
  public Stepper
{
public:

  /** \brief Ctor */
  Stepper_ExtrapolatedEuler(const std::shared_ptr<InitialSolution>& ic,
                            const std::shared_ptr<Assembler_Base>& assembler,
                            const std::shared_ptr<OptionHandler>& option_handler) : Stepper(option_handler,assembler)
  {
    // parse options
    option_handler->parse();
    // set the initial condition
    m_ic = ic;
    // get linear object builder 
    const auto lob = m_ic->getAppLOB(); assert(lob);
    // build the nonlinear solver
    m_nlsolver = std::make_shared<NewtonSolver>(option_handler,m_assembler);
    // initialize nonlinear solver
    m_nlsolver->initialize(); assert(m_nlsolver);
    // build u new
    m_u_new = lob->createVector();
    // build u half
    m_u_half = lob->createVector();
    // build u old
    m_u_old = lob->createVector();
    // build u dot
    m_u_dot = lob->createVector();

  }

  /** \brief Dtor */
  virtual ~Stepper_ExtrapolatedEuler() {}

  /** \brief options to be parsed for stepper */
  virtual void addOptionsCallback(po::options_description& /*config*/)
  {
    // nothing to parse
  }  

  /** \brief Main solve routine */
  void solveStep() override;

  std::shared_ptr<const FEVector> getCurrentSolution() const override
  { return m_u_new; }

private:
  
  //! Initial condition
  std::shared_ptr<InitialSolution> m_ic;
  //! Nonlinear solver
  std::shared_ptr<NewtonSolver> m_nlsolver;
  //! solution at t^n+1
  std::shared_ptr<FEVector> m_u_new;
  //! solution at t^n+0.5
  std::shared_ptr<FEVector> m_u_half;
  //! solution at t^n
  std::shared_ptr<FEVector> m_u_old;
  //! time derivative
  std::shared_ptr<FEVector> m_u_dot;

};

}
// end namespace hydrofem 

#endif /** __Hydrofem_Stepper_ExtrapolatedEuler_HPP__ */