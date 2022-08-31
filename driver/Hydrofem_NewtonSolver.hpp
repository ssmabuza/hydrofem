// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_NewtonSolver_HPP__
#define __Hydrofem_NewtonSolver_HPP__


#include "Hydrofem_LinearSolvers.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_OptionHandler.hpp"
#include "Hydrofem_Assembler_Base.hpp"


namespace hydrofem
{
  
class NewtonSolver
  :
  public Optionable
{
public:
  
  /**  \brief  */
  NewtonSolver(const std::shared_ptr<OptionHandler>& option_handler,
               const std::shared_ptr<Assembler_Base>& nlp)
    :
    Optionable(option_handler), m_initialized(false)
  {
    // get nonlinear problem
    m_nlp = nlp;
    option_handler->parse();
    m_lob = m_nlp->lob();
  }
  
  /**  \brief  */
  virtual ~NewtonSolver() {}
  
  /**  \brief  */
  void initialize();
  
  /**  \brief  */
  void solve(const std::shared_ptr<FEVector>& U_result,
             const std::shared_ptr<const FEVector>& U_guess) const;
  
  /**  \brief  */
  double solveStep(const std::shared_ptr<FEVector>& U_result,
                   const std::shared_ptr<const FEVector>& U_guess) const;
  
private: 
  
  /** \brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config);
  
  // linear obj builder type
  using LOB = LinearObjectBuilder;
  // actual lob
  std::shared_ptr<LOB> m_lob;
  // linear solver interface
  std::shared_ptr<LinearSolverInterface> m_linear_solver;
  // nonlinear problem
  std::shared_ptr<Assembler_Base>        m_nlp;
  //! increment in NLS
  std::shared_ptr<FEVector>              m_delta_u;
  //! Jac op
  std::shared_ptr<FEMatrix>              m_jac;
  //! total residual
  std::shared_ptr<FEVector>              m_residual;
  //! current solution
  std::shared_ptr<FEVector>              m_u_new;
  //! old solution t^n
  std::shared_ptr<FEVector>              m_u_old;
  //! max number of fixed point its
  int                                    m_its;
  //! tolerance for fixed point its 
  double                                 m_tol;
  //! accept inexact solution
  bool                                   m_inexact;
  //! is initialized
  bool                                   m_initialized;
  //! 
  std::string                            m_linear_solver_name;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_NewtonSolver_HPP__ */
