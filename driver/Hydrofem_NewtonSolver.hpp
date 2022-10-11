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

  /** \brief */
  void finalize() const;
  
  /**  \brief Full solve for steady problems  */
  void solve(const std::shared_ptr<FEVector>& U_result,
             const std::shared_ptr<const FEVector>& U_guess) const;
  
  /**  \brief  */
  void solveStep(const std::shared_ptr<FEVector>& U_result,
                 const std::shared_ptr<const FEVector>& U_guess) const;

  /**  \brief  */
  void solveStep(const std::shared_ptr<FEVector>& U_result,
                 const std::shared_ptr<const FEVector>& dot_U_guess,
                 const std::shared_ptr<const FEVector>& U_guess) const;

  inline void set_dt(double dt) { m_delta_t = dt; }
  inline void set_time(double time) { m_time = time; }
  inline void set_beta(double beta) { m_beta = beta; }

  bool reachedEnd() const
  { 
    m_converged = (m_res <= m_tol);
    return ((m_num_its >= m_its) || m_converged); 
  }

private:
  
  mutable int m_num_its = 0;
  mutable double m_res;

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
  //! for transient problems
  double                                 m_delta_t = NAN;
  //! time 
  double                                 m_time = NAN;
  //! beta 
  double                                 m_beta = NAN;
  //! converged flag
  mutable bool m_converged = false;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_NewtonSolver_HPP__ */
