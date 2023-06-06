// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_NewtonSolver.hpp"
#include "Hydrofem_Assembler_Base.hpp"

namespace hydrofem
{

void NewtonSolver::
solve(const std::shared_ptr<FEVector>& U_result,
      const std::shared_ptr<const FEVector>& U_guess) const
{
  assert(m_initialized);
  //! current solution
  std::shared_ptr<FEVector> m_u_new = m_lob->createVector();
  //! old solution t^n
  std::shared_ptr<FEVector> m_u_old = m_lob->createVector();

  solveStep(m_u_new,U_guess);
  int m = 0;
  if (completed())
  {
    // assign the solution
    *U_result = *m_u_new;
    // print statement
    finalize();
    // return
    return;
  }
  
  while (!(completed()))
  {
    m++;
    *m_u_old = *m_u_new;
    solveStep(m_u_new,m_u_old);
  }
  *U_result = *m_u_new;
  finalize();
}

void NewtonSolver::
solveStep(const std::shared_ptr<FEVector>& U_result,
          const std::shared_ptr<const FEVector>& U_guess) const
{
  assert(m_initialized);
  // compute the residual
  m_nlp->buildResidualAndJacobian(U_guess,m_residual,m_jac,m_beta);
  m_nlp->finalizeAssembly(m_residual,m_jac);
  // get the norm of the residual
  m_res = m_residual->norm();
  // check if residual norm is not usual
  if (std::isnan(m_res) || std::isinf(m_res))
  {
    std::string behavior = std::isnan(m_res)? "NAN" : "INF"; 
    std::stringstream ss;
    ss << "Error in NewtonSolver::solveInitialStep!!! residual norm is " 
          + behavior + " which is equal to " << m_res << ", aborting!!";
    throw std::runtime_error(ss.str());
  }
  *m_residual *= -1.0;
  // get the first initial solution
  m_delta_u->setZero();
  // solve the system 
  m_linear_solver->solve(m_jac,m_residual,m_delta_u);
  // get the new solution
  *U_result = *m_delta_u + *U_guess;
  // increment number of iterations
  ++m_num_its;
}

void NewtonSolver::
solveStep(const std::shared_ptr<FEVector>& U_result,
          const std::shared_ptr<const FEVector>& dot_U_guess,
          const std::shared_ptr<const FEVector>& U_guess) const
{
  assert(m_initialized);
  assert(U_result);
  assert(dot_U_guess);
  assert(U_guess);
  assert(m_nlp);
  assert(m_delta_u);
  assert(m_linear_solver);

  // compute the residual
  m_nlp->buildResidualAndJacobian(U_guess,dot_U_guess,m_residual,m_jac,m_time,m_delta_t,m_beta);
  m_nlp->finalizeAssembly(m_residual,m_jac);
  // get the norm of the residual
  m_res = m_residual->norm();
  // check if residual norm is not usual
  if (std::isnan(m_res) || std::isinf(m_res))
  {
    std::stringstream ss;
    ss << "Error in NewtonSolver::solveInitialStep!!! residual norm is " << m_res << ", aborting!!";
    throw std::runtime_error(ss.str());
  }
  *m_residual *= -1.0;
  // get the first initial solution
  m_delta_u->setZero();
  // solve the system 
  m_linear_solver->solve(m_jac,m_residual,m_delta_u);
  // get the new solution
  *U_result = *m_delta_u + *U_guess;
  // increment number of iterations
  ++m_num_its;
}

void NewtonSolver::addOptionsCallback(po::options_description &config)
{
  config.add_options()
    ("newton-linear-solver", po::value<std::string>(&m_linear_solver_name)->default_value("gmres"),"linear solver name")
    ("newton-accept-inexact", po::value<bool>(&m_inexact)->default_value(true),"flag for accepting inexact solution")
    ("newton-tol", po::value<double>(&m_tol)->default_value(1.0e-6),"tolerance level")
    ("newton-max-its", po::value<int>(&m_its)->default_value(50),"max number of iterations");
}

void NewtonSolver::initialize()
{
  // build the Jacobian 
  m_jac = m_lob->createSparseMatrix();
  // build all the other objs
  m_delta_u = m_lob->createVector();
  // build residual vector
  m_residual = m_lob->createVector();
  // build the linear solver
  m_linear_solver = std::make_shared<LinearSolverInterface>(m_linear_solver_name);
  // set initialized to true
  m_initialized = true;
}

void NewtonSolver::finalize() const
{
  if (m_converged)
  {
    std::cout << std::endl;
    std::cout << "Fixed point iteration converged  *********** "<< std::endl;
    std::cout << "Fixed point iteration number of iterations = " << m_num_its << std::endl;
    std::cout << "Fixed point iteration residual norm        = " << m_res << std::endl;
    std::cout << std::endl;
    std::cout << "=========== END Fixed Point Iteration =========================" << std::endl;
    
  } else { 
    std::cout << "Fixed point iteration did not converge in " << m_its;
    if (m_inexact)
    {
      std::cout << " iterations. Accepting solution!!" << std::endl;
      std::cout << "Fixed point iteration number of iterations = " << m_num_its << std::endl;
      std::cout << "Fixed point iteration residual norm        = " << m_res << std::endl;
    } else {
      std::cout << " iterations. Aborting!!" << std::endl;
      std::cout << std::endl;
      std::cout << "=========== END Fixed Point Iteration ========================" << std::endl;
      exit(1);
    }
  }
}

}
// end namespace hydrofem
