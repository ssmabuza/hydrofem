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

void converged_statement(int,double);  
  
void NewtonSolver::
solve(const std::shared_ptr<FEVector>& U_result,
      const std::shared_ptr<const FEVector>& U_guess) const
{
  assert(m_initialized);
  double res = solveStep(m_u_new,U_guess);
  int m = 0;
  if (res < m_tol)
  {
    // assign the solution
    *U_result = *m_u_new;
    // print statement
    converged_statement(m,res);
    // return
    return;
  }
  
  while ((res >= m_tol) && (m <= m_its))
  {
    m++;
    *m_u_old = *m_u_new;
    res = solveStep(m_u_new,m_u_old);
  }
  *U_result = *m_u_new;
  converged_statement(m,res);
}

double NewtonSolver::
solveStep(const std::shared_ptr<FEVector>& U_result,
          const std::shared_ptr<const FEVector>& U_guess) const
{
  assert(m_initialized);
  // compute the residual
  m_nlp->buildResidualAndJacobian(U_guess,m_residual,m_jac,0.0);
  m_nlp->finalizeAssembly(m_residual,m_jac);
  // get the norm of the residual
  const auto res = m_residual->norm();
  // check if residual norm is not usual
  if (std::isnan(res) || std::isinf(res))
  {
    std::stringstream ss;
    ss << "Error in NewtonSolver::solveInitialStep!!! residual norm is " << res << ", aborting!!";
    throw std::runtime_error(ss.str());
  }
  *m_residual *= -1.0;
  // get the first initial solution
  m_delta_u->setZero();
  // solve the system 
  m_linear_solver->solve(m_jac,m_residual,m_delta_u);
  // get the new solution
  *U_result = *m_delta_u + *U_guess;
  // return the residual 
  return res;
}

void converged_statement(int m, double res)
{
  std::cout << std::endl;
  std::cout << "************  Fixed point iteration converged  *********** "<< std::endl;
  std::cout << "Fixed point iteration number of iterations = " << m << std::endl;
  std::cout << "Fixed point iteration residual norm        = " << res << std::endl;
  std::cout << std::endl;
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
  m_residual = m_lob->createVector();
  m_u_new = m_lob->createVector();
  m_u_old = m_lob->createVector();
  // build the linear solver
  m_linear_solver = std::make_shared<LinearSolverInterface>(m_linear_solver_name);
  // set initialized to true
  m_initialized = true;
}

}
// end namespace hydrofem
