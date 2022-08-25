// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Stepper_Theta.hpp"

namespace hydrofem
{

void Stepper_Theta::solveStep()
{
  // set initial time
  if (m_stepnum == 0) m_time = m_t0;
  // initialize
  if (std::fabs(m_time - m_t0) < m_delta_t/2.0)
    *m_u_new = *(m_ic->get_evaluatedField());
  // safe copy of u
  *m_u_old = *m_u_new;
  // safe copy of u_dot
  *m_u_dot_old = *m_u_dot;
  // build initial residual and Jacobian
  m_assembler->buildResidualAndJacobian(m_u_old,m_u_dot_old,m_residual,m_jac,m_time,m_delta_t,1.0/(m_theta*m_delta_t));
  *m_residual *= -1.0;
  // get the first initial solution
  m_delta_u->setZero();
  // solve the system 
  LinearSolvers::solveSystemGMRES(m_jac,m_residual,m_delta_u);
  // update u
  *m_u_new = *m_u_old + *m_delta_u;
  // update u_dot
  (*m_u_dot) = (1.0/(m_theta*m_delta_t))*(*m_u_new)-(1.0/(m_theta*m_delta_t))*(*m_u_old)-((1.0-m_theta)/m_theta)*(*m_u_dot_old);
  // residual norm
  double res (m_tol+1);
  // update the time
  m_time += m_delta_t;
  // fixed point iteration
  int m = 0;
  // force fixed-point iteration to do at least one solve for implicit stepping
  while (((res >= m_tol) && (m <= m_its)) || (m < 1))
  {
    // build the residual and Jacobian
    m_assembler->buildResidualAndJacobian(m_u_new,m_u_dot,m_residual,m_jac,m_time,m_delta_t,1.0/(m_theta*m_delta_t));
    *m_residual *= -1.0;
    // check if fixed point iteration has converged
    res = m_residual->norm();
    std::cout << "End of fixed point iteration step = " << m << ", residual norm = " << res << std::endl;
    if (std::isnan(res) || std::isinf(res))
    {
      std::stringstream ss;
      ss << "Error!!! residual norm is " << res << ", aborting!!";
      throw std::runtime_error(ss.str());
    }
    // get the first initial solution
    m_delta_u->setZero();
    // solve the system
    LinearSolvers::solveSystemGMRES(m_jac,m_residual,m_delta_u);
    // update u
    *m_u_new = *m_u_new + *m_delta_u;
    // update u_dot
    *m_u_dot = (1.0/(m_theta*m_delta_t))*(*m_u_new) - (1.0/(m_theta*m_delta_t))*(*m_u_old) - ((1.0-m_theta)/m_theta)*(*m_u_dot_old);
    // increase the fixed point iteration count
    m++;
  }
  
  // write out the convergence status
  if (m < m_its)
  {
    
    std::cout << std::endl;
    std::cout << "Fixed point iteration converged  *********** "<< std::endl;
    std::cout << "Fixed point iteration number of iterations = " << m << std::endl;
    std::cout << "Fixed point iteration residual norm        = " << m_residual->norm() << std::endl;
    std::cout << std::endl;
    std::cout << "=========== END Fixed Point Iteration =========================" << std::endl;
    
  } else {
    
    std::cout << "Fixed point iteration did not converge in " << m_its;
    if (m_inexact)
    {
      std::cout << " iterations. Accepting solution!!" << std::endl;
    } else {
      std::cout << " iterations. Aborting!!" << std::endl;
      std::cout << std::endl;
      std::cout << "=========== END Fixed Point Iteration ========================" << std::endl;
      exit(1);
    }
  }
  
  m_stepnum++;
}

}
// end namespace hydrofem

