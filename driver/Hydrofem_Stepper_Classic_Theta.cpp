// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Stepper_Classic_Theta.hpp"

namespace hydrofem
{

void Stepper_Classic_Theta::solveStep()
{
  // update the time
  m_time += m_delta_t;
  // set initial time
  if (m_stepnum == 0)
  {
    // set time to initial time
    m_time = m_t0;
    // updates solver delta t 
    m_nlsolver->set_dt(m_delta_t);
    // updates solve time
    m_nlsolver->set_time(m_time);
    // update solver beta
    m_nlsolver->set_beta(1.0/(m_theta*m_delta_t));
  }

  // safe copy of u
  *m_u_old = *m_u_new;
  // safe copy of u_dot
  *m_u_dot_old = *m_u_dot;
  // update u_dot
  (*m_u_dot) = (1.0/(m_theta*m_delta_t))*(*m_u_new)-(1.0/(m_theta*m_delta_t))*(*m_u_old)-((1.0-m_theta)/m_theta)*(*m_u_dot_old);
  // solve the problem
  m_nlsolver->solveStep(m_u_new,m_u_dot_old,m_u_old);
  // compute delta t
  computeDeltaT();
  // residual norm
  //double res (m_tol+1);
  // updates solver delta t 
  m_nlsolver->set_dt(m_delta_t);
  // updates solve time
  m_nlsolver->set_time(m_time);
  // update solver beta
  m_nlsolver->set_beta(1.0/(m_theta*m_delta_t));
  // fixed point iteration
  int m = 0;
  // force fixed-point iteration to do at least one solve for implicit stepping
  while (!(m_nlsolver->completed()))
  {
    // solve the problem
    m_nlsolver->solveStep(m_u_new,m_u_dot,m_u_new);
    // compute delta t
    computeDeltaT();
    // update u_dot 
    //TODO: debug, try on the heat equation
    *m_u_dot = (1.0/(m_theta*m_delta_t))*(*m_u_new) - (1.0/(m_theta*m_delta_t))*(*m_u_old) - ((1.0-m_theta)/m_theta)*(*m_u_dot_old);
    // increase the fixed point iteration count
    m++;
  }
  
  // write out the convergence status
  m_nlsolver->finalize();
  // increase time step count  
  m_stepnum++;
  
}

}
// end namespace hydrofem

 