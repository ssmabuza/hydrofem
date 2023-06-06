// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Stepper_ExtrapolatedEuler.hpp"

namespace hydrofem
{

void Stepper_ExtrapolatedEuler::solveStep() 
{
  /// step 1: solve backward Euler half step problem
  // update the time
  m_time += m_delta_t/2.0;
  // set initial time
  if (m_stepnum == 0)
  {
    // set time to initial time
    m_time = m_t0;
    // updates solver delta t 
    m_nlsolver->set_dt(m_delta_t/2.0);
    // updates solve time
    m_nlsolver->set_time(m_time);
    // update solver beta
    m_nlsolver->set_beta(2.0/m_delta_t);
  }

  // safe copy of u
  *m_u_old = *m_u_new;
  // update u_dot
  (*m_u_dot) = (2.0/m_delta_t)*(*m_u_half-*m_u_old);
  // solve the problem
  m_nlsolver->solveStep(m_u_half,m_u_dot,m_u_old);
  // compute delta t
  computeDeltaT();
  // updates solver delta t 
  m_nlsolver->set_dt(m_delta_t/2.0);
  // updates solve time
  m_nlsolver->set_time(m_time);
  // update solver beta
  m_nlsolver->set_beta(2.0/m_delta_t);
  // fixed point iteration
  // force fixed-point iteration to do at least one solve for implicit stepping
  while (!(m_nlsolver->completed()))
  {
    // solve the problem
    m_nlsolver->solveStep(m_u_half,m_u_dot,m_u_half);
    // compute delta t
    computeDeltaT();
    // update u_dot 
    //TODO: debug, try on the heat equation
    *m_u_dot = (2.0/m_delta_t)*(*m_u_half-*m_u_old);
  }
  
  // write out the convergence status
  m_nlsolver->finalize();
  // increase time step count  
  m_stepnum++;

  /// step 2: update solution
  m_time += m_delta_t/2.0;
  *m_u_new = 2.0 * (*m_u_half) - *m_u_old;
}

}
// end namespace hydrofem



