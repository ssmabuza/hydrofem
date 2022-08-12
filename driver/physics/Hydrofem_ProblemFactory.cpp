// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_ProblemFactory.hpp"
#include "Hydrofem_Problem_ReactiveFlow.hpp"
#include "Hydrofem_Problem_TaylorDispersion.hpp"
#include "Hydrofem_Problem_Bioseparation.hpp"
#include "Hydrofem_Problem_Poisson.hpp"

namespace hydrofem
{

ProblemFactory::
ProblemFactory(const std::shared_ptr<OptionHandler>& option_handler)
  :
  Optionable(option_handler)
{
  m_option_handler = option_handler;
}

std::shared_ptr<Problem>
ProblemFactory::build() const
{
  std::shared_ptr<Problem> problem;
  // main problem
  if (m_problem_name=="taylor-dispersion")
    // a simple advection-diffusion equation
    problem = std::make_shared<Problem_TaylorDispersion>(m_problem_name);
  else if (m_problem_name=="bioseparation")
    // the bioseparation problem flow from bottom to top
    problem = std::make_shared<Problem_Bioseparation>(m_problem_name,m_option_handler);
  else if (m_problem_name=="poisson")
    // simple Poisson equation with zero BC
    problem = std::make_shared<Problem_Poisson>(m_problem_name);
  else {
    throw std::logic_error("Error in ProblemFactory::build(), invalid problem name provided.");
  }
  problem->init();
  return problem;
}

}
// end namespace hydrofem
