// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_StepperFactory.hpp"

#include "Hydrofem_Stepper_Theta.hpp"
#include "Hydrofem_Stepper_ExtrapolatedEuler.hpp"

namespace hydrofem
{

std::shared_ptr<Stepper> StepperFactory::build() const
{
  std::shared_ptr<Stepper> stepper;
  if (m_name == "theta") {
    stepper = std::make_shared<Stepper_Theta>(m_ic,m_semi_discrete_problem,m_option_handler);
  } else if (m_name == "extrapolated-euler") {
    stepper = std::make_shared<Stepper_ExtrapolatedEuler>(m_ic,m_semi_discrete_problem,m_option_handler);
  } else {
    throw std::runtime_error("Error in \"StepperFactory::build()\" invalid stepper name given, stepper-name = " + m_name);
  }
  assert(stepper);
  return stepper;
}

}
// end namespace hydrofem
