// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_StepperFactory_HPP__
#define __Hydrofem_StepperFactory_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Stepper.hpp"
#include "Hydrofem_InitialSolution.hpp"

namespace hydrofem
{

using SemiDiscreteProblem = Assembler_Base;

class StepperFactory
  : public Optionable
{
public:

  StepperFactory(const std::shared_ptr<InitialSolution>& ic,
                 const std::shared_ptr<OptionHandler>& option_handler,
                 const std::shared_ptr<SemiDiscreteProblem>& semi_discrete_problem)
    :
    Optionable(option_handler)
  {
    option_handler->parse();
    m_ic = ic;
    m_option_handler = option_handler;
    m_semi_discrete_problem = semi_discrete_problem;
  }

  virtual ~StepperFactory() = default;

  std::shared_ptr<Stepper> build() const;

private:

  virtual void addOptionsCallback(po::options_description& config)
  {
    // nothing to parse yet
    config.add_options()
      ("stepper-name",po::value<std::string>(&m_name)->default_value("theta"),"Stepper name")
      ("stepper-type",po::value<std::string>(&m_type)->default_value("classic"),"Stepper type");
  }

  // name of stepper 
  std::string m_name;
  // type of stepper (classic or not classic/modern)
  std::string m_type;
  //
  std::shared_ptr<InitialSolution> m_ic;
  // parameters and options
  std::shared_ptr<OptionHandler> m_option_handler;
  // the semi-discrete problem
  std::shared_ptr<SemiDiscreteProblem> m_semi_discrete_problem;

};

} 
// end namespace hydrofem

#endif /** __Hydrofem_StepperFactory_HPP__ */