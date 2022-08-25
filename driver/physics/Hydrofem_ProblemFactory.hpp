// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_ProblemFactory_HPP__
#define __Hydrofem_ProblemFactory_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

/**
 * \brief The class that stores all new problems and contructs a solver at run time
 */
class ProblemFactory
  :
  public Optionable
{
public:
  
  /** \brief Ctor */
  explicit ProblemFactory(const std::shared_ptr<OptionHandler>& option_handler);
  
  /** \brief Dtor */
  virtual ~ProblemFactory() = default;
  
  /** \brief Solver builder  */
  std::shared_ptr<Problem> build() const;

private:

  /** \brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config)
  {
    // nothing to parse
    config.add_options()
      ("problem-name",po::value<std::string>(&m_problem_name)->default_value("poisson"),
       "Continuous problem name.");
  }
  
  std::string m_problem_name;
  
  std::shared_ptr<OptionHandler> m_option_handler;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_ProblemFactory_HPP__ */
