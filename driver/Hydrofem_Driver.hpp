// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Driver_HPP__
#define __Hydrofem_Driver_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

/**
 * \brief A pure base class for the driver 
 */
class Driver
  :
  public Optionable
{
public:
  
  /** \brief Ctor */
  explicit Driver(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Optionable(option_handler)
  {
    m_option_handler = option_handler;
    option_handler->parse();
  }
  
  /** \brief Dtor */
  ~Driver() override = default;

  /** \brief The main routine that calls all solvers */
  virtual void solve() = 0;

  /** \brief The setup for the solvers */
  virtual void setup() = 0;
  
protected:
  
  /** \brief options to be parsed for solver */
  void addOptionsCallback(po::options_description &/*config*/) override
  {
    // nothing to parse
  }
  
  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Driver_HPP__ */
